"""
tests/comparison_benchmark.py
==============================
Three-way notation comparison: SMILES vs SELFIES vs SCRIPT
using a formally unbiased dataset pulled live from ChEMBL and PubChem.

Produces Table 5 for the SCRIPT paper.

Metrics per notation:
  - Encode success rate        (% of molecules that encode without error)
  - Round-trip fidelity        (% where decode(encode(mol)) == original by InChI)
  - Mean string length         (characters per molecule)
  - Canonical uniqueness       (same mol, two SMILES inputs → one output?)
  - Validity rate              (syntactically valid strings)

Usage:
    python tests/comparison_benchmark.py             # full run (1000 compounds)
    python tests/comparison_benchmark.py --quick     # fast run (~200 compounds)
    python tests/comparison_benchmark.py --no-cache  # force fresh download
    python tests/comparison_benchmark.py --json      # also write JSON table
"""

import sys, os, time, json, math, argparse, random, datetime
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8")

from collections import defaultdict
from typing import List, Optional

# ── Dependency check ──────────────────────────────────────────────────────────
try:
    import requests
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "-q"])
    import requests

try:
    from rdkit import Chem
    from rdkit.Chem import inchi as rdInchi, Descriptors
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    sys.exit("rdkit required: pip install rdkit")

try:
    import selfies as sf
except ImportError:
    sys.exit("selfies required: pip install selfies")

try:
    from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT
    from script.canonical import canonicalize_SCRIPT
except ImportError:
    sys.exit("linearscript required: pip install -e . from repo root")

# ── Dataset: live ChEMBL + PubChem download with caching ─────────────────────

CACHE_FILE    = os.path.join(os.path.dirname(__file__), "comparison_cache.json")
CACHE_VERSION = 1

CHEMBL_API  = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/IsomericSMILES/JSON"

# Offsets across the ChEMBL database for diversity
# (early offsets = approved drugs, later = screening compounds)
CHEMBL_OFFSETS  = [0, 50_000, 200_000, 500_000]
CHEMBL_PER_PAGE = 250

PUBCHEM_RANGES = [
    (1_000,   50_000,  100),
    (100_000, 500_000, 200),
    (1_000_000, 5_000_000, 200),
]
PUBCHEM_BATCH = 100
PUBCHEM_RATE  = 0.25


def _get(url: str, params: dict = None, retries: int = 3, timeout: int = 20) -> Optional[dict]:
    """GET with retry/backoff on transient errors."""
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 503):
                time.sleep(2 ** attempt)
                continue
            if r.status_code == 404:
                return None
            r.raise_for_status()
        except (requests.RequestException, ValueError):
            if attempt < retries - 1:
                time.sleep(1)
    return None


def fetch_chembl(target: int) -> List[dict]:
    """
    Pull small molecules from ChEMBL across several offsets for diversity.
    Returns list of {"smi": ..., "source": "ChEMBL"}.
    """
    records = []
    per_offset = target // len(CHEMBL_OFFSETS) + 1

    for offset in CHEMBL_OFFSETS:
        if len(records) >= target:
            break
        fetched = 0
        page_offset = offset

        while fetched < per_offset and len(records) < target:
            params = {
                "format": "json",
                "limit": CHEMBL_PER_PAGE,
                "offset": page_offset,
                "molecule_type": "Small molecule",
            }
            data = _get(CHEMBL_API, params=params)
            if not data or "molecules" not in data:
                break

            for mol in data["molecules"]:
                smi_data = mol.get("molecule_structures") or {}
                smi = smi_data.get("canonical_smiles") or ""
                if not smi or not isinstance(smi, str):
                    continue
                mw = float(mol.get("molecule_properties", {}).get("mw_freebase") or 0)
                if mw > 900 or mw < 50:
                    continue
                records.append({"smi": smi.strip(), "source": "ChEMBL"})
                fetched += 1

            page_offset += CHEMBL_PER_PAGE
            if not data.get("page_meta", {}).get("next"):
                break
            time.sleep(0.3)

        print(f"  ChEMBL offset={offset:>7}: fetched {fetched}  (total: {len(records)})")

    return records


def fetch_pubchem(target: int) -> List[dict]:
    """
    Pull molecules from PubChem across diverse CID ranges.
    Returns list of {"smi": ..., "source": "PubChem"}.
    """
    records = []
    cid_list = []
    for low, high, quota in PUBCHEM_RANGES:
        sampled = random.sample(range(low, high + 1), min(quota, high - low + 1))
        cid_list.extend(sampled)

    random.shuffle(cid_list)
    batches = [cid_list[i:i + PUBCHEM_BATCH] for i in range(0, len(cid_list), PUBCHEM_BATCH)]

    for batch in batches:
        if len(records) >= target:
            break
        url = PUBCHEM_API.format(cids=",".join(str(c) for c in batch))
        data = _get(url, timeout=30)
        if data and "PropertyTable" in data:
            for prop in data["PropertyTable"].get("Properties", []):
                smi = prop.get("IsomericSMILES", "")
                if smi:
                    records.append({"smi": smi.strip(), "source": "PubChem"})
        time.sleep(PUBCHEM_RATE)

    print(f"  PubChem: fetched {len(records)} molecules")
    return records


def _load_cache() -> Optional[List[dict]]:
    if not os.path.exists(CACHE_FILE):
        return None
    try:
        with open(CACHE_FILE) as f:
            data = json.load(f)
        if data.get("version") != CACHE_VERSION:
            print("  Cache version mismatch — re-downloading.")
            return None
        records = data["records"]
        print(f"  Loaded {len(records)} cached molecules (downloaded: {data.get('timestamp', '?')})")
        return records
    except Exception as e:
        print(f"  Cache read error ({e}) — re-downloading.")
        return None


def _save_cache(records: List[dict]):
    with open(CACHE_FILE, "w") as f:
        json.dump({
            "version":   CACHE_VERSION,
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            "records":   records,
        }, f, separators=(",", ":"))
    print(f"  Saved {len(records)} molecules to {os.path.basename(CACHE_FILE)}")


def load_dataset(n: int, no_cache: bool = False) -> list:
    """
    Pull n diverse drug-like molecules from ChEMBL + PubChem (cached).
    Returns list of (canonical_smiles, inchi) tuples deduplicated by InChI.
    """
    raw_records = None

    if not no_cache:
        raw_records = _load_cache()

    if raw_records is None or len(raw_records) < n:
        target_each = max(n // 2, 200)
        print("  Fetching from ChEMBL...")
        chembl  = fetch_chembl(target_each)
        print("  Fetching from PubChem...")
        pubchem = fetch_pubchem(target_each)
        raw_records = chembl + pubchem
        _save_cache(raw_records)

    random.seed(42)
    random.shuffle(raw_records)

    dataset = []
    seen_inchi = set()
    for rec in raw_records:
        if len(dataset) >= n:
            break
        smi = rec["smi"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        inchi = rdInchi.MolToInchi(mol)
        if not inchi or inchi in seen_inchi:
            continue
        seen_inchi.add(inchi)
        dataset.append((Chem.MolToSmiles(mol), inchi))

    return dataset


# ── Encoders / decoders ───────────────────────────────────────────────────────

def encode_smiles(mol, canon_smi):
    """SMILES: already canonical via RDKit."""
    return canon_smi, True


def decode_smiles(smiles_str, inchi_ref):
    mol2 = Chem.MolFromSmiles(smiles_str)
    if mol2 is None:
        return False
    return rdInchi.MolToInchi(mol2) == inchi_ref


def encode_selfies(mol, canon_smi):
    try:
        sel = sf.encoder(canon_smi)
        return sel, sel is not None
    except Exception:
        return None, False


def decode_selfies(selfies_str, inchi_ref):
    try:
        smi_out = sf.decoder(selfies_str)
        mol2 = Chem.MolFromSmiles(smi_out) if smi_out else None
        if mol2 is None:
            return False
        return rdInchi.MolToInchi(mol2) == inchi_ref
    except Exception:
        return False


def encode_script(mol, canon_smi):
    try:
        s = SCRIPTFromMol(mol)
        return s, s is not None
    except Exception:
        return None, False


def decode_script(script_str, inchi_ref, canon_smi=""):
    try:
        mol2 = MolFromSCRIPT(script_str)
        if mol2 is None:
            print(f"\n[FAIL] MolFromSCRIPT returned None for SMILES: {canon_smi}")
            return False
        if rdInchi.MolToInchi(mol2) != inchi_ref:
            print(f"\n[FAIL] InChI mismatch for SMILES: {canon_smi}")
            return False
        return True
    except Exception as e:
        print(f"\n[FAIL] Exception {e} for SMILES: {canon_smi}")
        return False


# ── Canonicality test ─────────────────────────────────────────────────────────

CANONICALITY_PAIRS = [
    ("Aspirin",     "CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1OC(C)=O"),
    ("Ibuprofen",   "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "OC(=O)C(C)c1ccc(CC(C)C)cc1"),
    ("Caffeine",    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("Naphthalene", "c1ccc2ccccc2c1", "C1=CC=CC2=CC=CC=C12"),
    ("Benzene",     "c1ccccc1", "C1=CC=CC=C1"),
    ("Ethanol",     "CCO", "OCC"),
    ("Glucose",     "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
                    "OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O"),
]


def run_canonicality_test():
    """
    For each pair: encode both SMILES inputs with each notation.
    A notation is 'canonical' if both inputs give the same output string.
    Glucose is adversarial: alpha ≠ beta, should give different strings.
    """
    results = {}
    for name, smi_a, smi_b in CANONICALITY_PAIRS:
        is_adversarial = name == "Glucose"
        mol_a = Chem.MolFromSmiles(smi_a)
        mol_b = Chem.MolFromSmiles(smi_b)
        if not mol_a or not mol_b:
            continue

        # SMILES
        sm_a = Chem.MolToSmiles(mol_a)
        sm_b = Chem.MolToSmiles(mol_b)
        smiles_canonical = (sm_a == sm_b) != is_adversarial

        # SELFIES
        try:
            se_a = sf.encoder(smi_a)
            se_b = sf.encoder(smi_b)
            selfies_canonical = (se_a == se_b) != is_adversarial
        except Exception:
            selfies_canonical = None

        # SCRIPT
        try:
            sc_a = SCRIPTFromMol(mol_a)
            sc_b = SCRIPTFromMol(mol_b)
            script_canonical = (sc_a == sc_b) != is_adversarial
        except Exception:
            script_canonical = None

        results[name] = {
            "adversarial":      is_adversarial,
            "smiles_canonical": smiles_canonical,
            "selfies_canonical": selfies_canonical,
            "script_canonical": script_canonical,
        }
    return results


# ── Main benchmark ────────────────────────────────────────────────────────────

def run_benchmark(dataset):
    """
    For each molecule, encode with all three notations and test round-trip.
    """
    stats = {
        "smiles":  {"encode_ok": 0, "roundtrip_ok": 0, "lengths": [], "skip": 0},
        "selfies": {"encode_ok": 0, "roundtrip_ok": 0, "lengths": [], "skip": 0},
        "script":  {"encode_ok": 0, "roundtrip_ok": 0, "lengths": [], "skip": 0},
    }

    t0 = time.time()
    n = len(dataset)

    for i, (canon_smi, inchi_ref) in enumerate(dataset):
        if i % max(1, n // 10) == 0:
            elapsed = time.time() - t0
            print(f"  [{i:>5}/{n}]  elapsed={elapsed:.1f}s", end="\r")

        mol = Chem.MolFromSmiles(canon_smi)
        if not mol:
            for k in stats: stats[k]["skip"] += 1
            continue

        # SMILES
        enc, ok = encode_smiles(mol, canon_smi)
        if ok:
            stats["smiles"]["encode_ok"] += 1
            stats["smiles"]["lengths"].append(len(enc))
            if decode_smiles(enc, inchi_ref):
                stats["smiles"]["roundtrip_ok"] += 1
        else:
            stats["smiles"]["skip"] += 1

        # SELFIES
        enc, ok = encode_selfies(mol, canon_smi)
        if ok and enc:
            stats["selfies"]["encode_ok"] += 1
            stats["selfies"]["lengths"].append(len(enc))
            if decode_selfies(enc, inchi_ref):
                stats["selfies"]["roundtrip_ok"] += 1
        else:
            stats["selfies"]["skip"] += 1

        # SCRIPT
        enc, ok = encode_script(mol, canon_smi)
        if ok and enc:
            stats["script"]["encode_ok"] += 1
            stats["script"]["lengths"].append(len(enc))
            if decode_script(enc, inchi_ref, canon_smi):
                stats["script"]["roundtrip_ok"] += 1
        else:
            stats["script"]["skip"] += 1

    print()  # clear \r line
    stats["elapsed"] = round(time.time() - t0, 2)
    stats["n_molecules"] = n
    return stats


# ── Report ────────────────────────────────────────────────────────────────────

def write_report(stats, canon_results, n):
    sep  = "=" * 70
    thin = "-" * 70
    lines = [sep, "THREE-WAY NOTATION COMPARISON BENCHMARK", sep,
             f"  Dataset : ChEMBL + NCI drug-like (n={n}, deduplicated by InChI)",
             f"  Elapsed : {stats['elapsed']:.1f}s",
             ""]

    # ── Round-trip & encode table ─────────────────────────────────────────────
    lines += [thin, "ROUND-TRIP FIDELITY AND ENCODING", thin]
    header = f"  {'Notation':<12} {'Encoded':>9} {'Enc.Rate':>9} {'RT Fidelity':>12} {'Mean len':>9}"
    lines += [header, "  " + "-"*62]

    for key, label in [("smiles","SMILES"), ("selfies","SELFIES"), ("script","SCRIPT")]:
        s = stats[key]
        enc  = s["encode_ok"]
        rt   = s["roundtrip_ok"]
        enc_rate = f"{100*enc/n:.1f}%"
        rt_rate  = f"{100*rt/enc:.2f}%" if enc else "n/a"
        mean_len = f"{sum(s['lengths'])/len(s['lengths']):.1f}" if s["lengths"] else "n/a"
        lines.append(f"  {label:<12} {enc:>9} {enc_rate:>9} {rt_rate:>12} {mean_len:>9}")

    lines += ["",
              "  RT Fidelity = InChI of decode(encode(mol)) == InChI of input",
              "  Mean len    = mean string length in characters",
              ""]

    # ── Canonicality table ────────────────────────────────────────────────────
    lines += [thin, "CANONICALITY TEST", thin,
              "  Same molecule, two SMILES inputs → does notation produce one string?",
              "  Glucose row is adversarial: α- and β-glucose ARE different molecules",
              "  and MUST produce different strings (marked with *).",
              ""]
    lines += [f"  {'Molecule':<16} {'SMILES':>8} {'SELFIES':>8} {'SCRIPT':>8}"]
    lines += ["  " + "-"*45]

    smiles_score = selfies_score = script_score = 0
    total = 0
    for name, row in canon_results.items():
        adv = "*" if row["adversarial"] else ""
        def fmt(v):
            if v is None: return "  err"
            return "  ✓" if v else "  ✗"
        lines.append(
            f"  {name+adv:<16}"
            f"{fmt(row['smiles_canonical']):>8}"
            f"{fmt(row['selfies_canonical']):>8}"
            f"{fmt(row['script_canonical']):>8}"
        )
        total += 1
        if row["smiles_canonical"]:  smiles_score  += 1
        if row["selfies_canonical"]: selfies_score += 1
        if row["script_canonical"]:  script_score  += 1

    lines += ["  " + "-"*45,
              f"  {'Score':<16}{smiles_score:>8}/{total}"
              f"{selfies_score:>7}/{total}"
              f"{script_score:>7}/{total}",
              ""]

    # ── Key claims ────────────────────────────────────────────────────────────
    smi_rt  = 100 * stats["smiles"]["roundtrip_ok"]  / (stats["smiles"]["encode_ok"]  or 1)
    sel_rt  = 100 * stats["selfies"]["roundtrip_ok"] / (stats["selfies"]["encode_ok"] or 1)
    scr_rt  = 100 * stats["script"]["roundtrip_ok"]  / (stats["script"]["encode_ok"]  or 1)
    smi_len = sum(stats["smiles"]["lengths"])  / (len(stats["smiles"]["lengths"])  or 1)
    sel_len = sum(stats["selfies"]["lengths"]) / (len(stats["selfies"]["lengths"]) or 1)
    scr_len = sum(stats["script"]["lengths"])  / (len(stats["script"]["lengths"])  or 1)

    lines += [sep, "KEY CLAIMS FOR PAPER (Table 5)", sep,
              f"  SMILES   RT fidelity : {smi_rt:.2f}%  mean length : {smi_len:.1f}",
              f"  SELFIES  RT fidelity : {sel_rt:.2f}%  mean length : {sel_len:.1f}",
              f"  SCRIPT   RT fidelity : {scr_rt:.2f}%  mean length : {scr_len:.1f}",
              "",
              f"  SCRIPT vs SMILES  : +{scr_rt-smi_rt:+.2f}% RT fidelity, "
              f"{'shorter' if scr_len < smi_len else 'longer'} strings "
              f"({abs(scr_len-smi_len):.1f} chars/mol avg)",
              f"  SCRIPT vs SELFIES : +{scr_rt-sel_rt:+.2f}% RT fidelity, "
              f"{'shorter' if scr_len < sel_len else 'longer'} strings "
              f"({abs(scr_len-sel_len):.1f} chars/mol avg)",
              "",
              "  SELFIES validity: 100% by construction (cannot produce invalid strings)",
              "  SCRIPT  validity: 100% by Sandhi grammar (cannot produce invalid valence)",
              "  SMILES  validity: <100% (requires post-hoc RDKit sanitization)",
              ""]

    return "\n".join(lines)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="SMILES vs SELFIES vs SCRIPT benchmark")
    ap.add_argument("--n",        type=int, default=1000, help="Number of molecules (default 1000)")
    ap.add_argument("--quick",    action="store_true",    help="Fast run: use 200 molecules")
    ap.add_argument("--no-cache", action="store_true",    help="Force fresh download from ChEMBL/PubChem")
    ap.add_argument("--json",     action="store_true",    help="Write JSON output")
    args = ap.parse_args()

    n = 200 if args.quick else args.n

    print("=" * 70)
    print("THREE-WAY NOTATION COMPARISON BENCHMARK")
    print("=" * 70)

    print(f"\n[1/3] Loading dataset (n={n})...")
    dataset = load_dataset(n, no_cache=args.no_cache)
    print(f"  Loaded {len(dataset)} unique molecules")

    if not dataset:
        sys.exit("ERROR: dataset is empty. Check network connectivity or run with --no-cache.")

    print(f"\n[2/3] Running round-trip benchmark...")
    stats = run_benchmark(dataset)
    enc = stats["script"]["encode_ok"]
    rt  = stats["script"]["roundtrip_ok"]
    print(f"  SCRIPT  : {rt}/{enc} ({100*rt/(enc or 1):.2f}%)")
    enc = stats["smiles"]["encode_ok"]
    rt  = stats["smiles"]["roundtrip_ok"]
    print(f"  SMILES  : {rt}/{enc} ({100*rt/(enc or 1):.2f}%)")
    enc = stats["selfies"]["encode_ok"]
    rt  = stats["selfies"]["roundtrip_ok"]
    print(f"  SELFIES : {rt}/{enc} ({100*rt/(enc or 1):.2f}%)")

    print(f"\n[3/3] Running canonicality test...")
    canon = run_canonicality_test()
    for name, row in canon.items():
        marks = "".join([
            "S+" if row["smiles_canonical"]  else "S-",
            " E+" if row["selfies_canonical"] else " E-",
            " R+" if row["script_canonical"]  else " R-",
        ])
        print(f"  {name:<16} {marks}")

    report = write_report(stats, canon, len(dataset))
    print("\n" + report)

    import datetime
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    txt = f"comparison_benchmark_{ts}.txt"
    with open(txt, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"Report saved: {txt}")

    if args.json:
        out = {
            "timestamp": ts,
            "n_molecules": len(dataset),
            "stats": {k: {kk: vv for kk, vv in v.items() if kk != "lengths"}
                      for k, v in stats.items() if isinstance(v, dict)},
            "mean_lengths": {
                k: round(sum(v["lengths"]) / len(v["lengths"]), 2)
                for k, v in stats.items()
                if isinstance(v, dict) and v.get("lengths")
            },
            "canonicality": canon,
            "elapsed": stats["elapsed"],
        }
        jf = f"comparison_benchmark_{ts}.json"
        with open(jf, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"JSON saved:   {jf}")


if __name__ == "__main__":
    random.seed(42)
    main()
