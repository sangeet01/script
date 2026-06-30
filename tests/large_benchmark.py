"""
SCRIPT Large-Scale Benchmark
==============================
Pulls 5000+ diverse compounds from ChEMBL and PubChem automatically,
runs full SMILES -> SCRIPT -> InChI round-trip benchmark,
and writes a paper-quality report.

Usage:
    python large_benchmark.py               # full run (~5000 compounds)
    python large_benchmark.py --quick       # quick run (~500 compounds, no download wait)
    python large_benchmark.py --no-cache    # force re-download even if cache exists

Requirements:
    pip install rdkit lark requests

Output:
    large_benchmark_results.txt   human-readable report
    large_benchmark_results.json  machine-readable data for paper tables
    benchmark_cache.json          downloaded SMILES cache (reused on next run)
"""

import os
import sys
import json
import time
import random
import argparse
import datetime
from collections import defaultdict
from typing import List, Dict, Optional, Tuple

# ── Dependency check ─────────────────────────────────────────────────────────

try:
    import requests
except ImportError:
    print("Installing requests...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "-q"])
    import requests

try:
    from rdkit import Chem
    from rdkit.Chem import inchi as rdInchi
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    sys.exit("ERROR: rdkit not found. Install with: pip install rdkit")

try:
    from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT
except ImportError:
    # Try adding parent directory to path
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT
    except ImportError:
        sys.exit("ERROR: SCRIPT not found. Run from the repo root, or: pip install -e .")


# ── Constants ─────────────────────────────────────────────────────────────────

CACHE_FILE   = "benchmark_cache.json"
CACHE_VERSION = 3                    # bump to invalidate old caches

CHEMBL_API   = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PUBCHEM_API  = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/IsomericSMILES/JSON"

# How many to pull from each source
TARGET_CHEMBL  = 2500
TARGET_PUBCHEM = 2500

# ChEMBL offsets to sample across the full database for diversity
# Database has ~2.4M compounds; offset at 0, 50k, 200k, 500k, 1M gives
# different chemical classes (early = approved drugs, later = screening cmpds)
CHEMBL_OFFSETS = [0, 50_000, 200_000, 500_000, 1_000_000]
CHEMBL_PER_PAGE = 500                # records per API call

# PubChem CID ranges — sample across the full CID space for diversity
# Low CIDs = simple molecules/solvents; high CIDs = complex/drug-like
PUBCHEM_RANGES = [
    (2,       500,      50),    # very simple: amino acids, sugars, solvents
    (1_000,   5_000,   150),    # simple organics: building blocks
    (10_000,  50_000,  300),    # drug fragments
    (100_000, 500_000, 500),    # drug-like leads
    (1_000_000, 5_000_000, 800),# complex drug-like / bioactives
    (5_000_000, 50_000_000, 700),# diverse screening compounds
]
PUBCHEM_BATCH = 100                  # CIDs per API call (PubChem limit: 200)
PUBCHEM_RATE  = 0.25                 # seconds between PubChem requests


# ── Download helpers ──────────────────────────────────────────────────────────

def _get(url: str, params: dict = None, retries: int = 3, timeout: int = 20) -> Optional[dict]:
    """GET with retry on transient errors."""
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 503):          # rate-limited
                time.sleep(2 ** attempt)
                continue
            if r.status_code == 404:
                return None                           # not found, not an error
            r.raise_for_status()
        except (requests.RequestException, ValueError):
            if attempt < retries - 1:
                time.sleep(1)
    return None


def fetch_chembl(target: int) -> List[Dict]:
    """
    Pull small molecules from ChEMBL across several offsets.
    Returns list of {"smi": ..., "source": "ChEMBL", "name": ...}
    """
    records = []
    per_offset = target // len(CHEMBL_OFFSETS) + 1

    for offset in CHEMBL_OFFSETS:
        if len(records) >= target:
            break
        fetched_this_offset = 0
        page_offset = offset

        while fetched_this_offset < per_offset and len(records) < target:
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
                smi = smi_data.get("canonical_smiles") or smi_data.get("molfile")
                if not smi or not isinstance(smi, str):
                    continue
                mw = float(mol.get("molecule_properties", {}).get("mw_freebase") or 0)
                if mw > 900 or mw < 50:              # skip extreme sizes
                    continue
                records.append({
                    "smi":    smi.strip(),
                    "source": "ChEMBL",
                    "name":   mol.get("pref_name") or mol.get("chembl_id") or "",
                })
                fetched_this_offset += 1

            page_offset += CHEMBL_PER_PAGE
            if not data.get("page_meta", {}).get("next"):
                break
            time.sleep(0.3)

        print(f"  ChEMBL offset={offset:>8}: {fetched_this_offset} molecules  "
              f"(total so far: {len(records)})")

    return records


def fetch_pubchem(target: int) -> List[Dict]:
    """
    Pull molecules from PubChem across diverse CID ranges.
    Returns list of {"smi": ..., "source": "PubChem", "name": ""}
    """
    records = []

    # Build CID sample list respecting per-range quotas
    cid_list = []
    for low, high, quota in PUBCHEM_RANGES:
        sampled = random.sample(range(low, high + 1), min(quota, high - low + 1))
        cid_list.extend(sampled)
        if len(cid_list) >= target * 1.2:     # oversample, some CIDs won't exist
            break

    random.shuffle(cid_list)

    # Batch-fetch
    batches = [cid_list[i:i + PUBCHEM_BATCH] for i in range(0, len(cid_list), PUBCHEM_BATCH)]
    batches_needed = (target // PUBCHEM_BATCH) + 2

    for batch in batches[:batches_needed]:
        if len(records) >= target:
            break
        url = PUBCHEM_API.format(cids=",".join(str(c) for c in batch))
        data = _get(url, timeout=30)
        if data and "PropertyTable" in data:
            for prop in data["PropertyTable"].get("Properties", []):
                smi = prop.get("IsomericSMILES", "")
                if smi:
                    records.append({
                        "smi":    smi.strip(),
                        "source": "PubChem",
                        "name":   str(prop.get("CID", "")),
                    })
        time.sleep(PUBCHEM_RATE)

    print(f"  PubChem: {len(records)} molecules from {min(len(batches), batches_needed)} batches")
    return records


def load_rdkit_local() -> List[Dict]:
    """Load the RDKit-bundled datasets as a network-free fallback."""
    records = []
    paths = [
        ("/usr/local/lib/python3.12/dist-packages/rdkit/Contrib/FreeWilson/data/CHEMBL2321810.smi", "RDKit/ChEMBL"),
        ("/usr/local/lib/python3.12/dist-packages/rdkit/Data/NCI/first_5K.smi", "RDKit/NCI"),
        ("/usr/local/lib/python3.12/dist-packages/rdkit/Contrib/Fastcluster/cdk2.smi", "RDKit/CDK2"),
        ("/usr/local/lib/python3.12/dist-packages/rdkit/Contrib/fraggle/data/ChEMBL_11265_actives.smi", "RDKit/ChEMBL_frag"),
    ]
    for path, src in paths:
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    records.append({"smi": parts[0], "source": src, "name": ""})
    return records


# ── Cache ─────────────────────────────────────────────────────────────────────

def load_cache() -> Optional[List[Dict]]:
    if not os.path.exists(CACHE_FILE):
        return None
    try:
        with open(CACHE_FILE) as f:
            data = json.load(f)
        if data.get("version") != CACHE_VERSION:
            print("  Cache version mismatch, will re-download.")
            return None
        records = data["records"]
        ts = data.get("timestamp", "unknown")
        print(f"  Loaded {len(records)} cached molecules (downloaded: {ts})")
        return records
    except Exception as e:
        print(f"  Cache read error ({e}), will re-download.")
        return None


def save_cache(records: List[Dict]):
    with open(CACHE_FILE, "w") as f:
        json.dump({
            "version":   CACHE_VERSION,
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            "records":   records,
        }, f, separators=(",", ":"))
    print(f"  Saved {len(records)} molecules to {CACHE_FILE}")


# ── Data collection ───────────────────────────────────────────────────────────

def collect_dataset(quick: bool = False, no_cache: bool = False) -> List[Dict]:
    """Download and deduplicate a diverse 5000+ molecule dataset."""

    if not no_cache:
        cached = load_cache()
        if cached and len(cached) >= (500 if quick else 5000):
            if quick:
                return random.sample(cached, 500)
            return cached

    print("\n  Fetching from ChEMBL...")
    chembl = fetch_chembl(200 if quick else TARGET_CHEMBL)

    print("  Fetching from PubChem...")
    pubchem = fetch_pubchem(300 if quick else TARGET_PUBCHEM)

    print("  Loading RDKit local fallback...")
    local = load_rdkit_local()

    all_records = chembl + pubchem + local

    # Deduplicate by canonical SMILES string (fast dedup before InChI dedup)
    seen_smi = set()
    deduped = []
    for rec in all_records:
        smi = rec["smi"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        canon = Chem.MolToSmiles(mol)
        if canon not in seen_smi:
            seen_smi.add(canon)
            rec["canon_smi"] = canon
            deduped.append(rec)

    print(f"  After deduplication: {len(deduped)} unique molecules "
          f"(from {len(all_records)} raw records)")

    if not quick:
        save_cache(deduped)

    if quick:
        return random.sample(deduped, min(500, len(deduped)))
    return deduped


# ── Benchmark ─────────────────────────────────────────────────────────────────

def run_benchmark(records: List[Dict]) -> Dict:
    """
    For each record:
      1. Parse canonical SMILES with RDKit
      2. Convert RDKit mol -> SCRIPT string
      3. Convert SCRIPT string -> RDKit mol
      4. Compare InChI of original vs round-tripped
    """
    stats = {
        "total":       len(records),
        "pass":        0,
        "fail":        0,
        "skip_rdkit":  0,
        "skip_script": 0,
        "by_source":   defaultdict(lambda: {"pass": 0, "fail": 0, "skip": 0}),
        "failures":    [],
    }

    t0 = time.time()
    report_every = max(1, len(records) // 20)     # progress every 5%

    for i, rec in enumerate(records):
        src  = rec["source"]
        smi  = rec.get("canon_smi") or rec["smi"]

        if i % report_every == 0:
            elapsed = time.time() - t0
            done    = stats["pass"] + stats["fail"]
            rate    = done / elapsed if elapsed > 0 else 0
            eta     = (len(records) - i) / rate if rate > 0 else 0
            print(f"  [{i:>5}/{len(records)}]  pass={stats['pass']}  "
                  f"fail={stats['fail']}  "
                  f"rate={rate:.1f}/s  ETA={eta:.0f}s", end="\r")

        # Step 1: RDKit parse
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            stats["skip_rdkit"] += 1
            stats["by_source"][src]["skip"] += 1
            continue
        inchi_in = rdInchi.MolToInchi(mol)
        if not inchi_in:
            stats["skip_rdkit"] += 1
            stats["by_source"][src]["skip"] += 1
            continue

        # Step 2: RDKit -> SCRIPT
        try:
            script_str = SCRIPTFromMol(mol)
        except Exception as e:
            stats["skip_script"] += 1
            stats["by_source"][src]["skip"] += 1
            if len(stats["failures"]) < 50:
                stats["failures"].append({
                    "smi": smi, "source": src,
                    "stage": "to_script", "error": str(e)[:120],
                })
            continue

        if not script_str:
            stats["skip_script"] += 1
            stats["by_source"][src]["skip"] += 1
            continue

        # Step 3: SCRIPT -> RDKit
        try:
            mol2 = MolFromSCRIPT(script_str)
        except Exception as e:
            stats["fail"] += 1
            stats["by_source"][src]["fail"] += 1
            if len(stats["failures"]) < 50:
                stats["failures"].append({
                    "smi": smi, "script": script_str, "source": src,
                    "stage": "from_script", "error": str(e)[:120],
                })
            continue

        if mol2 is None:
            stats["fail"] += 1
            stats["by_source"][src]["fail"] += 1
            if len(stats["failures"]) < 50:
                stats["failures"].append({
                    "smi": smi, "script": script_str, "source": src,
                    "stage": "from_script", "error": "MolFromSCRIPT returned None",
                })
            continue

        # Step 4: InChI comparison
        inchi_out = rdInchi.MolToInchi(mol2)
        if inchi_in == inchi_out:
            stats["pass"] += 1
            stats["by_source"][src]["pass"] += 1
        else:
            stats["fail"] += 1
            stats["by_source"][src]["fail"] += 1
            if len(stats["failures"]) < 50:
                stats["failures"].append({
                    "smi": smi, "script": script_str, "source": src,
                    "stage": "inchi_mismatch",
                })

    print()   # clear the \r progress line
    stats["elapsed"] = round(time.time() - t0, 2)
    stats["tested"]  = stats["pass"] + stats["fail"]
    stats["pass_rate"] = (
        round(100 * stats["pass"] / stats["tested"], 4)
        if stats["tested"] > 0 else 0.0
    )
    return stats


# ── Failure analysis ──────────────────────────────────────────────────────────

def analyse_failures(failures: List[Dict]) -> Dict:
    """Categorise failures by stage and flag organometallic patterns."""
    by_stage  = defaultdict(int)
    metal_pat = {"Cu", "Fe", "Zn", "Co", "Ni", "Pd", "Pt", "Ag", "Au",
                 "Rh", "Ir", "Ru", "Os", "Mn", "Cr", "V", "Ti", "Mo", "W"}
    n_metal   = 0
    n_organic = 0

    for f in failures:
        by_stage[f.get("stage", "unknown")] += 1
        smi = f.get("smi", "")
        if any(m in smi for m in metal_pat):
            n_metal += 1
        else:
            n_organic += 1

    return {
        "by_stage":  dict(by_stage),
        "metal":     n_metal,
        "organic":   n_organic,
    }


# ── Report ────────────────────────────────────────────────────────────────────

def write_report(stats: Dict, dataset_meta: Dict) -> str:
    sep = "=" * 72
    thin = "-" * 72
    lines = [sep, "SCRIPT LARGE-SCALE BENCHMARK REPORT", sep, ""]

    # metadata
    lines += [
        f"  Date        : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"  Molecules   : {stats['total']} downloaded, "
                        f"{stats['tested']} tested "
                        f"({stats['skip_rdkit']+stats['skip_script']} skipped)",
        f"  Sources     : {', '.join(dataset_meta.get('sources', []))}",
        "",
    ]

    # headline result
    lines += [
        thin,
        "ROUND-TRIP RESULTS  (SMILES -> SCRIPT -> InChI)",
        thin,
        f"  Passed      : {stats['pass']:>6}",
        f"  Failed      : {stats['fail']:>6}",
        f"  Skipped*    : {stats['skip_rdkit']+stats['skip_script']:>6}",
        f"  Tested      : {stats['tested']:>6}",
        f"  PASS RATE   : {stats['pass_rate']:.2f}%",
        f"  Wall time   : {stats['elapsed']:.1f}s  "
                        f"({stats['tested']/stats['elapsed']:.1f} mol/s)",
        "",
        "  * Skipped = RDKit rejected SMILES before SCRIPT was invoked.",
        "    These are not SCRIPT failures.",
        "",
    ]

    # per-source table
    lines += [
        thin,
        "BY SOURCE",
        thin,
        f"  {'Source':<25} {'Pass':>6} {'Fail':>6} {'Skip':>6} {'Rate':>8}",
        f"  {'-'*25} {'-'*6} {'-'*6} {'-'*6} {'-'*8}",
    ]
    for src, c in sorted(stats["by_source"].items()):
        tested = c["pass"] + c["fail"]
        rate   = f"{100*c['pass']/tested:.1f}%" if tested else "n/a"
        lines.append(
            f"  {src:<25} {c['pass']:>6} {c['fail']:>6} {c['skip']:>6} {rate:>8}"
        )
    lines.append("")

    # failure analysis
    if stats["failures"]:
        fa = analyse_failures(stats["failures"])
        lines += [
            thin,
            "FAILURE ANALYSIS",
            thin,
            f"  By stage: {fa['by_stage']}",
            f"  Organometallic failures : {fa['metal']}  "
                f"(RDKit bridge coordinate-bond limitation)",
            f"  Organic failures        : {fa['organic']}",
            "",
            f"  First 10 failures:",
        ]
        for f in stats["failures"][:10]:
            lines.append(
                f"    [{f.get('stage','?'):>15}]  {f.get('smi','?')[:55]}"
            )
            if "error" in f:
                lines.append(f"                         {f['error'][:65]}")
        lines.append("")

    # summary for paper
    lines += [
        sep,
        "SUMMARY  (for paper)",
        sep,
        f"  Organic round-trip pass rate : {stats['pass_rate']:.2f}%  "
                        f"({stats['pass']}/{stats['tested']})",
        f"  Dataset diversity            : "
                        f"{len(stats['by_source'])} sources, "
                        f"{stats['tested']} molecules",
        f"  Largest failure category     : "
                        f"{max(analyse_failures(stats['failures'])['by_stage'].items(), key=lambda x: x[1])[0] if stats['failures'] else 'none'}",
        "",
    ]

    return "\n".join(lines)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="SCRIPT large-scale benchmark (5000+ compounds)"
    )
    parser.add_argument("--quick",    action="store_true",
                        help="Quick mode: ~500 compounds, no new downloads")
    parser.add_argument("--no-cache", action="store_true",
                        help="Force re-download even if cache exists")
    args = parser.parse_args()

    print("=" * 72)
    print("SCRIPT LARGE-SCALE BENCHMARK")
    print("=" * 72)

    # ── 1. Collect dataset ────────────────────────────────────────────────────
    print("\n[1/3] Collecting dataset...")
    records = collect_dataset(quick=args.quick, no_cache=args.no_cache)
    sources = sorted(set(r["source"] for r in records))
    print(f"  Final dataset: {len(records)} molecules from {sources}")

    # ── 2. Run benchmark ──────────────────────────────────────────────────────
    print(f"\n[2/3] Running benchmark on {len(records)} molecules...")
    stats = run_benchmark(records)
    print(f"  Done: {stats['pass']}/{stats['tested']} passed "
          f"({stats['pass_rate']:.2f}%) in {stats['elapsed']:.1f}s")

    # ── 3. Write report ───────────────────────────────────────────────────────
    print("\n[3/3] Writing report...")
    dataset_meta = {"sources": sources}
    report_text = write_report(stats, dataset_meta)
    print(report_text)

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    txt_path  = f"large_benchmark_{ts}.txt"
    json_path = f"large_benchmark_{ts}.json"

    with open(txt_path, "w") as f:
        f.write(report_text)

    json_out = {
        "timestamp":    ts,
        "dataset_meta": dataset_meta,
        "stats": {
            **stats,
            "by_source": dict(stats["by_source"]),
        },
    }
    with open(json_path, "w") as f:
        json.dump(json_out, f, indent=2, default=str)

    print(f"\nReports written:")
    print(f"  {txt_path}")
    print(f"  {json_path}")
    if not args.quick and not args.no_cache:
        print(f"  {CACHE_FILE}  (molecule cache — reused on next run)")


if __name__ == "__main__":
    random.seed(42)
    main()
