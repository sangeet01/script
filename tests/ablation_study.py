#!/usr/bin/env python3
"""
SCRIPT Ablation Study
======================
Component-wise validation of SCRIPT's four core mechanisms:
  1. Sandhi valence validation
  2. Deterministic canonical ordering
  3. CIP-anchored stereochemistry
  4. Topological ring closures

Produces a structured report matching Section 4 (Ablation Study) of the paper.

Usage:
    python ablation_study.py              # full run
    python ablation_study.py --json       # also write ablation_results.json
    python ablation_study.py --verbose    # show all intermediate values
"""

import sys
import os
import json
import argparse
import datetime
import io

# ── path setup ────────────────────────────────────────────────────────────────
# Fix Windows console encoding for Unicode characters
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from rdkit import Chem
    from rdkit.Chem import inchi as rdInchi
    from rdkit import RDLogger
except ImportError:
    sys.exit("ERROR: rdkit not found. pip install rdkit")

try:
    from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT
    from script.parser import SCRIPTParser
except ImportError:
    sys.exit("ERROR: SCRIPT not found. Run from repo root or: pip install -e .")

# ── helpers ───────────────────────────────────────────────────────────────────

SEP  = "=" * 70
THIN = "-" * 70

def header(title):
    print(f"\n{SEP}")
    print(f"  {title}")
    print(SEP)

def subheader(title):
    print(f"\n{THIN}")
    print(f"  {title}")
    print(THIN)


# ─────────────────────────────────────────────────────────────────────────────
# ABLATION 1 — SANDHI VALENCE VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

INVALID_STRUCTURES = [
    ("6-valent carbon",     "C(C)(C)(C)(C)(C)C",      6,  "C"),
    ("5-valent nitrogen",   "[N](C)(C)(C)(C)C",        5,  "N"),
    ("3-valent oxygen",     "[O](C)(C)C",              3,  "O"),
    ("7-valent phosphorus", "[P](C)(C)(C)(C)(C)(C)C",  7,  "P"),
    ("5-valent boron",      "[B](C)(C)(C)(C)C",        5,  "B"),
    ("9-valent sulfur",     "[S](C)(C)(C)(C)(C)(C)(C)(C)C", 9, "S"),
]

def run_sandhi_ablation(parser, verbose=False):
    """
    Submit chemically impossible structures to both RDKit and SCRIPT.
    SCRIPT should reject all via Sandhi; document the error message.
    """
    results = []
    passed = 0

    for name, smiles, valence, element in INVALID_STRUCTURES:
        # RDKit behaviour
        mol = Chem.MolFromSmiles(smiles)
        rdkit_accepts = mol is not None

        # SCRIPT behaviour
        r = parser.parse(smiles)
        script_accepts = r.get("success", False)
        error_msg = r.get("error", "")

        correct = not script_accepts          # correct = SCRIPT rejects
        if correct:
            passed += 1

        results.append({
            "name":           name,
            "smiles":         smiles,
            "valence":        valence,
            "element":        element,
            "rdkit_accepts":  rdkit_accepts,
            "script_accepts": script_accepts,
            "script_error":   error_msg,
            "correct":        correct,
        })

        if verbose:
            rdkit_str  = "accept" if rdkit_accepts  else "reject"
            script_str = "ACCEPT (BUG)" if script_accepts else "REJECT ✓"
            print(f"  {name:<28} RDKit={rdkit_str:<7}  SCRIPT={script_str}"
                  f"  [{error_msg[:45]}]")

    return results, passed


# ─────────────────────────────────────────────────────────────────────────────
# ABLATION 2 — DETERMINISTIC CANONICAL ORDERING
# ─────────────────────────────────────────────────────────────────────────────

CANONICALITY_CASES = [
    ("Aspirin", [
        "CC(=O)Oc1ccccc1C(=O)O",
        "OC(=O)c1ccccc1OC(C)=O",
        "O=C(C)Oc1ccccc1C(=O)O",
    ]),
    ("Ibuprofen", [
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "OC(=O)C(C)c1ccc(CC(C)C)cc1",
        "CC(C)CC1=CC=C(C=C1)C(C)C(O)=O",
    ]),
    ("Naphthalene", [
        "c1ccc2ccccc2c1",
        "c1cccc2ccccc12",
        "C1=CC=CC2=CC=CC=C12",
    ]),
    ("Caffeine", [
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        "CN1C(=O)N(C)c2ncn(C)c2C1=O",
    ]),
    ("Testosterone", [
        "[C@@H]12(CC[C@H](O)C1)CC[C@@H]1[C@H]2CCC2=CC(=O)CC[C@@]12C",
        "O[C@@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@@]12C",
    ]),
    # Deliberate adversarial case: α vs β glucose ARE different molecules
    ("Glucose (α vs β — should differ)", [
        "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",   # alpha
        "OC[C@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",   # beta
    ]),
]

def run_canonicality_ablation(verbose=False):
    """
    For each molecule, convert multiple SMILES to SCRIPT and check uniqueness.
    The glucose case is expected to produce 2 distinct strings (they ARE different).
    """
    results = []
    passed = 0

    for name, smiles_list in CANONICALITY_CASES:
        scripts = []
        rdkit_canons = set()
        is_adversarial = "should differ" in name.lower()

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                continue
            rdkit_canons.add(Chem.MolToSmiles(mol))
            try:
                s = SCRIPTFromMol(mol)
                if s:
                    scripts.append(s)
            except Exception as e:
                scripts.append(f"ERROR: {e}")

        unique_scripts = set(s for s in scripts if not s.startswith("ERROR"))
        unique_rdkit   = len(rdkit_canons)

        if is_adversarial:
            # Correct = produces 2 distinct strings (different molecules)
            correct = len(unique_scripts) == len(smiles_list)
            expected_unique = len(smiles_list)
        else:
            # Correct = all inputs produce 1 canonical string
            correct = len(unique_scripts) == 1
            expected_unique = 1

        if correct:
            passed += 1

        results.append({
            "name":            name,
            "n_inputs":        len(smiles_list),
            "rdkit_unique":    unique_rdkit,
            "script_unique":   len(unique_scripts),
            "expected_unique": expected_unique,
            "script_strings":  list(unique_scripts),
            "correct":         correct,
            "adversarial":     is_adversarial,
        })

        if verbose:
            status = "PASS ✓" if correct else "FAIL ✗"
            note = " [adversarial — correct to differ]" if is_adversarial else ""
            print(f"  {name:<40} inputs={len(smiles_list)}"
                  f"  SMILES_unique={unique_rdkit}"
                  f"  SCRIPT_unique={len(unique_scripts)}"
                  f"  {status}{note}")

    return results, passed


# ─────────────────────────────────────────────────────────────────────────────
# ABLATION 3 — CIP-ANCHORED STEREOCHEMISTRY
# ─────────────────────────────────────────────────────────────────────────────

STEREO_CASES = [
    ("L-Alanine",          "N[C@@H](C)C(=O)O"),
    ("D-Alanine",          "N[C@H](C)C(=O)O"),
    ("alpha-D-Glucose",    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("beta-D-Glucose",     "OC[C@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O"),
    ("(E)-2-Butene",       "C/C=C/C"),
    ("(Z)-2-Butene",       "C/C=C\\C"),
    ("Thalidomide-(R)",    "O=C1CCC(=O)N1[C@@H]1C(=O)Nc2ccccc21"),
    ("Thalidomide-(S)",    "O=C1CCC(=O)N1[C@H]1C(=O)Nc2ccccc21"),
    ("(R)-Ibuprofen",      "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"),
    ("(S)-Ibuprofen",      "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
]

ENANTIOMERIC_PAIRS = [
    ("L-Alanine",       "D-Alanine"),
    ("alpha-D-Glucose", "beta-D-Glucose"),
    ("Thalidomide-(R)", "Thalidomide-(S)"),
    ("(E)-2-Butene",    "(Z)-2-Butene"),
    ("(R)-Ibuprofen",   "(S)-Ibuprofen"),
]

def run_stereo_ablation(verbose=False):
    """
    Round-trip each stereoisomer through SCRIPT and verify:
      (a) InChI is preserved (no stereo information lost)
      (b) Enantiomeric pairs produce distinct SCRIPT strings
    """
    roundtrip_results = []
    script_map = {}
    rt_passed = 0

    for name, smiles in STEREO_CASES:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            roundtrip_results.append({
                "name": name, "smiles": smiles,
                "error": "RDKit rejected SMILES", "correct": False,
            })
            continue

        inchi_in = rdInchi.MolToInchi(mol)

        try:
            script_str = SCRIPTFromMol(mol)
            mol2       = MolFromSCRIPT(script_str) if script_str else None
            inchi_out  = rdInchi.MolToInchi(mol2) if mol2 else None
            match      = (inchi_in == inchi_out) if inchi_out else False
        except Exception as e:
            script_str = None
            inchi_out  = None
            match      = False

        if match:
            rt_passed += 1
        if script_str:
            script_map[name] = script_str

        roundtrip_results.append({
            "name":       name,
            "smiles":     smiles,
            "script":     script_str or "ERROR",
            "inchi_match": match,
            "correct":    match,
        })

        if verbose:
            status = "PASS ✓" if match else "FAIL ✗"
            s_excerpt = (script_str or "ERROR")[:55]
            print(f"  {name:<22} {status}  SCRIPT={s_excerpt}")

    # Enantiomeric pair distinctness
    pair_results = []
    pair_passed  = 0

    for name_a, name_b in ENANTIOMERIC_PAIRS:
        s_a = script_map.get(name_a)
        s_b = script_map.get(name_b)

        if s_a is None or s_b is None:
            pair_results.append({
                "pair": f"{name_a} / {name_b}",
                "distinct": None,
                "correct": False,
                "note": "one or both SCRIPT strings missing",
            })
            continue

        distinct = s_a != s_b
        if distinct:
            pair_passed += 1

        pair_results.append({
            "pair":     f"{name_a} / {name_b}",
            "script_a": s_a,
            "script_b": s_b,
            "distinct": distinct,
            "correct":  distinct,
        })

        if verbose:
            status = "DISTINCT ✓" if distinct else "COLLAPSED ✗ (BUG)"
            print(f"  {name_a:<22} vs {name_b:<22} {status}")

    return roundtrip_results, rt_passed, pair_results, pair_passed


# ─────────────────────────────────────────────────────────────────────────────
# ABLATION 4 — TOPOLOGICAL RING CLOSURES
# ─────────────────────────────────────────────────────────────────────────────

RING_CASES = [
    ("Cyclopropane",  "C1CC1",                  3,  "3-membered ring"),
    ("Cyclobutane",   "C1CCC1",                 4,  "4-membered ring"),
    ("Cyclopentane",  "C1CCCC1",                5,  "5-membered ring"),
    ("Cyclohexane",   "C1CCCCC1",               6,  "6-membered ring"),
    ("Benzene",       "c1ccccc1",               6,  "aromatic 6-ring"),
    ("Naphthalene",   "c1ccc2ccccc2c1",         10, "fused bicyclic"),
    ("Indole",        "c1ccc2[nH]ccc2c1",       9,  "heteroaromatic fused"),
    ("Cubane",        "C12C3C4C1C1C2C3C41",     8,  "cage polycyclic"),
    ("Adamantane",    "C1C2CC3CC1CC(C2)C3",     10, "cage polycyclic"),
    ("Steroid core",  "C1CCC2CCCCC2C1",         10, "fused decalin"),
]

def run_ring_closure_ablation(verbose=False):
    """
    For each cyclic system, show the SMILES ring label notation vs
    SCRIPT's topological &N notation. Verify round-trip fidelity.
    """
    results = []
    passed  = 0

    for name, smiles, ring_size, desc in RING_CASES:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            results.append({
                "name": name, "error": "RDKit rejected", "correct": False
            })
            continue

        canon_smiles = Chem.MolToSmiles(mol)
        inchi_in     = rdInchi.MolToInchi(mol)

        try:
            script_str = SCRIPTFromMol(mol)
            mol2       = MolFromSCRIPT(script_str) if script_str else None
            inchi_out  = rdInchi.MolToInchi(mol2) if mol2 else None
            match      = (inchi_in == inchi_out) if inchi_out else False
        except Exception as e:
            script_str = None
            match      = False

        # Check that SCRIPT uses &N notation, not SMILES-style numeric labels
        uses_local_notation = "&" in (script_str or "")

        correct = match
        if correct:
            passed += 1

        results.append({
            "name":               name,
            "ring_size":          ring_size,
            "description":        desc,
            "smiles":             canon_smiles,
            "script":             script_str or "ERROR",
            "inchi_match":        match,
            "uses_local_notation": uses_local_notation,
            "correct":            correct,
        })

        if verbose:
            status = "PASS ✓" if match else "FAIL ✗"
            s = (script_str or "ERROR")[:45]
            print(f"  {name:<15} {canon_smiles:<28} -> {s:<45} {status}")

    return results, passed


# ─────────────────────────────────────────────────────────────────────────────
# REPORT WRITER
# ─────────────────────────────────────────────────────────────────────────────

def write_report(sandhi, sandhi_p,
                 canon, canon_p,
                 stereo_rt, stereo_rt_p, stereo_pairs, stereo_pair_p,
                 rings, rings_p):

    lines = []
    lines += [SEP, "SCRIPT ABLATION STUDY — RESULTS REPORT", SEP, ""]
    lines += [
        f"  Generated : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"  Purpose   : Component-wise validation for paper Section 4",
        "",
    ]

    # ── Ablation 1 ────────────────────────────────────────────────────────────
    lines += [
        THIN,
        "ABLATION 1 — SANDHI VALENCE VALIDATION",
        THIN,
        f"  {'Structure':<28} {'RDKit':>7} {'SCRIPT':>14} {'Error (excerpt)'}",
        f"  {'-'*28} {'-'*7} {'-'*14} {'-'*35}",
    ]
    for r in sandhi:
        rdkit_str  = "accept" if r["rdkit_accepts"]  else "reject"
        script_str = "accept" if r["script_accepts"] else "REJECT"
        err = r["script_error"][:38] if r["script_error"] else ""
        lines.append(
            f"  {r['name']:<28} {rdkit_str:>7} {script_str:>14} {err}"
        )
    lines += [
        "",
        f"  Result: {sandhi_p}/{len(sandhi)} structures correctly rejected by SCRIPT",
        "",
    ]

    # ── Ablation 2 ────────────────────────────────────────────────────────────
    lines += [
        THIN,
        "ABLATION 2 — DETERMINISTIC CANONICAL ORDERING",
        THIN,
        f"  {'Molecule':<42} {'Inputs':>6} {'SMILES uniq':>12} {'SCRIPT uniq':>12} {'Result':>8}",
        f"  {'-'*42} {'-'*6} {'-'*12} {'-'*12} {'-'*8}",
    ]
    for r in canon:
        note = " *" if r["adversarial"] else ""
        status = "PASS" if r["correct"] else "FAIL"
        lines.append(
            f"  {r['name']+note:<42} {r['n_inputs']:>6} "
            f"{r['rdkit_unique']:>12} {r['script_unique']:>12} {status:>8}"
        )
    lines += [
        "",
        "  * Glucose α vs β are genuinely different molecules (correct to differ).",
        f"  Result: {canon_p}/{len(canon)} cases correct",
        "",
    ]

    # ── Ablation 3 ────────────────────────────────────────────────────────────
    lines += [
        THIN,
        "ABLATION 3 — CIP-ANCHORED STEREOCHEMISTRY",
        THIN,
        "  3a. Round-trip fidelity (SMILES -> SCRIPT -> InChI)",
        "",
        f"  {'Molecule':<22} {'InChI match':>12} {'SCRIPT (excerpt)'}",
        f"  {'-'*22} {'-'*12} {'-'*45}",
    ]
    for r in stereo_rt:
        match_str = "YES ✓" if r.get("inchi_match") else "NO ✗"
        s = r.get("script", "ERROR")[:45]
        lines.append(f"  {r['name']:<22} {match_str:>12}   {s}")
    lines += [
        "",
        f"  Round-trip result: {stereo_rt_p}/{len(stereo_rt)} stereoisomers preserved",
        "",
        "  3b. Enantiomeric pair distinctness",
        "",
        f"  {'Pair':<48} {'Distinct':>10}",
        f"  {'-'*48} {'-'*10}",
    ]
    for r in stereo_pairs:
        dist_str = "YES ✓" if r.get("distinct") else ("NO ✗" if r.get("distinct") is False else "N/A")
        lines.append(f"  {r['pair']:<48} {dist_str:>10}")
    lines += [
        "",
        f"  Distinctness result: {stereo_pair_p}/{len(stereo_pairs)} pairs correctly distinct",
        "",
    ]

    # ── Ablation 4 ────────────────────────────────────────────────────────────
    lines += [
        THIN,
        "ABLATION 4 — TOPOLOGICAL RING CLOSURES",
        THIN,
        f"  {'Molecule':<15} {'SMILES':<30} {'SCRIPT':<45} {'RT':>5}",
        f"  {'-'*15} {'-'*30} {'-'*45} {'-'*5}",
    ]
    for r in rings:
        if "error" in r:
            lines.append(f"  {r['name']:<15} ERROR")
            continue
        match_str = "OK" if r["inchi_match"] else "FAIL"
        lines.append(
            f"  {r['name']:<15} {r['smiles']:<30} {r['script']:<45} {match_str:>5}"
        )
    lines += [
        "",
        f"  Result: {rings_p}/{len(rings)} ring systems round-trip correctly",
        "",
    ]

    # ── Summary ───────────────────────────────────────────────────────────────
    total_tests = (len(sandhi) + len(canon) +
                   len(stereo_rt) + len(stereo_pairs) + len(rings))
    total_passed = sandhi_p + canon_p + stereo_rt_p + stereo_pair_p + rings_p

    lines += [
        SEP,
        "SUMMARY",
        SEP,
        f"  {'Component':<40} {'Passed':>8} {'Total':>8}",
        f"  {'-'*40} {'-'*8} {'-'*8}",
        f"  {'Sandhi valence validation':<40} {sandhi_p:>8} {len(sandhi):>8}",
        f"  {'Canonical ordering':<40} {canon_p:>8} {len(canon):>8}",
        f"  {'Stereo round-trip':<40} {stereo_rt_p:>8} {len(stereo_rt):>8}",
        f"  {'Enantiomeric distinctness':<40} {stereo_pair_p:>8} {len(stereo_pairs):>8}",
        f"  {'Ring closure round-trip':<40} {rings_p:>8} {len(rings):>8}",
        f"  {'-'*40} {'-'*8} {'-'*8}",
        f"  {'TOTAL':<40} {total_passed:>8} {total_tests:>8}",
        "",
        f"  Overall pass rate: {100*total_passed/total_tests:.1f}% ({total_passed}/{total_tests})",
        "",
    ]

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="SCRIPT ablation study — validates four core components"
    )
    ap.add_argument("--json",    action="store_true", help="Write JSON results file")
    ap.add_argument("--verbose", action="store_true", help="Show per-case detail during run")
    args = ap.parse_args()

    parser = SCRIPTParser()

    print(SEP)
    print("  SCRIPT ABLATION STUDY")
    print(SEP)

    # ── Run all four ablations ────────────────────────────────────────────────

    print("\n[1/4] Sandhi valence validation...")
    sandhi, sandhi_p = run_sandhi_ablation(parser, verbose=args.verbose)
    print(f"      {sandhi_p}/{len(sandhi)} correctly rejected")

    print("\n[2/4] Canonical ordering...")
    canon, canon_p = run_canonicality_ablation(verbose=args.verbose)
    print(f"      {canon_p}/{len(canon)} canonical")

    print("\n[3/4] CIP stereo round-trip + enantiomeric distinctness...")
    stereo_rt, stereo_rt_p, stereo_pairs, stereo_pair_p = run_stereo_ablation(
        verbose=args.verbose
    )
    print(f"      Round-trip: {stereo_rt_p}/{len(stereo_rt)}   "
          f"Pairs distinct: {stereo_pair_p}/{len(stereo_pairs)}")

    print("\n[4/4] Topological ring closures...")
    rings, rings_p = run_ring_closure_ablation(verbose=args.verbose)
    print(f"      {rings_p}/{len(rings)} ring systems correct")

    # ── Print and save report ────────────────────────────────────────────────

    report = write_report(
        sandhi, sandhi_p,
        canon, canon_p,
        stereo_rt, stereo_rt_p, stereo_pairs, stereo_pair_p,
        rings, rings_p,
    )

    print("\n" + report)

    txt_path = "ablation_results.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"Report saved to {txt_path}")

    if args.json:
        json_path = "ablation_results.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump({
                "timestamp":        datetime.datetime.now().isoformat(),
                "sandhi":           sandhi,
                "canonicality":     canon,
                "stereo_roundtrip": stereo_rt,
                "stereo_pairs":     stereo_pairs,
                "ring_closures":    rings,
            }, f, indent=2, default=str)
        print(f"JSON saved to {json_path}")


if __name__ == "__main__":
    main()
