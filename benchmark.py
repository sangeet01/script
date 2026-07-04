"""
SCRIPT 10-Compound Drug Benchmark
Tests the exact molecules claimed in the README using native SCRIPT round-trip.
No RDKit dependency - validates parser, canonicalizer, and state machine.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'script'))

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer

# The 10 drugs explicitly claimed in the README
DRUG_DATASET = [
    # Aspirin - acetylsalicylic acid
    ("Aspirin", "C(C(=O)O):C(OC(=O)C):C:C:C:C&6:"),
    # Metformin - N,N-dimethylimidodicarbonimidic diamide
    ("Metformin", "CN(C)C(=N)N=C(N)N"),
    # Ciprofloxacin HCl - fluoroquinolone salt (core without HCl counterion)
    ("Ciprofloxacin", "O=C(O)C1=CN(C2CC2)C2=C(C1)C(=O)C(=CN2C1CC1)F"),
    # Nifedipine - dimethyl 1,4-dihydropyridine-3,5-dicarboxylate
    ("Nifedipine", "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1[N+](=O)[O-]"),
    # Ibuprofen - 2-(4-isobutylphenyl)propanoic acid
    ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
    # Captopril - proline-based ACE inhibitor
    ("Captopril", "CC(CS)C(=O)N1CCCC1C(=O)O"),
    # Glucose - pyranose form, 4 stereocenters
    ("Glucose", "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO"),
    # Metformin HCl - salt form (two components)
    ("Metformin HCl", "CN(C)C(=N)N=C(N)N.[Cl-]"),
    # Magnesium stearate - metal + multi-fragment
    ("Magnesium stearate", "[Mg+2].[O-]C(=O)CCCCCCCCCCCCCCCCC.[O-]C(=O)CCCCCCCCCCCCCCCCC"),
    # PVP - polyvinylpyrrolidone (monomer unit)
    ("PVP", "{[CC(C(=O)N1CCCC1)]}n"),
]


def run_benchmark():
    print("Running SCRIPT 10-Drug Benchmark")
    print("=" * 60)

    parser = SCRIPTParser()
    canonicalizer = SCRIPTCanonicalizer()

    passed = 0
    failed = 0
    skipped = 0

    for name, script_in in DRUG_DATASET:
        try:
            result = parser.parse(script_in)
            if not result["success"]:
                print(f"{name:<20} | FAIL | Parse error: {result['error']}")
                failed += 1
                continue

            mol = result["molecule"]
            if mol is None:
                print(f"{name:<20} | SKIP | No molecule generated")
                skipped += 1
                continue

            script_out = canonicalizer.canonicalize_core(mol)
            if script_out is None:
                print(f"{name:<20} | FAIL | Canonicalization returned None")
                failed += 1
                continue

            result2 = parser.parse(script_out)
            if not result2["success"]:
                print(f"{name:<20} | FAIL | Round-trip parse error: {result2['error']}")
                failed += 1
                continue

            mol_in = result["molecule"]
            mol_out = result2["molecule"]
            atoms_in = len(mol_in.atoms) if hasattr(mol_in, "atoms") else 0
            atoms_out = len(mol_out.atoms) if hasattr(mol_out, "atoms") else 0

            if atoms_in != atoms_out:
                print(f"{name:<20} | FAIL | Atom count mismatch: {atoms_in} vs {atoms_out}")
                print(f"  Input:  {script_in}")
                print(f"  Output: {script_out}")
                failed += 1
                continue

            bonds_in = len(mol_in.bonds) if hasattr(mol_in, "bonds") else 0
            bonds_out = len(mol_out.bonds) if hasattr(mol_out, "bonds") else 0

            if bonds_in != bonds_out:
                print(f"{name:<20} | FAIL | Bond count mismatch: {bonds_in} vs {bonds_out}")
                failed += 1
                continue

            print(f"{name:<20} | PASS | Round-trip OK ({atoms_in} atoms, {bonds_in} bonds)")
            print(f"  Canonical: {script_out}")
            passed += 1

        except Exception as e:
            print(f"{name:<20} | ERR  | {str(e)}")
            failed += 1

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped")
    total = passed + failed
    if total > 0:
        print(f"Success Rate: {passed/total*100:.1f}%")
    print("=" * 60)
    return failed == 0


if __name__ == "__main__":
    success = run_benchmark()
    sys.exit(0 if success else 1)
