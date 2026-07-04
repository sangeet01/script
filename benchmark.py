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
from script.mol import CoreMolecule


def _iter_molecules(mol_or_list):
    """Yield CoreMolecule instances from either a single mol or a list.

    The SCRIPT parser returns a `list[CoreMolecule]` for `.` / `~` separated
    multi-fragment inputs (salts, solvates, ionic pairs) and a single
    CoreMolecule otherwise. Benchmark accounting needs to aggregate atom and
    bond counts across all fragments, so we normalise both shapes to an
    iterable here.
    """
    if mol_or_list is None:
        return []
    if isinstance(mol_or_list, list):
        return [m for m in mol_or_list if isinstance(m, CoreMolecule)]
    if isinstance(mol_or_list, CoreMolecule):
        return [mol_or_list]
    return []


def _atom_bond_counts(mol_or_list):
    """Return (total_atoms, total_bonds) across all fragments."""
    mols = _iter_molecules(mol_or_list)
    n_atoms = sum(len(m.atoms) for m in mols)
    n_bonds = sum(len(m.bonds) for m in mols)
    return n_atoms, n_bonds

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
            mols_in = _iter_molecules(mol)
            if not mols_in:
                print(f"{name:<20} | SKIP | No molecule generated")
                skipped += 1
                continue

            # canonicalize_mols handles both a single CoreMolecule and a
            # list[CoreMolecule] (multi-fragment salts / solvates / ionic
            # pairs). Using canonicalize_core here would crash on the list
            # shape — that was the original benchmark bug.
            script_out = canonicalizer.canonicalize_mols(mol)
            if script_out is None:
                print(f"{name:<20} | FAIL | Canonicalization returned None")
                failed += 1
                continue

            result2 = parser.parse(script_out)
            if not result2["success"]:
                print(f"{name:<20} | FAIL | Round-trip parse error: {result2['error']}")
                failed += 1
                continue

            atoms_in, bonds_in = _atom_bond_counts(result["molecule"])
            atoms_out, bonds_out = _atom_bond_counts(result2["molecule"])

            if atoms_in != atoms_out:
                print(f"{name:<20} | FAIL | Atom count mismatch: {atoms_in} vs {atoms_out}")
                print(f"  Input:  {script_in}")
                print(f"  Output: {script_out}")
                failed += 1
                continue

            if bonds_in != bonds_out:
                print(f"{name:<20} | FAIL | Bond count mismatch: {bonds_in} vs {bonds_out}")
                failed += 1
                continue

            frag_note = f" ({len(mols_in)} fragments)" if len(mols_in) > 1 else ""
            print(f"{name:<20} | PASS | Round-trip OK ({atoms_in} atoms, {bonds_in} bonds{frag_note})")
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
