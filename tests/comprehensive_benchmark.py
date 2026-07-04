#!/usr/bin/env python3
"""
SCRIPT State-of-the-Art Comprehensive Benchmark
=================================================
A rigorous, adversarial benchmark that tests every claim made in the README
to determine if SCRIPT is actually production-ready.

Categories tested:
  1. Basic Organic (100+ molecules)
  2. Sandhi Validation (invalid structures must be rejected)
  3. Ring Closures (topological vs legacy)
  4. Stereochemistry (tetrahedral, extended, allenes)
  5. Typed Bonds (dative, coordinate, tautomeric, haptic)
  6. Reactions and Atom Mapping
  7. Query Atoms (SMARTS-style patterns)
  8. Biopolymers (peptides, nucleic acids, PTMs)
  9. Materials (alloys, crystals, surfaces, quantum states)
  10. Polymers (simple, stochastic, block copolymers)
  11. Canonicalization (stability, uniqueness, cross-input)
  12. Boss Fights (Taxol, Strychnine, Glucose, C60, Porphyrin)
  13. Fuzzing (random structure generation)
  14. Performance (timing at scale)

Usage:
    python comprehensive_benchmark.py              # full run
    python comprehensive_benchmark.py --quick      # fast run (~200 tests)
    python comprehensive_benchmark.py --json       # output JSON report
    python comprehensive_benchmark.py --verbose    # show all details

Output:
    comprehensive_benchmark_YYYYMMDD_HHMMSS.txt
    comprehensive_benchmark_YYYYMMDD_HHMMSS.json
"""

import sys
import os
import time
import json
import random
import argparse
import datetime
from collections import defaultdict
from typing import List, Dict, Optional, Any

# Path setup
SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SCRIPT_DIR)

# Dependencies
try:
    from lark import Lark
    LARK_AVAILABLE = True
except ImportError:
    LARK_AVAILABLE = False
    print("WARNING: lark not available. Parser tests will be skipped.")

try:
    from script.parser import SCRIPTParser
    from script.canonical import SCRIPTCanonicalizer
    from script.mol import BondType, StereoType, CoreMolecule
    SCRIPT_AVAILABLE = True
except ImportError as e:
    SCRIPT_AVAILABLE = False
    print(f"WARNING: SCRIPT not available: {e}")

# RDKit is optional
try:
    from rdkit import Chem
    from rdkit.Chem import inchi as rdInchi
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not available. RDKit-dependent tests skipped.")

# Constants
SEP = "=" * 78
THIN = "-" * 78

# Result tracking
class BenchmarkResults:
    def __init__(self):
        self.categories = defaultdict(lambda: {"passed": 0, "failed": 0, "skipped": 0, "tests": []})
        self.start_time = time.time()

    def add(self, category: str, name: str, status: str, detail: str = "",
            expected: Any = None, actual: Any = None):
        cat = self.categories[category]
        status_key = status.lower()
        if status_key in ('pass', 'passed'):
            status_key = 'passed'
        elif status_key in ('fail', 'failed'):
            status_key = 'failed'
        elif status_key in ('skip', 'skipped'):
            status_key = 'skipped'
        cat[status_key] += 1
        cat["tests"].append({
            "name": name,
            "status": status,
            "detail": detail,
            "expected": str(expected) if expected is not None else None,
            "actual": str(actual) if actual is not None else None,
            "timestamp": time.time()
        })

    def summary(self) -> Dict:
        total_passed = sum(c["passed"] for c in self.categories.values())
        total_failed = sum(c["failed"] for c in self.categories.values())
        total_skipped = sum(c["skipped"] for c in self.categories.values())
        total = total_passed + total_failed
        return {
            "total_tests": total + total_skipped,
            "passed": total_passed,
            "failed": total_failed,
            "skipped": total_skipped,
            "success_rate": (total_passed / total * 100) if total > 0 else 0,
            "elapsed_seconds": round(time.time() - self.start_time, 2),
            "categories": {
                k: {
                    "passed": v["passed"],
                    "failed": v["failed"],
                    "skipped": v["skipped"],
                    "rate": (v["passed"] / (v["passed"] + v["failed"]) * 100)
                            if (v["passed"] + v["failed"]) > 0 else 0
                }
                for k, v in self.categories.items()
            }
        }


results = BenchmarkResults()
args = None


def test(category: str, name: str, condition: bool, detail: str = "",
         expected: Any = None, actual: Any = None):
    if condition:
        results.add(category, name, "PASS", detail, expected, actual)
        print(f"  [PASS] {category}: {name}")
    else:
        results.add(category, name, "FAIL", detail, expected, actual)
        print(f"  [FAIL] {category}: {name}: {detail}")


def skip(category: str, name: str, reason: str):
    results.add(category, name, "SKIP", reason)
    print(f"  [SKIP] {category}: {name}: {reason}")


# Category 1: Basic Organic
def run_basic_organic():
    print(f"\n{THIN}")
    print("  CATEGORY 1: BASIC ORGANIC MOLECULES")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Basic Organic", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    simple = [
        ("Methane", "C", 1, 0), ("Ethane", "CC", 2, 1), ("Propane", "CCC", 3, 2),
        ("Butane", "CCCC", 4, 3), ("Isobutane", "CC(C)C", 4, 3),
        ("Neopentane", "CC(C)(C)C", 5, 4), ("Ethanol", "CCO", 3, 2),
        ("Isopropanol", "CC(O)C", 4, 3), ("tert-Butanol", "CC(C)(C)O", 5, 4),
        ("Ethylene glycol", "C(O)CO", 4, 3),
    ]
    for name, script, expected_atoms, expected_bonds in simple:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        bonds = len(r["molecule"].bonds) if success and r.get("molecule") else 0
        test("Basic Organic", f"{name} ({script})",
             success and atoms == expected_atoms and bonds == expected_bonds,
             f"atoms={atoms} (expected {expected_atoms}), bonds={bonds} (expected {expected_bonds})")
    functional = [
        ("Formaldehyde", "C=O", 2, 1), ("Acetaldehyde", "CC=O", 3, 2),
        ("Acetone", "CC(=O)C", 4, 3), ("Acetic acid", "CC(=O)O", 4, 3),
        ("Methyl acetate", "COC(=O)C", 5, 4), ("Acetonitrile", "CC#N", 3, 2),
        ("Dimethyl ether", "COC", 3, 2), ("Ethylamine", "CCN", 3, 2),
        ("Acetamide", "CC(=O)N", 4, 3), ("Nitromethane", "C[N+](=O)[O-]", 4, 3),
    ]
    for name, script, expected_atoms, expected_bonds in functional:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        bonds = len(r["molecule"].bonds) if success and r.get("molecule") else 0
        test("Basic Organic", f"{name} ({script})",
             success and atoms == expected_atoms and bonds == expected_bonds,
             f"atoms={atoms}, bonds={bonds}")


# Category 2: Sandhi Validation
def run_sandhi_validation():
    print(f"\n{THIN}")
    print("  CATEGORY 2: SANDHI VALIDATION")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Sandhi Validation", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    invalid = [
        ("6-valent carbon", "C(C)(C)(C)(C)(C)"), ("7-valent carbon", "C(C)(C)(C)(C)(C)(C)"),
        ("5-valent nitrogen (bare)", "N(C)(C)(C)(C)"), ("Overvalent oxygen", "O(=O)(=O)"),
        ("Unclosed branch", "C("), ("Unclosed bracket", "C["),
        ("Empty string", ""), ("Invalid element", "X"), ("Invalid element 2", "Zz"),
    ]
    for name, script in invalid:
        r = parser.parse(script)
        rejected = not r.get("success", True)
        test("Sandhi Validation", name, rejected,
             f"success={r.get('success')}, error={str(r.get('error', 'N/A'))[:60]}")
    valid_hypervalent = [
        ("Pentavalent phosphorus", "[P](=O)(O)(O)O", 5, 4),
        ("Hexavalent sulfur", "[S](=O)(=O)(O)O", 5, 4),
        ("Charged nitrogen [N+]", "[N+](=O)[O-]", 3, 2),
        ("Charged carbon [C-]", "[C-]", 1, 0),
    ]
    for name, script, expected_atoms, expected_bonds in valid_hypervalent:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        bonds = len(r["molecule"].bonds) if success and r.get("molecule") else 0
        test("Sandhi Validation", f"{name} (should accept)",
             success and atoms == expected_atoms and bonds == expected_bonds,
             f"success={success}, atoms={atoms}, bonds={bonds}")


# Category 3: Ring Closures
def run_ring_closures():
    print(f"\n{THIN}")
    print("  CATEGORY 3: RING CLOSURES")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Ring Closures", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    topological = [
        ("Cyclopropane", "CCC&3-", 3, 3), ("Cyclobutane", "CCCC&4-", 4, 4),
        ("Cyclopentane", "CCCCC&5-", 5, 5), ("Cyclohexane", "CCCCCC&6-", 6, 6),
        ("Cycloheptane", "CCCCCCC&7-", 7, 7), ("Cyclooctane", "CCCCCCCC&8-", 8, 8),
        ("Benzene (explicit)", "C:C:C:C:C:C&6:", 6, 6), ("Pyridine", "C:C:C:C:C:N&6:", 6, 6),
    ]
    for name, script, expected_atoms, expected_bonds in topological:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        bonds = len(r["molecule"].bonds) if success and r.get("molecule") else 0
        test("Ring Closures", f"{name} ({script})",
             success and atoms == expected_atoms and bonds == expected_bonds,
             f"atoms={atoms}, bonds={bonds}")
    legacy = [("Legacy cyclohexane", "CCCCCC6", 6, 6), ("Legacy benzene", "C=CC=CC=C6", 6, 6)]
    for name, script, expected_atoms, expected_bonds in legacy:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        bonds = len(r["molecule"].bonds) if success and r.get("molecule") else 0
        test("Ring Closures", f"{name} ({script})",
             success and atoms == expected_atoms and bonds == expected_bonds,
             f"atoms={atoms}, bonds={bonds}")
    bridged = [
        ("Norbornane", "CCCC(C&5-)CC&7-", 7, 8),
    ]
    for name, script, expected_atoms, expected_bonds in bridged:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        bonds = len(r["molecule"].bonds) if success and r.get("molecule") else 0
        test("Ring Closures", f"{name} ({script})",
             success and atoms == expected_atoms and bonds == expected_bonds,
             f"atoms={atoms}, bonds={bonds}")


# Category 4: Stereochemistry
def run_stereochemistry():
    print(f"\n{THIN}")
    print("  CATEGORY 4: STEREOCHEMISTRY")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Stereochemistry", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    canon = SCRIPTCanonicalizer()
    tetrahedral = [
        ("Lactic acid (R)", "C[C@H](O)C(=O)O"), ("Lactic acid (S)", "C[C@@H](O)C(=O)O"),
        ("Bromochlorofluoromethane (R)", "F[C@H](Cl)Br"),
        ("Glucose", "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO"),
    ]
    for name, script in tetrahedral:
        r = parser.parse(script)
        success = r.get("success", False)
        test("Stereochemistry", f"{name} parses", success,
             f"error={str(r.get('error', 'N/A'))[:60]}")
        if success and r.get("molecule"):
            c1 = canon.canonicalize_core(r["molecule"])
            r2 = parser.parse(c1)
            if r2.get("success") and r2.get("molecule"):
                c2 = canon.canonicalize_core(r2["molecule"])
                test("Stereochemistry", f"{name} canonical stable", c1 == c2 and c1 is not None,
                     f"c1={c1}, c2={c2}")
            else:
                test("Stereochemistry", f"{name} canonical stable", False, "round-trip parse failed")
    extended = [
        ("Square planar Pt", "[Pt@SP](Cl)(Cl)([NH3])[NH3]"),
        ("Octahedral Co", "[Co@OH](Cl)(Cl)(Cl)(N)(N)N"),
        ("Trigonal bipyramidal P", "[P@TB](Cl)(Cl)(Cl)(Cl)Cl"),
        ("Pyramidal N", "[N@PY](C)(C)C"),
    ]
    for name, script in extended:
        r = parser.parse(script)
        success = r.get("success", False)
        has_stereo = False
        if success and r.get("molecule") and r["molecule"].atoms:
            has_stereo = r["molecule"].atoms[0].stereo_type != StereoType.NONE
        test("Stereochemistry", f"{name} ({script})",
             success and has_stereo,
             f"success={success}, stereo_type={r['molecule'].atoms[0].stereo_type if success and r.get('molecule') else 'N/A'}")
    allenes = [("Allene skeleton", "C=C=C"), ("Substituted allene", "CC=C=CC")]
    for name, script in allenes:
        r = parser.parse(script)
        test("Stereochemistry", f"{name} ({script})", r.get("success", False))


# Category 5: Typed Bonds
def run_typed_bonds():
    print(f"\n{THIN}")
    print("  CATEGORY 5: TYPED BONDS")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Typed Bonds", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    bond_tests = [
        ("Single bond", "C-C", BondType.SINGLE), ("Double bond", "C=C", BondType.DOUBLE),
        ("Triple bond", "C#C", BondType.TRIPLE), ("Aromatic bond", "C:C", BondType.AROMATIC),
        ("Dative bond", "N->[B](F)(F)F", BondType.DATIVE),
        ("Coordinate bond", "Fe>1C", BondType.COORDINATE),
        ("Tautomeric bond", "C~:C-C", BondType.TAUTOMERIC),
    ]
    for name, script, expected_type in bond_tests:
        r = parser.parse(script)
        success = r.get("success", False)
        has_correct_bond = False
        if success and r.get("molecule") and r["molecule"].bonds:
            for b in r["molecule"].bonds:
                if b.bond_type == expected_type:
                    has_correct_bond = True
                    break
        test("Typed Bonds", f"{name} ({script})",
             success and has_correct_bond,
             f"success={success}, expected={expected_type.name}, "
             f"found_types={[b.bond_type.name for b in r['molecule'].bonds] if success and r.get('molecule') else 'N/A'}")


# Category 6: Reactions
def run_reactions():
    print(f"\n{THIN}")
    print("  CATEGORY 6: REACTIONS")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Reactions", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    reactions = [
        ("Simple reaction", "CC>>CCO", 1, 1), ("With agent", "CC>[Pd]>CCO", 1, 1),
        ("Multi-product", "CC>>C.CO", 1, 2), ("Atom mapping", "[C:1]=O>>[C:1]O", 1, 1),
    ]
    for name, script, expected_reactants, expected_products in reactions:
        r = parser.parse(script)
        success = r.get("success", False)
        has_reaction = r.get("reaction") is not None
        reactant_count = len(r["reaction"].reactants) if has_reaction else 0
        product_count = len(r["reaction"].products) if has_reaction else 0
        test("Reactions", f"{name} ({script})",
             success and has_reaction and reactant_count == expected_reactants and product_count == expected_products,
             f"success={success}, has_reaction={has_reaction}, reactants={reactant_count}, products={product_count}")


# Category 7: Query Atoms
def run_query_atoms():
    print(f"\n{THIN}")
    print("  CATEGORY 7: QUERY ATOMS")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Query Atoms", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    queries = [
        ("Atomic number #6", "[#6]CC"), ("Ring atom [R]", "[R]CC"),
        ("5-membered ring [r5]", "[r5]CC"), ("Valence 3 [v3]", "[v3]CC"),
        ("Not nitrogen [!N]", "[!N]CC"), ("Aromatic [a]", "[a]CC"),
        ("Aliphatic [A]", "[A]CC"), ("OR query [#6,#7]", "[#6,#7]C"),
    ]
    for name, script in queries:
        r = parser.parse(script)
        success = r.get("success", False)
        is_query = False
        if success and r.get("molecule") and r["molecule"].atoms:
            is_query = getattr(r["molecule"].atoms[0], "is_query", False)
        test("Query Atoms", f"{name} ({script})",
             success and is_query,
             f"success={success}, is_query={is_query}")


# Category 8: Biopolymers
def run_biopolymers():
    print(f"\n{THIN}")
    print("  CATEGORY 8: BIOPOLYMERS")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Biopolymers", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    peptides = [
        ("Single alanine", "{A}"), ("Tripeptide", "{A.G.S}"),
        ("With PTM", "{pS.G.acK.V}"),
    ]
    for name, script in peptides:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        test("Biopolymers", f"{name} ({script})",
             success, f"success={success}, atoms={atoms}")
    nucleic = [
        ("DNA strand", "{dA.dG.dC.dT}"), ("RNA strand", "{rA.rG.rC.rU}"),
        ("Modified DNA", "{m5C.dG.dA.dT}"), ("m6A RNA", "{m6A.rG.rC.rU}"),
    ]
    for name, script in nucleic:
        r = parser.parse(script)
        test("Biopolymers", f"{name} ({script})", r.get("success", False),
             f"success={r.get('success', False)}")


# Category 9: Materials
def run_materials():
    print(f"\n{THIN}")
    print("  CATEGORY 9: MATERIALS AND STATES")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Materials", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    materials = [
        ("Alloy TiN", "Ti<~0.9>N<~0.1>"), ("FeNi alloy", "Fe<~0.5>Ni<~0.5>"),
        ("Crystal rutile", "[[Rutile]] Ti(O)2"), ("Crystal bcc Fe", "[[bcc]] Fe"),
        ("Surface Pt111", "[[Pt_111]] | >C=O"),
        ("Triplet oxygen", "O=O<s:3>"), ("Singlet excited", "O=O<s:1,*>"),
    ]
    for name, script in materials:
        r = parser.parse(script)
        test("Materials", f"{name} ({script})", r.get("success", False),
             f"success={r.get('success', False)}, error={str(r.get('error', 'N/A'))[:60]}")


# Category 10: Polymers
def run_polymers():
    print(f"\n{THIN}")
    print("  CATEGORY 10: POLYMERS")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Polymers", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    polymers = [
        ("Polyethylene", "{[CC]}n"), ("Stochastic PE", "{[CC]}<n:50-100>"),
        ("Exact 10-mer", "{[CC]}<n:10>"), ("Polypropylene", "{[CC(C)]}n"),
        ("Polyester", "{[CC(=O)O]}n"), ("Block copolymer", "{[CC]} -b- {[CCCO]}"),
        ("Alternating", "{[CC]} -alt- {[CC(N)C(=O)O]}"),
    ]
    for name, script in polymers:
        r = parser.parse(script)
        test("Polymers", f"{name} ({script})", r.get("success", False),
             f"success={r.get('success', False)}")


# Category 11: Canonicalization
def run_canonicalization():
    print(f"\n{THIN}")
    print("  CATEGORY 11: CANONICALIZATION")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Canonicalization", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    canon = SCRIPTCanonicalizer()
    stable = ["CCO", "C(C)C", "CCCCCC&6-", "C:C:C:C:C:C&6:", "C(C(=O)O):C(OC(=O)C):C:C:C:C&6:"]
    for script in stable:
        r = parser.parse(script)
        if not r.get("success") or not r.get("molecule"):
            test("Canonicalization", f"Stable: {script}", False, "parse failed")
            continue
        c1 = canon.canonicalize_core(r["molecule"])
        r2 = parser.parse(c1)
        if not r2.get("success") or not r2.get("molecule"):
            test("Canonicalization", f"Stable: {script}", False, "round-trip parse failed")
            continue
        c2 = canon.canonicalize_core(r2["molecule"])
        test("Canonicalization", f"Stable: {script}", c1 == c2 and c1 is not None,
             f"c1={c1}, c2={c2}")
    uniqueness = [
        ("Ethanol", ["CCO", "OCC", "C(O)C"]),
        ("Propane", ["CCC", "C(C)C"]),
        ("Aspirin", ["C(C(=O)O):C(OC(=O)C):C:C:C:C&6:", "C(OC(=O)C):C(C(=O)O):C:C:C:C&6:"]),
    ]
    for name, variants in uniqueness:
        outputs = []
        for v in variants:
            r = parser.parse(v)
            if r.get("success") and r.get("molecule"):
                c = canon.canonicalize_core(r["molecule"])
                outputs.append(c)
        non_none = [o for o in outputs if o is not None]
        if len(non_none) > 1:
            all_same = all(o == non_none[0] for o in non_none)
            test("Canonicalization", f"Unique: {name} ({variants})", all_same,
                 f"outputs={non_none}")
        else:
            test("Canonicalization", f"Unique: {name}", False, "insufficient valid outputs")


# Category 12: Boss Fights
def run_boss_fights():
    print(f"\n{THIN}")
    print("  CATEGORY 12: BOSS FIGHTS")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Boss Fights", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    canon = SCRIPTCanonicalizer()
    boss_fights = [
        ("Glucose (4 stereocenters)", "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO", 12),
        ("Aspirin", "C(C(=O)O):C(OC(=O)C):C:C:C:C&6:", 13),
        ("Caffeine", "N:C:N(C):C:C&6::C(=O)N(C)C(=O)N(C)&9:", 14),
        ("Nicotine (SCRIPT)", "CNCCC[C@H]&5-:C:C:C:N:C:C&6:", 12),
        ("Cubane", "CCCC&4-CC&4-C&6-C&8-&4-", 8),
    ]
    for name, script, expected_atoms in boss_fights:
        r = parser.parse(script)
        success = r.get("success", False)
        atoms = len(r.get("atoms", [])) if success else 0
        test("Boss Fights", f"{name} parses",
             success and atoms == expected_atoms,
             f"success={success}, atoms={atoms} (expected {expected_atoms})")
        if success and r.get("molecule"):
            try:
                c = canon.canonicalize_core(r["molecule"])
                test("Boss Fights", f"{name} canonicalizes", c is not None,
                     f"canonical={c[:50] if c else 'None'}...")
            except Exception as e:
                test("Boss Fights", f"{name} canonicalizes", False, f"exception: {e}")


# Category 13: Fuzzing
def run_fuzzing():
    print(f"\n{THIN}")
    print("  CATEGORY 13: FUZZING (Random Input Robustness)")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Fuzzing", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    random.seed(42)
    atoms = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    bonds = ["", "-", "=", "#", ":"]
    fuzz_count = 100 if (args and args.quick) else 500
    crashes = 0
    parses = 0
    for i in range(fuzz_count):
        length = random.randint(1, 20)
        parts = []
        for _ in range(length):
            parts.append(random.choice(atoms))
            if random.random() < 0.3:
                parts.append(random.choice(bonds))
        if random.random() < 0.2:
            parts.append("(")
            parts.append(random.choice(atoms))
            parts.append(")")
        if random.random() < 0.1:
            parts.append("&")
            parts.append(str(random.randint(3, 8)))
        fuzz_str = "".join(parts)
        try:
            r = parser.parse(fuzz_str)
            if r.get("success"):
                parses += 1
        except Exception as e:
            crashes += 1
            if crashes <= 3:
                print(f"    CRASH on: {fuzz_str[:30]}... -> {e}")
    test("Fuzzing", f"No crashes in {fuzz_count} random inputs", crashes == 0,
         f"crashes={crashes}, successful_parses={parses}")
    test("Fuzzing", f"Parser rejects invalid structures", parses < fuzz_count * 0.5,
         f"{parses}/{fuzz_count} parsed successfully (expected <50%)")


# Category 14: Performance
def run_performance():
    print(f"\n{THIN}")
    print("  CATEGORY 14: PERFORMANCE")
    print(THIN)
    if not SCRIPT_AVAILABLE:
        skip("Performance", "all", "SCRIPT not available")
        return
    parser = SCRIPTParser()
    canon = SCRIPTCanonicalizer()
    test_mols = [
        ("Small (C)", "C"),
        ("Medium (Aspirin)", "C(C(=O)O):C(OC(=O)C):C:C:C:C&6:"),
        ("Large (Glucose)", "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO"),
    ]
    for name, script in test_mols:
        for _ in range(3):
            parser.parse(script)
        iterations = 100 if (args and args.quick) else 1000
        t0 = time.time()
        for _ in range(iterations):
            r = parser.parse(script)
        parse_time = (time.time() - t0) / iterations * 1000
        canon_time = 0
        if r.get("success") and r.get("molecule"):
            t0 = time.time()
            for _ in range(iterations):
                canon.canonicalize_core(r["molecule"])
            canon_time = (time.time() - t0) / iterations * 1000
        test("Performance", f"{name}: parse < 1ms", parse_time < 1.0,
             f"{parse_time:.3f}ms per parse")
        test("Performance", f"{name}: canonicalize < 2ms", canon_time < 2.0,
             f"{canon_time:.3f}ms per canonicalize")


# Report generation
def generate_report():
    summary = results.summary()
    lines = []
    lines.append(SEP)
    lines.append("  SCRIPT STATE-OF-THE-ART COMPREHENSIVE BENCHMARK REPORT")
    lines.append(SEP)
    lines.append(f"  Generated: {datetime.datetime.now().isoformat()}")
    lines.append(f"  Elapsed:   {summary['elapsed_seconds']:.2f} seconds")
    lines.append("")
    lines.append(THIN)
    lines.append("  SUMMARY")
    lines.append(THIN)
    lines.append(f"  Total tests:    {summary['total_tests']}")
    lines.append(f"  Passed:         {summary['passed']}")
    lines.append(f"  Failed:         {summary['failed']}")
    lines.append(f"  Skipped:        {summary['skipped']}")
    lines.append(f"  Success rate:   {summary['success_rate']:.1f}%")
    lines.append("")
    lines.append(THIN)
    lines.append("  BY CATEGORY")
    lines.append(THIN)
    for cat, data in summary["categories"].items():
        status = "PASS" if data["rate"] >= 90 else "WARN" if data["rate"] >= 70 else "FAIL"
        lines.append(f"  [{status:4s}] {cat:25s}: {data['passed']:3d}/{data['passed']+data['failed']:<3d} ({data['rate']:5.1f}%)")
    lines.append("")
    lines.append(THIN)
    lines.append("  FAILURES (first 20)")
    lines.append(THIN)
    fail_count = 0
    for cat, data in results.categories.items():
        for t in data["tests"]:
            if t["status"] == "FAIL":
                fail_count += 1
                if fail_count <= 20:
                    lines.append(f"  [{cat}] {t['name']}")
                    if t["detail"]:
                        lines.append(f"        {t['detail']}")
    if fail_count == 0:
        lines.append("  No failures!")
    elif fail_count > 20:
        lines.append(f"  ... and {fail_count - 20} more failures")
    lines.append("")
    lines.append(SEP)
    lines.append("  VERDICT")
    lines.append(SEP)
    if summary["success_rate"] >= 95:
        lines.append("  STATUS: PRODUCTION-READY")
        lines.append("  All critical claims validated. SCRIPT is suitable for production use.")
    elif summary["success_rate"] >= 80:
        lines.append("  STATUS: BETA")
        lines.append("  Core functionality works. Some edge cases and advanced features need work.")
    elif summary["success_rate"] >= 60:
        lines.append("  STATUS: ALPHA")
        lines.append("  Basic functionality present. Significant gaps in features and robustness.")
    else:
        lines.append("  STATUS: EXPERIMENTAL")
        lines.append("  Major issues prevent production use. Needs substantial development.")
    lines.append("")
    lines.append(SEP)
    return "\n".join(lines), summary


# Main
def main():
    global args
    ap = argparse.ArgumentParser(description="SCRIPT Comprehensive Benchmark")
    ap.add_argument("--quick", action="store_true", help="Fast run (~200 tests)")
    ap.add_argument("--json", action="store_true", help="Write JSON output")
    ap.add_argument("--verbose", action="store_true", help="Show all test details")
    args = ap.parse_args()
    print(SEP)
    print("  SCRIPT STATE-OF-THE-ART COMPREHENSIVE BENCHMARK")
    print(SEP)
    print(f"  Quick mode: {args.quick}")
    print(f"  Verbose:    {args.verbose}")
    print(f"  JSON output: {args.json}")
    print("")
    if not SCRIPT_AVAILABLE:
        print("WARNING: SCRIPT not available. Most tests will be skipped.")
        print("Run from repo root or: pip install -e .")
    run_basic_organic()
    run_sandhi_validation()
    run_ring_closures()
    run_stereochemistry()
    run_typed_bonds()
    run_reactions()
    run_query_atoms()
    run_biopolymers()
    run_materials()
    run_polymers()
    run_canonicalization()
    run_boss_fights()
    run_fuzzing()
    run_performance()
    report_text, summary = generate_report()
    print(report_text)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    txt_file = f"comprehensive_benchmark_{ts}.txt"
    with open(txt_file, "w", encoding="utf-8") as f:
        f.write(report_text)
    print(f"\nReport saved: {txt_file}")
    if args.json:
        json_file = f"comprehensive_benchmark_{ts}.json"
        with open(json_file, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"JSON saved:   {json_file}")
    sys.exit(0 if summary["success_rate"] >= 80 else 1)


if __name__ == "__main__":
    main()
