import sys
sys.path.insert(0, 's:/mver33/script')

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer
from script.mol import BondType

p = SCRIPTParser()
c = SCRIPTCanonicalizer()

tests = [
    ("Ethanol",       "CCO",                                  3,  2),
    ("Nitromethane",  "C[N+](=O)[O-]",                        4,  3),
    ("Dative N->B",   "N->[B](F)(F)F",                       5,  4),
    ("Rev dative",    "[B]<-N",                               2,  1),
    ("Aspirin",       "C(C(=O)O):C(OC(=O)C):C:C:C:C&6:",    13, 13),
    ("Glucose",       "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO", 12, 12),
    ("Benzene",       "C:C:C:C:C:C&6:",                       6,  6),
    ("Perchloric",    "[Cl](=O)(=O)(=O)O",                    5,  4),
    ("N+ nitro",      "[N+](=O)[O-]",                         3,  2),

]

print("=" * 62)
print("  SCRIPT SMOKE TEST")
print("=" * 62)

all_ok = True
for name, script, ea, eb in tests:
    r = p.parse(script)
    ok = r.get("success", False)
    atoms = len(r.get("atoms", [])) if ok else 0
    bonds = len(r["molecule"].bonds) if ok and r.get("molecule") else 0
    atom_ok = atoms == ea
    bond_ok = bonds == eb
    status = "PASS" if (ok and atom_ok and bond_ok) else "FAIL"
    if status == "FAIL":
        all_ok = False
    err = r.get("error", "")[:50] if not ok else ""
    print(f"  [{status}] {name:<20} atoms={atoms}/{ea} bonds={bonds}/{eb} {err}")

print("=" * 62)
if all_ok:
    print("  All tests passed.")
else:
    print("  Some tests failed.")
