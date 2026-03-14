"""
SCRIPT V3 - Materials Science Test Suite
Tests all 5 domains from remaining.md:
  1. Alloys & Non-Stoichiometry     (Vikalpa)
  2. Crystallography & Polymorphs   (Sthana)
  3. Extended Polymers              (Abhyasa)
  4. Electronic & Excited States    (Svara)
  5. Surface & Interface Chemistry  (Adhikarana)
"""
from script.parser import SCRIPTParser
from script.mol import CoreMolecule

def test_v3():
    parser = SCRIPTParser()
    passed = 0
    failed = 0

    def check(label, s, expect_success=True, check_fn=None):
        nonlocal passed, failed
        res = parser.parse(s)
        ok = res["success"] == expect_success
        if ok and check_fn:
            try:
                check_fn(res["molecule"])
            except AssertionError as e:
                print(f"  FAIL  [{label}] {s}")
                print(f"        Assert failed: {e}")
                failed += 1
                return
        if ok:
            print(f"  PASS  [{label}] {s}")
            passed += 1
        else:
            print(f"  FAIL  [{label}] {s}")
            print(f"        Error: {res['error']}")
            failed += 1

    print("=" * 60)
    print("V3.1 Alloys & Non-Stoichiometry (Vikalpa ~FLOAT)")
    print("=" * 60)

    def check_occupancy(mol):
        assert mol.atoms[0].occupancy == 0.9, f"Expected 0.9 got {mol.atoms[0].occupancy}"
        assert mol.atoms[1].occupancy == 0.1, f"Expected 0.1 got {mol.atoms[1].occupancy}"

    check("Alloy TiN doped", "Ti<~0.9>N<~0.1>(O)2", check_fn=check_occupancy)
    check("LiCoO2 like", "Li<~0.7>(O)2", check_fn=lambda m: None)
    check("Fractional Fe", "Fe<~0.5>Ni<~0.5>", check_fn=lambda m: None)

    print()
    print("=" * 60)
    print("V3.2 Crystallography & Polymorphs (Sthana [[label]])")
    print("=" * 60)

    def check_rutile(mol):
        assert mol.macroscopic_context == "Rutile", f"Expected Rutile got {mol.macroscopic_context}"

    def check_bcc(mol):
        assert mol.macroscopic_context == "bcc", f"Expected bcc got {mol.macroscopic_context}"

    check("Rutile TiO2", "[[Rutile]] Ti(O)2", check_fn=check_rutile)
    check("Anatase TiO2", "[[Anatase]] Ti(O)2", check_fn=lambda m: None)
    check("BCC Iron", "[[bcc]] Fe", check_fn=check_bcc)
    check("FCC Iron", "[[fcc]] Fe", check_fn=lambda m: None)
    check("Diamond", "[[diamond]] C", check_fn=lambda m: None)

    print()
    print("=" * 60)
    print("V3.3 Extended Polymers (Abhyasa)")
    print("=" * 60)

    check("Simple repeat", "{[CC]}n", check_fn=lambda m: None)
    check("Defined blocks", "{[CC]}<n:50>", check_fn=lambda m: None)

    print()
    print("=" * 60)
    print("V3.4 Electronic & Excited States (Svara s:INT, *)")
    print("=" * 60)

    def check_triplet(mol):
        last = mol.atoms[-1]
        assert last.spin == 3, f"Expected spin=3 got {last.spin}"

    def check_excited(mol):
        last = mol.atoms[-1]
        assert last.is_excited, f"Expected is_excited=True"

    def check_singlet_excited(mol):
        last = mol.atoms[-1]
        assert last.spin == 1, f"Expected spin=1 got {last.spin}"
        assert last.is_excited, f"Expected is_excited=True"

    check("Triplet O2", "O=O<s:3>", check_fn=check_triplet)
    check("Singlet excited O2", "O=O<s:1,*>", check_fn=check_singlet_excited)
    check("Excited azobenzene", "N=N<*>", check_fn=check_excited)

    print()
    print("=" * 60)
    print("V3.5 Surface & Interface Chemistry (Adhikarana |)")
    print("=" * 60)

    check("CO on Pt(111)", "[[Pt_111]] | >C=O")
    check("Li ion in LiCoO2", "[[LiCoO2]] | Li<+>")
    check("H2O on TiO2", "[[Rutile]] | O")
    check("Two-phase system", "[[Au_100]] | C=O | N")

    print()
    print("=" * 60)
    print("V3 BONUS: Reaction with atom mapping + context")
    print("=" * 60)

    check("Reaction basic", "C>>CC")
    check("Mapped reaction", "[C:1]>>[C:1]")
    check("Salt", "[Na+].[Cl-]")
    check("Haptic bond eta-5", "[Fe]*5CCCCC")
    check("Bare star bond", "C*C")

    print()
    print("=" * 60)
    print(f"TOTAL: {passed} passed, {failed} failed out of {passed+failed}")
    print("=" * 60)
    return failed == 0

if __name__ == "__main__":
    ok = test_v3()
    exit(0 if ok else 1)
