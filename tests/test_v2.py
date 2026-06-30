from script.parser import SCRIPTParser
import sys

def test_v2():
    parser = SCRIPTParser()
    
    test_cases = [
        ("C<13>C", "Isotope in state block"),
        ("C<+>C", "Charge in state block"),
        ("C<h2>", "H-count in state block"),
        ("C<13,+>@2", "Multiple attributes + multiplier"),
        ("C*C", "Resonant bond (*)"),
        ("C>C", "Coordinate bond (>)"),
        ("C~C", "Ionic separator (fragment separator)"),
        ("C >> C", "Reaction support")
    ]
    
    for script, desc in test_cases:
        print(f"Testing: {desc} -> '{script}'")
        try:
            result = parser.parse(script)
            if result["success"]:
                mol = result["molecule"]
                if isinstance(mol, list):
                    print(f"  [OK] Reaction: Found {len(mol)} components")
                else:
                    print(f"  [OK] Molecule: Found {len(mol.atoms)} atoms")
            else:
                print(f"  [FAIL] Error: {result['error']}")
        except Exception as e:
            print(f"  [CRASH] {str(e)}")
        print("-" * 20)

if __name__ == "__main__":
    test_v2()
