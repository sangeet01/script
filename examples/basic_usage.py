"""
Basic SCRIPT Usage Examples
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script import SCRIPTParser, parse_script

def basic_parsing_examples():
    """Demonstrate basic SCRIPT parsing"""
    print("=== Basic SCRIPT Parsing Examples ===\n")
    
    parser = SCRIPTParser()
    
    # Simple molecules
    molecules = [
        ("C", "methane"),
        ("CC", "ethane"), 
        ("CCO", "ethanol"),
        ("CC(O)C", "isopropanol"),
        ("C1CCCCC6", "cyclohexane"),
        ("C=CC=CC=C6", "benzene (Kekulé)"),
        ("[Na+].[Cl-]", "sodium chloride"),
    ]
    
    for script, name in molecules:
        print(f"Parsing {name}: {script}")
        result = parser.parse(script)
        
        if result["success"]:
            print(f"  ✓ Success! Found {len(result['atoms'])} atoms")
            for i, atom in enumerate(result["atoms"]):
                print(f"    Atom {i}: {atom['element']}")
        else:
            print(f"  ✗ Failed: {result['error']}")
        print()

def peptide_examples():
    """Demonstrate peptide mode"""
    print("=== Peptide Mode Examples ===\n")
    
    parser = SCRIPTParser()
    
    peptides = [
        ("{A}", "alanine"),
        ("{A.G.S}", "Ala-Gly-Ser tripeptide"),
        ("{C[A].G.C[A]}", "disulfide-bridged dipeptide"),
    ]
    
    for script, name in peptides:
        print(f"Parsing {name}: {script}")
        result = parser.parse(script)
        
        if result["success"]:
            print(f"  ✓ Success!")
        else:
            print(f"  ✗ Failed: {result['error']}")
        print()

def polymer_examples():
    """Demonstrate polymer mode"""
    print("=== Polymer Mode Examples ===\n")
    
    parser = SCRIPTParser()
    
    polymers = [
        ("{CC}n", "polyethylene"),
        ("{CC}100", "100-mer polyethylene"),
        ("{CC(C1=CC=CC=C6)}n", "polystyrene"),
    ]
    
    for script, name in polymers:
        print(f"Parsing {name}: {script}")
        result = parser.parse(script)
        
        if result["success"]:
            print(f"  ✓ Success!")
        else:
            print(f"  ✗ Failed: {result['error']}")
        print()

def validation_examples():
    """Demonstrate validation"""
    print("=== Validation Examples ===\n")
    
    parser = SCRIPTParser()
    
    # Valid vs invalid examples
    test_cases = [
        ("CCO", True, "valid ethanol"),
        ("C(", False, "unclosed branch"),
        ("C[", False, "unclosed bracket"),
        ("[C+]", True, "valid carbocation"),
        ("", False, "empty string"),
    ]
    
    for script, should_be_valid, description in test_cases:
        is_valid = parser.is_valid(script)
        status = "✓" if is_valid == should_be_valid else "✗"
        print(f"{status} {description}: '{script}' -> {is_valid}")

if __name__ == "__main__":
    basic_parsing_examples()
    peptide_examples() 
    polymer_examples()
    validation_examples()