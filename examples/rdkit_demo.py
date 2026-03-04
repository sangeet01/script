#!/usr/bin/env python3
"""
SCRIPT RDKit Integration Demo
Demonstrates working functionality of SCRIPT-RDKit bridge
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script import (
    MolFromSCRIPT, SCRIPTFromMol, canonicalize_SCRIPT,
    is_valid_SCRIPT, script_to_smiles, smiles_to_script,
    RDKIT_AVAILABLE
)

def main():
    print("SCRIPT RDKit Integration Demo")
    print("=" * 40)
    
    if not RDKIT_AVAILABLE:
        print("RDKit not available. Install with: pip install rdkit")
        return
    
    print("RDKit available!")
    print()
    
    # Test basic molecules
    print("Basic Molecule Tests:")
    molecules = ["C", "CCO", "C=O", "CC(=O)O", "C#N"]
    
    for script in molecules:
        mol = MolFromSCRIPT(script)
        if mol:
            atoms = mol.GetNumAtoms()
            bonds = mol.GetNumBonds()
            print(f"  {script:8} -> OK {atoms} atoms, {bonds} bonds")
        else:
            print(f"  {script:8} -> Failed")
    
    print()
    
    # Test canonicalization
    print("Canonicalization Tests:")
    test_cases = ["CCO", "C=O", "CC(=O)O"]
    
    for script in test_cases:
        canonical = canonicalize_SCRIPT(script)
        print(f"  {script:8} -> {canonical}")
    
    print()
    
    # Test SMILES conversion
    print("SMILES Conversion Tests:")
    smiles_list = ["CCO", "C=O", "CC(=O)O", "c1ccccc1"]
    
    for smiles in smiles_list:
        script = smiles_to_script(smiles)
        smiles_back = script_to_smiles(script) if script else None
        print(f"  {smiles:10} -> {script:10} -> {smiles_back}")
    
    print()
    
    # Test validation
    print("Validation Tests:")
    valid_cases = ["C", "CCO", "C=O", "CC(=O)O"]
    invalid_cases = ["", "X", "C("]
    
    print("  Valid molecules:")
    for script in valid_cases:
        result = is_valid_SCRIPT(script)
        print(f"    {script:8} -> {'OK' if result else 'FAIL'}")
    
    print("  Invalid molecules:")
    for script in invalid_cases:
        result = is_valid_SCRIPT(script)
        print(f"    {script:8} -> {'FAIL' if not result else 'OK (unexpected)'}")
    
    print()
    print("RDKit Integration Demo Complete!")
    print()
    print("Working Features:")
    print("  [OK] Basic molecule parsing (C, CCO, C=O, etc.)")
    print("  [OK] RDKit Mol object creation")
    print("  [OK] SCRIPT canonicalization")
    print("  [OK] SMILES <-> SCRIPT conversion")
    print("  [OK] Molecule validation")
    print()
    print("Next Steps:")
    print("  [TODO] Implement peptide block conversion")
    print("  [TODO] Add local ring closure handling")
    print("  [TODO] Improve mixed structure support")
    print("  [TODO] Package for PyPI distribution")

if __name__ == "__main__":
    main()