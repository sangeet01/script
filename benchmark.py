import sys
import os
from rdkit import Chem
from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT

# Curated list of 100 diverse SMILES for benchmarking
# Includes chains, rings, stereo, branches, and functional groups
SMILES_DATASET = [
    # 1-10: Simple Alkanes & Alcohols
    "C", "CC", "CCC", "CCCC", "CCCCC", "CCO", "CC(C)O", "CCC(O)CC", "C(O)CO", "CC(O)C(O)C",
    # 11-20: Simple Branched
    "CC(C)C", "CC(C)(C)C", "CC(C)CC(C)C", "CC(C)C(C)C", "C(C)(C)(C)C", "CC(C)CCC(C)C", "CCC(C)C(C)CC", "CC(C)C(C)CC(C)C", "CCCC(C)C", "CCC(C)(C)C",
    # 21-30: Rings (3-8)
    "C1CC1", "C1CCC1", "C1CCCC1", "C1CCCCC1", "C1CCCCCC1", "C1CCCCCCC1", "CC1CCCCC1", "CC1(C)CCCCC1", "CCC1CCCCC1", "C1CCC(CC)CC1",
    # 31-40: Fused & Bridged Rings
    "c1ccccc1", "c1ccc2ccccc2c1", "c1ccc2c(c1)ccc3ccccc32", "C12CCCCC1CCCC2", "C1CC2CCC1CC2", "CC12CCC(CC1)CC2", "c1ccc2[nH]ccc2c1", "c1cc2ccccc2s1", "c1cc2ccccc2o1", "C1=CC2CCC1C2",
    # 41-50: Common Drugs & Bio
    "CC(=O)Oc1ccccc1C(=O)O", # Aspirin
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", # Ibuprofen
    "CC(=O)Nc1ccc(O)cc1", # Paracetamol
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine
    "CC1(C(N2C(S1)C(C2=O)NC(=O)Cc3ccccc3)C(=O)O)C", # Penicillin G
    "COc1ccc(cc1)C(C)C(=O)O", # Naproxen
    "CN(C)C(=N)N=C(N)N", # Metformin
    "CC1=C(C(=O)C2=C(C1=O)C(=C(C=C2)O)O)C", # Alizarin (complex)
    "CC1=C(C=C(C=C1)S(=O)(=O)N)N", # Sulfanilamide
    "CN1CCC[C@H]1c2cccnc2", # Nicotine
    # 51-60: Stereochemistry
    "C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O", "F[C@H](Cl)Br", "F[C@@H](Cl)Br", "C/C=C/C", "C/C=C\\C", "C1[C@H]2[C@@H]2C1", "C[C@H]1CCC[C@@H]1C", "O[C@H]1[C@@H](O)CCCC1", "C/C(=C/C)/C",
    # 61-70: Nitrogen & Oxygen groups
    "CCN", "CC(=O)N", "CC(=O)NC", "CN(C)C", "C(N)N", "C(=N)N", "CC(=N)O", "C(=O)(N)O", "CC#N", "CNC",
    # 71-80: Halogens & Sulfur
    "CCl", "CC(Cl)Cl", "C(Cl)(Cl)Cl", "C(F)(F)(F)F", "CCS", "CC(=O)S", "CC(=S)O", "CS(C)=O", "CS(=O)(=O)C", "CP(C)C",
    # 81-90: Aromatics & Heterocycles
    "c1ccncc1", "c1cncnc1", "c1ccoc1", "c1ccsc1", "c1ccc[nH]1", "c1cc2ccccc2c1", "c1ccc(cc1)c2ccccc2", "c1ccc(cc1)Oc2ccccc2", "c1ccc(cc1)Cc2ccccc2", "c1ccc(cc1)N",
    # 91-100: Complex Scaffolds
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", # Testosterone
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", # Glucose
    "CC12CC(C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)CCC34C)O", # Cortisol
    "CN1C2CCC1C(C(C2)OC(=O)c3ccccc3)C(=O)OC", # Cocaine
    "CC1(C)C2CCC1(C)C(=O)C2", # Camphor
    "CC1=C(C=C(C=C1)C)C(C)C", # p-Cymene
    "Cc1ccccc1", # Toluene
    "CC(C)c1ccc(cc1)O", # Propylphenol
    "CC1=CC(=O)CC(C1)(C)C", # Isophorone
]

def run_benchmark():
    print(f"Running SCRIPT Benchmark on {len(SMILES_DATASET)} compounds...")
    print("=" * 60)
    print(f"{'SMILES Input':<30} | {'Status':<8} | {'Details'}")
    print("-" * 60)
    
    passed = 0
    failed = 0
    
    for smi in SMILES_DATASET:
        try:
            mol_in = Chem.MolFromSmiles(smi)
            if mol_in is None:
                print(f"{smi:<30} | SKIP     | RDKit could not parse input")
                continue
                
            inchi_in = Chem.MolToInchi(mol_in)
            
            # SMILES -> SCRIPT
            script_out = SCRIPTFromMol(mol_in)
            if script_out is None:
                print(f"{smi:<30} | FAIL     | SCRIPT canonicalization failed")
                failed += 1
                continue
            
            # SCRIPT -> Mol
            mol_out = MolFromSCRIPT(script_out)
            if mol_out is None:
                print(f"{smi:<30} | FAIL     | MolFromSCRIPT failed for: {script_out}")
                failed += 1
                continue
                
            inchi_out = Chem.MolToInchi(mol_out)
            
            if inchi_in == inchi_out:
                print(f"{smi:<30} | PASS     | Round-trip successful")
                passed += 1
            else:
                print(f"{smi:<30} | FAIL     | InChI mismatch")
                print(f"  SCRIPT: {script_out}")
                print(f"  Input InChI:  {inchi_in[:50]}...")
                print(f"  Output InChI: {inchi_out[:50]}...")
                failed += 1
                
        except Exception as e:
            print(f"{smi:<30} | ERROR    | {str(e)}")
            failed += 1
            
    print("=" * 60)
    print(f"Benchmark Results:")
    print(f"  Passed: {passed}")
    print(f"  Failed: {failed}")
    if (passed + failed) > 0:
        print(f"  Success Rate: {passed/(passed+failed)*100:.1f}%")
    print("=" * 60)

if __name__ == "__main__":
    run_benchmark()
