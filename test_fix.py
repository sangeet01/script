from rdkit import Chem
from rdkit.Chem import inchi
from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT

# Test bridged benzoxazepine
smiles = 'CC1(C)Oc2ccc(C#N)cc2[C@H]2NC(=O)CO[C@@H]21'
print('Testing:', smiles)

mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print('ERROR: Cannot parse SMILES')
    exit(1)

script = SCRIPTFromMol(mol)
print('SCRIPT:', script)

if script is None:
    print('ERROR: SCRIPTFromMol returned None')
    exit(1)

mol2 = MolFromSCRIPT(script)
if mol2 is None:
    print('ERROR: MolFromSCRIPT returned None')
    exit(1)

inchi1 = inchi.MolToInchi(mol)
inchi2 = inchi.MolToInchi(mol2)

print('InChI match:', inchi1 == inchi2)
if inchi1 != inchi2:
    print('Original:', inchi1)
    print('Decoded: ', inchi2)
