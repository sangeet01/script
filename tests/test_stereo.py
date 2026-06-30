import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.rdkit_bridge import RDKIT_AVAILABLE, from_rdkit, SCRIPTFromMol
from script.mol import StereoType, CoreMolecule
from script.canonical import SCRIPTCanonicalizer

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem

@pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit is not available")
def test_allene_stereochemistry():
    # Build a simple chiral allene: (R)-penta-2,3-diene
    # Let's create it via SMILES and generate 3D coordinates.
    # SMILES with axial chirality: CC=[C]=CC
    mol = Chem.MolFromSmiles("CC=[C]=CC")
    assert mol is not None
    
    # Force 3D conformation
    AllChem.EmbedMolecule(mol, randomSeed=42)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Convert using from_rdkit
    core = from_rdkit(mol)
    
    # The central carbon (index 2) should be recognized as ATROPISOMER
    allene_atoms = [atom for atom in core.atoms if getattr(atom, '_is_allene_centre', False)]
    assert len(allene_atoms) == 1
    assert allene_atoms[0].stereo_type == StereoType.ATROPISOMER
    
    # Check if a bit (CW/CCW) is assigned in core.chiral_centers
    allene_idx = core.atoms.index(allene_atoms[0])
    assert allene_idx in core.chiral_centers
    assert core.chiral_centers[allene_idx] in (0, 1)

def test_periodic_adjacency():
    # Construct a periodic structure representation manually to verify DFS traversal/safety
    from script.mol import CoreAtom
    core = CoreMolecule()
    
    # Add two carbon atoms
    atom1 = CoreAtom(atomic_num=6, symbol="C")
    atom2 = CoreAtom(atomic_num=6, symbol="C")
    a1 = core.add_atom(atom1)
    a2 = core.add_atom(atom2)
    
    # Only cross-cell bond — intracell + periodic on same pair deduplicates
    core.add_bond(a1, a2, 1, translation=(1, 0, 0))
    
    assert core.is_periodic
    
    canonicalizer = SCRIPTCanonicalizer()
    res = canonicalizer.canonicalize_mol(core)
    # V3.6: LQG canonicalization — ~P~ sentinel removed, @tx,ty,tz used instead
    assert "~P~" not in res, f"Unexpected ~P~ sentinel: {res!r}"
    assert "@" in res, f"Expected @tx,ty,tz in: {res!r}"
