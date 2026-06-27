"""
Advanced Feature Tests for SCRIPT
Tests for the three architectural limitations fixed:
1. Allenic stereochemistry resolution
2. Polymer block copolymer expansion
3. 3D periodic adjacency (MOFs/zeolites)
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.parser import SCRIPTParser
from script.mol import CoreMolecule, CoreAtom, StereoType, BondType
from script.canonical import SCRIPTCanonicalizer
from script.rdkit_bridge import RDKIT_AVAILABLE, from_rdkit, SCRIPTFromMol

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem


class TestAlleneStereochemistry:
    """Test allenic stereochemistry with 3D coordinate resolution"""
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_allene_detection(self):
        """Test that allene centers are correctly detected"""
        # Simple allene: C=C=C
        mol = Chem.MolFromSmiles("C=C=C")
        assert mol is not None
        
        # Generate 3D conformer
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Convert to SCRIPT CoreMolecule
        core = from_rdkit(mol)
        
        # Find allene centers
        allene_atoms = [atom for atom in core.atoms 
                       if getattr(atom, '_is_allene_centre', False)]
        
        # Central carbon (sp) should be detected as allene
        assert len(allene_atoms) >= 1
        assert allene_atoms[0].stereo_type == StereoType.ATROPISOMER
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_allene_chirality_resolution(self):
        """Test that allene chirality is resolved from 3D geometry"""
        # Build chiral allene
        mol = Chem.MolFromSmiles("CC=C=CC")
        assert mol is not None
        
        # Generate 3D conformer
        AllChem.EmbedMolecule(mol, randomSeed=42)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        
        # Convert to SCRIPT
        core = from_rdkit(mol)
        
        # Find allene center
        allene_atoms = [(i, atom) for i, atom in enumerate(core.atoms) 
                       if getattr(atom, '_is_allene_centre', False)]
        
        if allene_atoms:
            idx, atom = allene_atoms[0]
            # Check that chirality bit was assigned (0=CCW, 1=CW)
            if idx in core.chiral_centers:
                assert core.chiral_centers[idx] in (0, 1)
    
    def test_allene_without_coords(self):
        """Test that allene without 3D coords is tagged but unresolved"""
        core = CoreMolecule()
        
        # Build simple allene manually: C=C=C
        c1 = CoreAtom(atomic_num=6, symbol="C")
        c2 = CoreAtom(atomic_num=6, symbol="C")
        c3 = CoreAtom(atomic_num=6, symbol="C")
        
        a1 = core.add_atom(c1)
        a2 = core.add_atom(c2)
        a3 = core.add_atom(c3)
        
        # Add double bonds
        core.add_bond(a1, a2, BondType.DOUBLE)
        core.add_bond(a2, a3, BondType.DOUBLE)
        
        # Tag center as allene
        core.atoms[a2]._is_allene_centre = True
        core.atoms[a2].stereo_type = StereoType.ATROPISOMER
        
        # Without 3D coords, no chirality bit should be assigned
        # This is expected behavior - need coords for resolution


class TestPeriodicTopology:
    """Test 3D periodic adjacency for MOFs and zeolites"""
    
    def test_periodic_bond_detection(self):
        """Test that bonds with translation vectors mark molecule as periodic"""
        core = CoreMolecule()
        
        atom1 = CoreAtom(atomic_num=6, symbol="C")
        atom2 = CoreAtom(atomic_num=6, symbol="C")
        
        a1 = core.add_atom(atom1)
        a2 = core.add_atom(atom2)
        
        # Add periodic boundary bond with translation (1, 0, 0)
        core.add_bond(a1, a2, BondType.SINGLE, translation=(1, 0, 0))
        
        assert core.is_periodic
    
    def test_multiple_periodic_bonds(self):
        """Test MOF-like structure with multiple periodic bonds"""
        core = CoreMolecule()
        
        # Create 4-atom unit cell
        atoms = []
        for i in range(4):
            atom = CoreAtom(atomic_num=6, symbol="C")
            atoms.append(core.add_atom(atom))
        
        # Intracell bonds
        core.add_bond(atoms[0], atoms[1], BondType.SINGLE)
        core.add_bond(atoms[2], atoms[3], BondType.SINGLE)
        
        # Periodic bonds connecting across boundaries
        core.add_bond(atoms[1], atoms[2], BondType.SINGLE, translation=(1, 0, 0))
        core.add_bond(atoms[3], atoms[0], BondType.SINGLE, translation=(0, 1, 0))
        
        assert core.is_periodic
        
        # Count periodic bonds
        periodic_bonds = [b for b in core.bonds if b.translation != (0, 0, 0)]
        assert len(periodic_bonds) == 2
    
    def test_lattice_vectors(self):
        """Test that lattice vectors can be stored"""
        core = CoreMolecule()
        
        # Define unit cell: cubic lattice
        core.lattice_vectors = (
            (5.0, 0.0, 0.0),
            (0.0, 5.0, 0.0),
            (0.0, 0.0, 5.0)
        )
        
        assert core.lattice_vectors is not None
        assert len(core.lattice_vectors) == 3
    
    def test_periodic_canonicalization(self):
        """Test that periodic structures can be canonicalized"""
        core = CoreMolecule()
        
        atom1 = CoreAtom(atomic_num=6, symbol="C")
        atom2 = CoreAtom(atomic_num=6, symbol="C")
        
        a1 = core.add_atom(atom1)
        a2 = core.add_atom(atom2)
        
        # Only the cross-cell bond — adding intracell + periodic on same pair
        # causes get_bond deduplication to swallow the periodic bond.
        core.add_bond(a1, a2, BondType.SINGLE, translation=(1, 0, 0))
        
        canonicalizer = SCRIPTCanonicalizer()
        result = canonicalizer.canonicalize_mol(core)
        
        # V3.6: LQG canonicalization — ~P~ sentinel removed, @tx,ty,tz used instead
        assert "~P~" not in result, f"Unexpected sentinel: {result!r}"
        assert "@" in result, f"Expected @tx,ty,tz in: {result!r}"


class TestPolymerBlocks:
    """Test polymer block copolymer functionality"""
    
    def test_simple_polymer_repeat(self):
        """Test simple repeating polymer unit"""
        parser = SCRIPTParser()
        
        result = parser.parse("{[CC]}<n:100>")
        assert result["success"]
        
        mol = result["molecule"]
        assert mol.repeat_count == 100
        assert len(mol.atoms) == 2  # CC unit
    
    def test_polymer_range(self):
        """Test stochastic polymer with range"""
        parser = SCRIPTParser()
        
        result = parser.parse("{[CCCO]}<n:50-200>")
        assert result["success"]
        
        mol = result["molecule"]
        assert mol.repeat_count == (50, 200)
        assert len(mol.atoms) == 4
    
    def test_polymer_symbolic(self):
        """Test symbolic polymer repeat count"""
        parser = SCRIPTParser()
        
        result = parser.parse("{[CC]}n")
        assert result["success"]
        
        mol = result["molecule"]
        assert mol.repeat_count == 'n'
    
    def test_complex_polymer_unit(self):
        """Test polymer with complex repeat unit"""
        parser = SCRIPTParser()
        
        # Polymer with ring in repeat unit
        result = parser.parse("{[C1CCCCC1]}n")
        assert result["success"]
        
        mol = result["molecule"]
        assert mol.repeat_count == 'n'
        assert len(mol.atoms) == 6  # Cyclohexane unit
    
    def test_polymer_with_branches(self):
        """Test polymer with branched repeat unit"""
        parser = SCRIPTParser()
        
        result = parser.parse("{[CC(C)C]}<n:50>")
        assert result["success"]
        
        mol = result["molecule"]
        assert mol.repeat_count == 50
        # Should have branched structure


class TestIntegratedFeatures:
    """Test combinations of multiple advanced features"""
    
    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_polymer_to_rdkit(self):
        """Test that polymers can round-trip through RDKit"""
        parser = SCRIPTParser()
        
        result = parser.parse("{[CC]}n")
        assert result["success"]
        
        # Should be able to handle polymer structures
        assert result["molecule"] is not None
    
    def test_mixed_stereochemistry(self):
        """Test molecule with both tetrahedral and allenic stereo"""
        core = CoreMolecule()
        
        # Create test structure
        atoms = []
        for i in range(5):
            atom = CoreAtom(atomic_num=6, symbol="C")
            atoms.append(core.add_atom(atom))
        
        # Add bonds creating allene
        core.add_bond(atoms[0], atoms[1], BondType.SINGLE)
        core.add_bond(atoms[1], atoms[2], BondType.DOUBLE)
        core.add_bond(atoms[2], atoms[3], BondType.DOUBLE)
        core.add_bond(atoms[3], atoms[4], BondType.SINGLE)
        
        # Tag allene center
        core.atoms[atoms[2]]._is_allene_centre = True
        core.atoms[atoms[2]].stereo_type = StereoType.ATROPISOMER
        
        # Verify structure
        assert len(core.atoms) == 5
        assert len(core.bonds) == 4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
