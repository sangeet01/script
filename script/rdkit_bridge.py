"""
RDKit Bridge - Complete integration between SCRIPT and RDKit.
This file is the ONLY place where RDKit is allowed as a dependency.
"""

from typing import Optional, Dict, List, Tuple, Any
from .mol import CoreMolecule, CoreAtom
from .parser import SCRIPTParser
from .cip import compute_cip_priorities, permutation_parity, get_cip_chirality

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Install with: pip install rdkit")

def from_rdkit(rd_mol) -> CoreMolecule:
    """Converts an RDKit Mol to a standalone SCRIPT CoreMolecule."""
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for conversion.")
        
    core = CoreMolecule()
    
    # 1. Kekulization & Cleanup (Removed to support aromaticity natively)
    mol = Chem.Mol(rd_mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    Chem.SetDoubleBondNeighborDirections(mol)
    
    # 2. Generate 3D coordinates for geometry-based stereochemistry
    # This ensures SCRIPT's DFS-based perception works correctly
    has_stereo = any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms())
    if has_stereo and mol.GetNumConformers() == 0:
        try:
            # Try multiple times with different random seeds for difficult molecules
            for seed in [42, 123, 456, 789, 2023]:
                AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=False, maxAttempts=100)
                if mol.GetNumConformers() > 0:
                    break
            # If 3D fails, use 2D
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)
        except:
            pass
    
    # Add atoms
    for i, atom in enumerate(mol.GetAtoms()):
        core_atom = CoreAtom(
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            isotope=atom.GetIsotope(),
            radical_electrons=atom.GetNumRadicalElectrons(),
            symbol=atom.GetSymbol(),
            is_aromatic=atom.GetIsAromatic()
        )
        core_atom.implicit_hs = atom.GetTotalNumHs()
        
        # Capture 3D coordinates for geometry-based stereochemistry
        if mol.GetNumConformers() > 0:
            pos = mol.GetConformer().GetAtomPosition(i)
            core_atom.coords = (pos.x, pos.y, pos.z)
        
        # Store RDKit chiral tag for CIP reconciliation
        chiral_tag = atom.GetChiralTag()
        if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            core_atom._rdkit_chiral_tag = int(chiral_tag)
        else:
            core_atom._rdkit_chiral_tag = 0
            
        core.add_atom(core_atom)
        
    # Add bonds
    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        order = 1
        if bt == Chem.BondType.DOUBLE: order = 2
        elif bt == Chem.BondType.TRIPLE: order = 3
        elif bt == Chem.BondType.AROMATIC: order = 4
        
        core.add_bond(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            order,
            int(bond.GetBondDir())
        )
    
    # CIP-based stereochemistry reconciliation
    _reconcile_stereochemistry_cip(core, mol)
        
    return core

def _reconcile_stereochemistry_cip(core: CoreMolecule, rd_mol):
    """
    CIP-based stereochemistry reconciliation.
    Transforms RDKit's chiral tags to SCRIPT's geometry-based representation
    using CIP priorities as universal reference frame.
    """
    has_any_stereo = False
    
    for i, atom in enumerate(core.atoms):
        if not hasattr(atom, '_rdkit_chiral_tag') or atom._rdkit_chiral_tag == 0:
            continue
        
        has_any_stereo = True
        
        # Get RDKit neighbor order
        rd_atom = rd_mol.GetAtomWithIdx(i)
        rdkit_neighbors = [n.GetIdx() for n in rd_atom.GetNeighbors()]
        if atom.implicit_hs > 0:
            rdkit_neighbors.append(-1)
        
        if len(rdkit_neighbors) != 4:
            continue
        
        # Compute CIP priorities using CoreMolecule (graph-based, not RDKit order)
        cip_order = compute_cip_priorities(core, i)
        
        if len(cip_order) != 4:
            continue
        
        # Transform RDKit tag to CIP space
        # RDKit: 1=CW(@@), 2=CCW(@)
        rdkit_bit = 0 if atom._rdkit_chiral_tag == 2 else 1
        parity_rdkit = permutation_parity(rdkit_neighbors, cip_order)
        cip_chirality = rdkit_bit ^ parity_rdkit
        
        # Store CIP-space chirality and CIP order as reference
        # This will be transformed to DFS space during canonicalization
        core.chiral_centers[i] = cip_chirality
        core._chiral_ref_nbrs = getattr(core, '_chiral_ref_nbrs', {})
        core._chiral_ref_nbrs[i] = cip_order
    
    # Mark that this molecule has CIP-based stereochemistry
    if has_any_stereo:
        core._cip_based_stereo = True

def CoreToRDKit(core_mol: CoreMolecule) -> Optional[Chem.Mol]:
    """Converts a SCRIPT CoreMolecule to an RDKit Mol."""
    if not RDKIT_AVAILABLE:
        return None
        
    try:
        mol = Chem.RWMol()
        for i, atom_data in enumerate(core_mol.atoms):
            atom = Chem.Atom(atom_data.atomic_num)
            atom.SetFormalCharge(atom_data.formal_charge)
            atom.SetIsotope(atom_data.isotope)
            atom.SetIsAromatic(getattr(atom_data, 'is_aromatic', False))
            
            if atom_data.implicit_hs is not None:
                atom.SetNumExplicitHs(atom_data.implicit_hs)
                atom.SetNoImplicit(True)
            
            mol.AddAtom(atom)

        for bond_data in core_mol.bonds:
            order = bond_data.bond_type
            bt = Chem.BondType.SINGLE
            if order == 2: bt = Chem.BondType.DOUBLE
            elif order == 3: bt = Chem.BondType.TRIPLE
            elif order == 4: bt = Chem.BondType.AROMATIC
            
            b_idx = mol.AddBond(bond_data.begin_atom_idx, bond_data.end_atom_idx, bt)
            rd_bond = mol.GetBondWithIdx(b_idx - 1)
            
            if hasattr(bond_data, 'bond_dir') and bond_data.bond_dir > 0:
                if bond_data.bond_dir == 1: rd_bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif bond_data.bond_dir == 2: rd_bond.SetBondDir(Chem.BondDir.BEGINDASH)
                elif bond_data.bond_dir == 3: rd_bond.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                elif bond_data.bond_dir == 4: rd_bond.SetBondDir(Chem.BondDir.ENDUPRIGHT)

        mol.UpdatePropertyCache(strict=False)
        
        # FIX: Correct chiral tags based on neighbor order reconciliation
        for i, atom_data in enumerate(core_mol.atoms):
            if hasattr(atom_data, '_initial_tag') and atom_data._initial_tag > 0:
                rd_atom = mol.GetAtomWithIdx(i)
                # SCRIPT priority order
                script_order = _get_script_neighbor_order(core_mol, i)
                # Actual bond addition order in RDKit
                rdkit_order = [n.GetIdx() for n in rd_atom.GetNeighbors()]
                if atom_data.implicit_hs > 0:
                    rdkit_order.append(-1)
                
                if len(script_order) == 4 and len(rdkit_order) == 4:
                    # Stored bit: 0=@ (CCW), 1=@@ (CW)
                    stored_bit = 0 if atom_data._initial_tag == 2 else 1
                    # Calculate parity of permutation from SCRIPT to RDKit
                    p = permutation_parity(script_order, rdkit_order)
                    target_bit = stored_bit ^ p
                    tag = 2 if target_bit == 0 else 1
                    rd_atom.SetChiralTag(Chem.ChiralType(tag))
                else:
                    # Fallback
                    rd_atom.SetChiralTag(Chem.ChiralType(atom_data._initial_tag))

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        Chem.SetDoubleBondNeighborDirections(mol)
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return mol.GetMol()
    except Exception as e:
        # print(f"DEBUG: CoreToRDKit failed: {e}")
        return None

def _get_script_neighbor_order(mol, atom_idx):
    """Reconstructs the neighbor order used by SCRIPT stereochemistry."""
    atom = mol.atoms[atom_idx]
    adj = mol.adj.get(atom_idx, [])
    
    parent = -1
    ring_closures = [] # Back-counts initiated BY this atom
    ring_openings = []  # Back-counts pointing TO this atom
    children = [] # Normal tree bonds to descendants
    
    for nbr_idx, bond_idx in adj:
        bond = mol.bonds[bond_idx]
        if bond.end_atom_idx == atom_idx:
            if getattr(bond, 'is_rc', False):
                ring_openings.append(nbr_idx)
            else:
                parent = nbr_idx
        else:
            if getattr(bond, 'is_rc', False):
                ring_closures.append(nbr_idx)
            else:
                children.append(nbr_idx)
    
    # SCRIPT Priority: Parent < H < Ring-Closures < Ring-Openings < Branches/Chain
    order = []
    if parent != -1: order.append(parent)
    if atom.implicit_hs > 0: order.append(-1)
    order.extend(ring_closures)
    order.extend(ring_openings)
    order.extend(children)
    return order

def MolFromSCRIPT(script_string: str) -> Optional[Chem.Mol]:
    """Helper for testing: SCRIPT string -> RDKit Mol."""
    parser = SCRIPTParser()
    result = parser.parse(script_string)
    if not result["success"]:
        return None
    return CoreToRDKit(result["molecule"])

def SCRIPTFromMol(rd_mol) -> Optional[str]:
    """Helper for testing: RDKit Mol -> Canonical SCRIPT string."""
    from .canonical import SCRIPTCanonicalizer
    core = from_rdkit(rd_mol)
    canonicalizer = SCRIPTCanonicalizer()
    return canonicalizer.canonicalize_mol(core)