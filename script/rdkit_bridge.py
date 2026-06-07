"""
RDKit Bridge - Complete integration between SCRIPT and RDKit.
This file is the ONLY place where RDKit is allowed as a dependency.
"""

from typing import Optional, Dict, List, Tuple, Any
from .mol import CoreMolecule, CoreAtom, BondType, StereoType
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

        # Detect allene centre: atom with exactly 2 double bonds
        # Tag as @AX (axial chirality) even though RDKit loses the stereo bit —
        # at least the structural feature is recorded for downstream use.
        dbl_count = sum(1 for b in atom.GetBonds()
                        if b.GetBondType() == Chem.BondType.DOUBLE)
        if dbl_count == 2:
            core_atom.stereo_type = StereoType.ATROPISOMER  # @AX

        core.add_atom(core_atom)
        
    # Add bonds
    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        if bt == Chem.BondType.DOUBLE:
            order = BondType.DOUBLE
        elif bt == Chem.BondType.TRIPLE:
            order = BondType.TRIPLE
        elif bt == Chem.BondType.AROMATIC:
            order = BondType.AROMATIC
        elif bt == Chem.BondType.DATIVE:
            order = BondType.DATIVE
        else:
            order = BondType.SINGLE

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

    # BondType -> RDKit BondType mapping
    _BOND_MAP = {
        BondType.SINGLE:     Chem.BondType.SINGLE,
        BondType.DOUBLE:     Chem.BondType.DOUBLE,
        BondType.TRIPLE:     Chem.BondType.TRIPLE,
        BondType.AROMATIC:   Chem.BondType.AROMATIC,
        BondType.DATIVE:     Chem.BondType.DATIVE,
        BondType.REV_DATIVE: Chem.BondType.DATIVE,   # RDKit has one dative direction
        BondType.TAUTOMERIC: Chem.BondType.SINGLE,   # RDKit has no tautomeric; use single
        BondType.COORDINATE: Chem.BondType.DATIVE,
        BondType.STAR:       Chem.BondType.SINGLE,
        # Legacy int fallback
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.AROMATIC,
    }

    try:
        mol = Chem.RWMol()
        for i, atom_data in enumerate(core_mol.atoms):
            # Wildcard atom
            if getattr(atom_data, 'is_wildcard', False) or atom_data.atomic_num == 0:
                atom = Chem.Atom(0)   # atomic_num=0 = wildcard/query in RDKit
                atom.SetNoImplicit(True)
            else:
                atom = Chem.Atom(atom_data.atomic_num)
                atom.SetFormalCharge(atom_data.formal_charge)
                atom.SetIsotope(atom_data.isotope)
                atom.SetIsAromatic(getattr(atom_data, 'is_aromatic', False))
                if atom_data.implicit_hs is not None:
                    atom.SetNumExplicitHs(atom_data.implicit_hs)
                    atom.SetNoImplicit(True)
                # Radical electrons
                if getattr(atom_data, 'radical_electrons', 0):
                    atom.SetNumRadicalElectrons(atom_data.radical_electrons)
                # Atom map number (reaction tracking)
                if getattr(atom_data, 'mapping', 0):
                    atom.SetAtomMapNum(atom_data.mapping)
            mol.AddAtom(atom)

        for bond_data in core_mol.bonds:
            bt = _BOND_MAP.get(bond_data.bond_type, Chem.BondType.SINGLE)
            b_idx = mol.AddBond(bond_data.begin_atom_idx, bond_data.end_atom_idx, bt)
            rd_bond = mol.GetBondWithIdx(b_idx - 1)

            if hasattr(bond_data, 'bond_dir') and bond_data.bond_dir > 0:
                if bond_data.bond_dir == 1: rd_bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif bond_data.bond_dir == 2: rd_bond.SetBondDir(Chem.BondDir.BEGINDASH)
                elif bond_data.bond_dir == 3: rd_bond.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                elif bond_data.bond_dir == 4: rd_bond.SetBondDir(Chem.BondDir.ENDUPRIGHT)

        mol.UpdatePropertyCache(strict=False)

        # Wire tetrahedral stereo (existing logic)
        for i, atom_data in enumerate(core_mol.atoms):
            stereo_t = getattr(atom_data, 'stereo_type', StereoType.NONE)
            tag = getattr(atom_data, '_initial_tag', 0)

            if stereo_t in (StereoType.NONE, StereoType.TETRAHEDRAL) and tag > 0:
                rd_atom = mol.GetAtomWithIdx(i)
                script_order = _get_script_neighbor_order(core_mol, i)
                rdkit_order = [n.GetIdx() for n in rd_atom.GetNeighbors()]
                if atom_data.implicit_hs > 0:
                    rdkit_order.append(-1)
                if len(script_order) == 4 and len(rdkit_order) == 4:
                    stored_bit = 0 if tag == 2 else 1
                    p = permutation_parity(script_order, rdkit_order)
                    target_bit = stored_bit ^ p
                    # target_bit=0 means CCW (tag 2), target_bit=1 means CW (tag 1)
                    final_tag = 2 if target_bit == 0 else 1
                    rd_atom.SetChiralTag(Chem.ChiralType(final_tag))
                else:
                    rd_atom.SetChiralTag(Chem.ChiralType(tag))

            elif stereo_t == StereoType.SQUARE_PLANAR and tag > 0:
                # RDKit represents square planar as CHI_SQUAREPLANAR
                rd_atom = mol.GetAtomWithIdx(i)
                try:
                    rd_atom.SetChiralTag(Chem.ChiralType.CHI_SQUAREPLANAR)
                except AttributeError:
                    pass  # Older RDKit without CHI_SQUAREPLANAR

            elif stereo_t == StereoType.OCTAHEDRAL and tag > 0:
                rd_atom = mol.GetAtomWithIdx(i)
                try:
                    rd_atom.SetChiralTag(Chem.ChiralType.CHI_OCTAHEDRAL)
                except AttributeError:
                    pass

            elif stereo_t == StereoType.TRIG_BIPYRAMIDAL and tag > 0:
                rd_atom = mol.GetAtomWithIdx(i)
                try:
                    rd_atom.SetChiralTag(Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL)
                except AttributeError:
                    pass

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        Chem.SetDoubleBondNeighborDirections(mol)
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)

        # Re-apply radical electrons after sanitization (RDKit may override them)
        for i, atom_data in enumerate(core_mol.atoms):
            rad = getattr(atom_data, 'radical_electrons', 0)
            if rad:
                mol.GetAtomWithIdx(i).SetNumRadicalElectrons(rad)

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
    
    # SCRIPT Priority: H < Parent < Ring-Closures < Ring-Openings < Branches/Chain
    # Must match canonical.py ordered_neighbors build order exactly.
    order = []
    if atom.implicit_hs > 0: order.append(-1)
    if parent != -1: order.append(parent)
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
    
    molecule = result["molecule"]
    
    # Parser returns a list of CoreMolecules for multi-fragment strings
    if isinstance(molecule, list):
        if len(molecule) == 0:
            return None
        if len(molecule) == 1:
            return CoreToRDKit(molecule[0])
        # Convert each fragment and combine with RDKit
        frags = []
        for core_frag in molecule:
            frag_mol = CoreToRDKit(core_frag)
            if frag_mol is None:
                return None
            frags.append(frag_mol)
        combined = frags[0]
        for frag in frags[1:]:
            combined = Chem.CombineMols(combined, frag)
        return combined
    
    return CoreToRDKit(molecule)


def SCRIPTFromMol(rd_mol) -> Optional[str]:
    """Helper for testing: RDKit Mol -> Canonical SCRIPT string."""
    from .canonical import SCRIPTCanonicalizer
    # Normalize atom ordering to RDKit canonical ranks before DFS traversal.
    # RDKit guarantees a canonical SMILES string but NOT a canonical internal
    # atom index order — the same molecule parsed from two different SMILES
    # inputs may have different atom orderings, making SCRIPTFromMol sensitive
    # to the input SMILES rather than the molecular structure.
    # Renumbering to canonical ranks here ensures the same molecule always
    # enters from_rdkit with the same atom order.
    ranks = list(Chem.CanonicalRankAtoms(rd_mol))
    new_order = sorted(range(rd_mol.GetNumAtoms()), key=lambda i: ranks[i])
    mol_renum = Chem.RenumberAtoms(rd_mol, new_order)
    core = from_rdkit(mol_renum)
    canonicalizer = SCRIPTCanonicalizer()
    return canonicalizer.canonicalize_mol(core)


def script_to_smiles(script_string: str) -> Optional[str]:
    """Convenience: SCRIPT string -> RDKit SMILES."""
    mol = MolFromSCRIPT(script_string)
    if mol:
        return Chem.MolToSmiles(mol)
    return None


def smiles_to_script(smiles: str) -> Optional[str]:
    """Convenience: SMILES -> Canonical SCRIPT."""
    if not RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return SCRIPTFromMol(mol)
    return None