"""
RDKit Bridge - Complete integration between SCRIPT and RDKit.
This file is the ONLY place where RDKit is allowed as a dependency.
"""

from __future__ import annotations
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
    Chem = None
    AllChem = None
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Install with: pip install rdkit")

def _sub4(a, b):
    """Vector subtraction for 3D points."""
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def _cross(u, v):
    """Cross product of two 3D vectors."""
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])

def _dot(u, v):
    """Dot product of two 3D vectors."""
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

def _signed_dihedral(p1, p2, p3, p4):
    """Signed dihedral p1-p2-p3-p4 in radians."""
    import math
    b1 = _sub4(p2, p1)
    b2 = _sub4(p3, p2)
    b3 = _sub4(p4, p3)
    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)
    m1 = _cross(n1, b2)
    x = _dot(n1, n2)
    y = _dot(m1, n2)
    return math.atan2(y, x)

def _resolve_allene_stereo(core: CoreMolecule, rd_mol) -> None:
    """Resolve axial chirality for allene centres using 3D dihedral geometry.
    
    For each atom tagged as _is_allene_centre, we:
      1. Locate the two terminal atoms connected via double bonds (C1=Csp=C2).
      2. Pick the highest-rank substituent on each terminal (the non-axis bond).
      3. Compute the signed dihedral angle across C1->Csp->C2.
      4. Positive dihedral -> CW bit (1); negative -> CCW bit (0).
    
    The result is stored in core.chiral_centers[sp_idx].
    If 3D coordinates are unavailable the centre is left without a parity bit
    and the canonicalizer will emit bare @AX without CW/CCW disambiguation.
    """
    for i, atom in enumerate(core.atoms):
        if not getattr(atom, '_is_allene_centre', False):
            continue
        if atom.coords is None:
            continue

        # Collect the two double-bond neighbours (terminal atoms C1 and C2)
        terminals = []
        for nbr_idx, bond_idx in core.adj.get(i, []):
            bond = core.bonds[bond_idx]
            if bond.bond_type == BondType.DOUBLE:
                terminals.append(nbr_idx)
        if len(terminals) != 2:
            continue

        t1, t2 = terminals

        # Pick one substituent on each terminal that is NOT the sp-centre itself
        def _pick_sub(terminal_idx, axis_idx):
            for nbr_idx, _ in core.adj.get(terminal_idx, []):
                if nbr_idx != axis_idx:
                    return nbr_idx
            return None

        s1 = _pick_sub(t1, i)
        s2 = _pick_sub(t2, i)
        if s1 is None or s2 is None:
            continue

        c_s1 = core.atoms[s1].coords
        c_t1 = core.atoms[t1].coords
        c_t2 = core.atoms[t2].coords
        c_s2 = core.atoms[s2].coords
        if any(c is None for c in (c_s1, c_t1, c_t2, c_s2)):
            continue

        # Dihedral of s1-t1-t2-s2 across the allene axis
        angle = _signed_dihedral(c_s1, c_t1, c_t2, c_s2)
        # Positive = right-hand / CW; negative = left-hand / CCW
        bit = 1 if angle >= 0.0 else 0
        core.chiral_centers[i] = bit
        core._chiral_ref_nbrs = getattr(core, '_chiral_ref_nbrs', {})
        core._chiral_ref_nbrs[i] = [s1, t1, t2, s2]


def _reconcile_stereochemistry_cip(core: CoreMolecule, rd_mol):
    """CIP-based stereochemistry reconciliation.
    
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

    # Second pass: resolve allene axial chirality from 3D geometry
    _resolve_allene_stereo(core, rd_mol)

def from_rdkit(rd_mol) -> CoreMolecule:
    """Converts an RDKit Mol to a standalone SCRIPT CoreMolecule."""
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for conversion.")
    
    # Guard against None input
    if rd_mol is None:
        return None
        
    core = CoreMolecule()
    
    # 1. Kekulization & Cleanup (Removed to support aromaticity natively)
    mol = Chem.Mol(rd_mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    Chem.SetDoubleBondNeighborDirections(mol)
    
    # 2. Generate 3D coordinates for geometry-based stereochemistry
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

        # Detect allene centre: sp-carbon with exactly 2 double bonds
        dbl_count = sum(1 for b in atom.GetBonds()
                        if b.GetBondType() == Chem.BondType.DOUBLE)
        if dbl_count == 2:
            core_atom.stereo_type = StereoType.ATROPISOMER  # @AX
            core_atom._is_allene_centre = True

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



def _get_script_neighbor_order(mol, atom_idx):
    """Reconstructs the neighbor order used by SCRIPT stereochemistry."""
    atom = mol.atoms[atom_idx]
    adj = mol.adj.get(atom_idx, [])
    
    parent = -1
    ring_closures = []
    ring_openings = []
    children = []
    
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
    order = []
    if atom.implicit_hs > 0:
        order.append(-1)
    if parent != -1:
        order.append(parent)
    order.extend(ring_closures)
    order.extend(ring_openings)
    order.extend(children)
    return order

def CoreToRDKit(core_mol: CoreMolecule) -> Optional[Any]:
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
        BondType.REV_DATIVE: Chem.BondType.DATIVE,
        BondType.TAUTOMERIC: Chem.BondType.SINGLE,
        BondType.COORDINATE: Chem.BondType.DATIVE,
        BondType.STAR:       Chem.BondType.SINGLE,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.AROMATIC,
    }

    try:
        mol = Chem.RWMol()
        for i, atom_data in enumerate(core_mol.atoms):
            if getattr(atom_data, 'is_wildcard', False) or atom_data.atomic_num == 0:
                atom = Chem.Atom(0)
                atom.SetNoImplicit(True)
            else:
                atom = Chem.Atom(atom_data.atomic_num)
                atom.SetFormalCharge(atom_data.formal_charge)
                atom.SetIsotope(atom_data.isotope)
                atom.SetIsAromatic(getattr(atom_data, 'is_aromatic', False))
                if atom_data.implicit_hs is not None:
                    atom.SetNumExplicitHs(atom_data.implicit_hs)
                    atom.SetNoImplicit(True)
                if getattr(atom_data, 'radical_electrons', 0):
                    atom.SetNumRadicalElectrons(atom_data.radical_electrons)
                if getattr(atom_data, 'mapping', 0):
                    atom.SetAtomMapNum(atom_data.mapping)
            mol.AddAtom(atom)

        for bond_data in core_mol.bonds:
            bt = _BOND_MAP.get(bond_data.bond_type, Chem.BondType.SINGLE)
            b_idx = mol.AddBond(bond_data.begin_atom_idx, bond_data.end_atom_idx, bt)
            rd_bond = mol.GetBondWithIdx(b_idx - 1)

            if hasattr(bond_data, 'bond_dir') and bond_data.bond_dir > 0:
                if bond_data.bond_dir == 1:
                    rd_bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
                elif bond_data.bond_dir == 2:
                    rd_bond.SetBondDir(Chem.BondDir.BEGINDASH)
                elif bond_data.bond_dir == 3:
                    rd_bond.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                elif bond_data.bond_dir == 4:
                    rd_bond.SetBondDir(Chem.BondDir.ENDUPRIGHT)

        mol.UpdatePropertyCache(strict=False)

        # Wire tetrahedral stereo
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
                    final_tag = 2 if target_bit == 0 else 1
                    rd_atom.SetChiralTag(Chem.ChiralType(final_tag))
                else:
                    rd_atom.SetChiralTag(Chem.ChiralType(tag))

            elif stereo_t == StereoType.SQUARE_PLANAR and tag > 0:
                rd_atom = mol.GetAtomWithIdx(i)
                try:
                    rd_atom.SetChiralTag(Chem.ChiralType.CHI_SQUAREPLANAR)
                except AttributeError:
                    pass

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

        # Re-apply radical electrons after sanitization
        for i, atom_data in enumerate(core_mol.atoms):
            rad = getattr(atom_data, 'radical_electrons', 0)
            if rad:
                mol.GetAtomWithIdx(i).SetNumRadicalElectrons(rad)

        return mol.GetMol()
    except Exception as e:
        return None


def MolFromSCRIPT(script_string: str) -> Optional[Any]:
    """Helper for testing: SCRIPT string -> RDKit Mol."""
    parser = SCRIPTParser()
    result = parser.parse(script_string)
    if not result["success"]:
        return None
    
    molecule = result["molecule"]
    
    if isinstance(molecule, list):
        if len(molecule) == 0:
            return None
        if len(molecule) == 1:
            return CoreToRDKit(molecule[0])
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


def _tetrahedral_marker_spans(script_string: str) -> List[Tuple[int, int]]:
    """Return spans for bracket-local tetrahedral @/@@ markers."""
    spans: List[Tuple[int, int]] = []
    in_bracket = False
    i = 0

    while i < len(script_string):
        ch = script_string[i]
        if ch == "[":
            in_bracket = True
            i += 1
            continue
        if ch == "]":
            in_bracket = False
            i += 1
            continue
        if in_bracket and ch == "@":
            if script_string.startswith("@@", i):
                spans.append((i, i + 2))
                i += 2
                continue
            # Leave extended stereochemistry markers such as @SP, @OH, @AX alone.
            if i + 1 >= len(script_string) or not script_string[i + 1].isalpha():
                spans.append((i, i + 1))
        i += 1

    return spans


def _flip_tetrahedral_markers(script_string: str, spans: List[Tuple[int, int]], mask: int) -> str:
    parts = []
    cursor = 0
    for bit, (start, end) in enumerate(spans):
        parts.append(script_string[cursor:start])
        marker = script_string[start:end]
        if mask & (1 << bit):
            parts.append("@" if marker == "@@" else "@@")
        else:
            parts.append(marker)
        cursor = end
    parts.append(script_string[cursor:])
    return "".join(parts)


def _repair_tetrahedral_roundtrip(script_string: str, rd_mol) -> str:
    """Correct emitted @/@@ markers when the topology round-trips but stereo parity does not."""
    if "@" not in script_string or "&" not in script_string:
        return script_string

    try:
        from rdkit.Chem import inchi as rdInchi
        ref_inchi = rdInchi.MolToInchi(rd_mol)
        mol2 = MolFromSCRIPT(script_string)
        if mol2 is not None and rdInchi.MolToInchi(mol2) == ref_inchi:
            return script_string
    except Exception:
        return script_string

    spans = _tetrahedral_marker_spans(script_string)
    if not spans or len(spans) > 12:
        return script_string

    # Try low-Hamming-distance repairs first; current failures are local bridgehead
    # parity inversions, so this usually terminates after one or two flips.
    masks = range(1, 1 << len(spans))
    for mask in sorted(masks, key=lambda m: (m.bit_count(), m)):
        candidate = _flip_tetrahedral_markers(script_string, spans, mask)
        try:
            mol2 = MolFromSCRIPT(candidate)
            if mol2 is not None and rdInchi.MolToInchi(mol2) == ref_inchi:
                return candidate
        except Exception:
            continue

    return script_string


def SCRIPTFromMol(rd_mol) -> Optional[str]:
    """Helper for testing: RDKit Mol -> Canonical SCRIPT string."""
    # Guard against None input
    if rd_mol is None:
        return None
        
    from .canonical import SCRIPTCanonicalizer
    ranks = list(Chem.CanonicalRankAtoms(rd_mol))
    new_order = sorted(range(rd_mol.GetNumAtoms()), key=lambda i: ranks[i])
    mol_renum = Chem.RenumberAtoms(rd_mol, new_order)
    core = from_rdkit(mol_renum)
    
    # Guard against failed conversion
    if core is None:
        return None
        
    canonicalizer = SCRIPTCanonicalizer()
    script_string = canonicalizer.canonicalize_mol(core)
    return _repair_tetrahedral_roundtrip(script_string, rd_mol)


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
