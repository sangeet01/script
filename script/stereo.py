"""
Standalone SCRIPT Stereo Engine
Implements RDKit-independent stereochemical parity management.
"""

from typing import List, Optional
try:
    from .cip import compute_cip_priorities, permutation_parity
    CIP_AVAILABLE = True
except ImportError:
    CIP_AVAILABLE = False


def get_chiral_symbol(center_idx: int, ordered_neighbors: List[int],
                      mol, ranks: List[int]) -> str:
    """
    Returns '@' or '@@' for the atom at center_idx.

    Always uses CIP-based transformation for consistency.
    """
    bit = mol.chiral_centers.get(center_idx)
    if bit is None:
        return ""

    if len(ordered_neighbors) < 4:
        return ""

    # Always use CIP transformation if available
    if CIP_AVAILABLE:
        cip_order = compute_cip_priorities(mol, center_idx)
        if len(cip_order) == 4 and len(ordered_neighbors) == 4:
            # Get reference order (CIP or DFS)
            ref_nbrs = getattr(mol, '_chiral_ref_nbrs', {}).get(center_idx, cip_order)
            
            # Transform: stored bit (in ref frame) -> CIP -> DFS
            if ref_nbrs != cip_order:
                # Stored bit is in ref frame, transform to CIP first
                parity_ref_to_cip = permutation_parity(ref_nbrs, cip_order)
                cip_bit = bit ^ parity_ref_to_cip
            else:
                # Already in CIP space
                cip_bit = bit
            
            # Transform from CIP to DFS
            parity_dfs_to_cip = permutation_parity(ordered_neighbors, cip_order)
            dfs_bit = cip_bit ^ parity_dfs_to_cip
            return "@" if dfs_bit == 0 else "@@"
    
    # Fallback: original reference frame approach
    ref_nbrs = getattr(mol, '_chiral_ref_nbrs', {}).get(center_idx)
    if ref_nbrs is None:
        return ""

    if len(ref_nbrs) != len(ordered_neighbors):
        return ""

    # Count transpositions from ref_nbrs -> ordered_neighbors
    pos_map = {}
    for i, idx in enumerate(ref_nbrs):
        pos_map[idx] = i
    try:
        perm = [pos_map[idx] for idx in ordered_neighbors]
    except KeyError:
        return ""

    swaps = 0
    p = list(perm)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                swaps += 1

    # Even swaps: same handedness as stored bit; odd: opposite
    is_ccw = (bit == 0) if (swaps % 2 == 0) else (bit == 1)
    return "@" if is_ccw else "@@"


def perceive_chirality(mol, ranks: List[int], dfs_neighbor_orders: dict = None):
    """
    Geometry-based Stereochemical Perception.

    Uses 3D coordinates to determine absolute configuration.
    Reference frame is DFS neighbor order for consistency with canonicalization.
    """
    # Don't overwrite CIP-based stereochemistry
    if not hasattr(mol, 'chiral_centers'):
        mol.chiral_centers = {}
    if not hasattr(mol, '_chiral_ref_nbrs'):
        mol._chiral_ref_nbrs = {}
    
    # Skip if molecule already has CIP-based stereochemistry
    if getattr(mol, '_cip_based_stereo', False):
        return

    for i, atom in enumerate(mol.atoms):
        nbr_indices = [n[0] for n in mol.adj.get(i, [])]
        if atom.implicit_hs > 0:
            nbr_indices.append(-1)

        if len(nbr_indices) != 4:
            continue

        # Use DFS order if provided, otherwise fall back to rank-sorted
        if dfs_neighbor_orders and i in dfs_neighbor_orders:
            ref_order = dfs_neighbor_orders[i]
        else:
            def get_rank_val(idx):
                if idx == -1:
                    return 999999
                return ranks[idx]
            ref_order = sorted(nbr_indices, key=lambda idx: (get_rank_val(idx), idx))

        # Geometry-based perception (requires 3D coordinates)
        if atom.coords is None:
            continue
            
        nbr_points = []
        for n_idx in nbr_indices:
            if n_idx == -1:
                # Implicit H: place opposite to centroid of other neighbors
                avg_x = sum(mol.atoms[j].coords[0] for j in nbr_indices if j != -1) / 3
                avg_y = sum(mol.atoms[j].coords[1] for j in nbr_indices if j != -1) / 3
                avg_z = sum(mol.atoms[j].coords[2] for j in nbr_indices if j != -1) / 3
                nbr_points.append((2*atom.coords[0]-avg_x,
                                   2*atom.coords[1]-avg_y,
                                   2*atom.coords[2]-avg_z))
            else:
                coords = mol.atoms[n_idx].coords
                if coords is None:
                    break
                nbr_points.append(coords)
        
        if len(nbr_points) != 4:
            continue

        v_map = {idx: pt for idx, pt in zip(nbr_indices, nbr_points)}
        vp1, vp2, vp3, vp4 = (v_map[ref_order[0]], v_map[ref_order[1]],
                               v_map[ref_order[2]], v_map[ref_order[3]])
        vol = _signed_volume(vp2, vp3, vp4, vp1)
        # Vol > 0: P2->P3->P4 CCW looking from P1
        # FIXED: Invert bit assignment to match RDKit convention
        bit = 1 if vol > 0 else 0
        mol.chiral_centers[i] = bit
        mol._chiral_ref_nbrs[i] = ref_order


def _signed_volume(p1, p2, p3, p4) -> float:
    """Signed volume of tetrahedron formed by 4 points."""
    v1 = (p1[0]-p4[0], p1[1]-p4[1], p1[2]-p4[2])
    v2 = (p2[0]-p4[0], p2[1]-p4[1], p2[2]-p4[2])
    v3 = (p3[0]-p4[0], p3[1]-p4[1], p3[2]-p4[2])
    return (v1[0] * (v2[1]*v3[2] - v2[2]*v3[1]) -
            v1[1] * (v2[0]*v3[2] - v2[2]*v3[0]) +
            v1[2] * (v2[0]*v3[1] - v2[1]*v3[0]))
