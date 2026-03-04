"""
CIP Priority Engine (RDKit-Free)
Implements Cahn-Ingold-Prelog priority rules for stereochemistry.
"""

from typing import List, Tuple, Dict, Set


def compute_cip_priorities(mol, atom_idx: int) -> List[int]:
    """
    Returns neighbor indices sorted by CIP priority (highest first).
    
    CIP Rules:
    1. Higher atomic number = higher priority
    2. Heavier isotope = higher priority
    3. If tied, recurse to next-nearest neighbors
    4. Multiple bonds count as multiple single bonds
    """
    neighbors = mol.get_neighbors(atom_idx)
    
    # Include implicit H if present
    if mol.atoms[atom_idx].implicit_hs > 0:
        neighbors = neighbors + [-1]
    
    if len(neighbors) < 2:
        return neighbors
    
    # Build priority tuples for each neighbor
    priorities = []
    for nbr_idx in neighbors:
        priority = _compute_priority_tuple(mol, nbr_idx, atom_idx, depth=5)
        priorities.append((priority, nbr_idx))
    
    # Sort by priority (descending), then by index for stability
    priorities.sort(key=lambda x: (x[0], -x[1] if x[1] != -1 else 999999), reverse=True)
    
    return [idx for _, idx in priorities]


def _compute_priority_tuple(mol, atom_idx: int, parent_idx: int, 
                            depth: int, visited: Set[int] = None) -> Tuple:
    """
    Recursively compute CIP priority tuple for an atom.
    Returns tuple: (atomic_num, isotope, [neighbor_priorities...])
    """
    if depth == 0:
        return (0, 0)
    
    if visited is None:
        visited = set()
    
    # Implicit hydrogen
    if atom_idx == -1:
        return (1, 0)  # H has atomic number 1
    
    atom = mol.atoms[atom_idx]
    atomic_num = atom.atomic_num
    isotope = atom.isotope if atom.isotope else 0
    
    # Avoid cycles
    if atom_idx in visited:
        return (atomic_num, isotope)
    
    visited = visited | {atom_idx}
    
    # Get neighbors (excluding parent to avoid backtracking)
    neighbors = []
    for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
        if nbr_idx == parent_idx:
            continue
        
        bond = mol.bonds[bond_idx]
        bond_order = _get_bond_order(bond.bond_type)
        
        # Multiple bonds: add neighbor multiple times
        for _ in range(bond_order):
            nbr_priority = _compute_priority_tuple(mol, nbr_idx, atom_idx, 
                                                   depth - 1, visited)
            neighbors.append(nbr_priority)
    
    # Add implicit H
    if atom.implicit_hs and atom.implicit_hs > 0:
        for _ in range(atom.implicit_hs):
            neighbors.append((1, 0))  # H priority
    
    # Sort neighbors by priority (descending)
    neighbors.sort(reverse=True)
    
    return (atomic_num, isotope, tuple(neighbors))


def _get_bond_order(bond_type) -> int:
    """Convert bond type to integer order."""
    if hasattr(bond_type, 'value'):
        bond_type = bond_type.value
    
    bond_map = {
        1: 1,   # SINGLE
        2: 2,   # DOUBLE
        3: 3,   # TRIPLE
        4: 1,   # AROMATIC (treat as 1.5, but use 1 for simplicity)
    }
    return bond_map.get(bond_type, 1)


def permutation_parity(order_a: List[int], order_b: List[int]) -> int:
    """
    Compute parity of permutation from order_a to order_b.
    Returns 0 (even) or 1 (odd).
    """
    if len(order_a) != len(order_b):
        return 0
    
    # Build mapping: order_a[i] -> position in order_b
    pos_map = {val: i for i, val in enumerate(order_b)}
    
    try:
        perm = [pos_map[val] for val in order_a]
    except KeyError:
        return 0
    
    # Count inversions
    inversions = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inversions += 1
    
    return inversions % 2


def get_cip_chirality(mol, atom_idx: int, neighbor_order: List[int], 
                      chiral_bit: int) -> int:
    """
    Transform local chirality to CIP-space chirality.
    
    Args:
        mol: CoreMolecule
        atom_idx: Chiral center index
        neighbor_order: Neighbor ordering used to compute chiral_bit
        chiral_bit: 0 (CCW) or 1 (CW) in the given neighbor_order
    
    Returns:
        CIP-space chirality: 0 or 1
    """
    cip_order = compute_cip_priorities(mol, atom_idx)
    
    if len(cip_order) != len(neighbor_order):
        return chiral_bit
    
    parity = permutation_parity(neighbor_order, cip_order)
    
    # XOR: chiral_bit XOR parity
    return chiral_bit ^ parity
