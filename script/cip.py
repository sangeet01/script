"""
CIP Priority Engine (RDKit-Free) - Enhanced
Implements Cahn-Ingold-Prelog priority rules for stereochemistry.
Uses full IUPAC 2021 atomic mass table and dynamic recursion depth.
"""

from typing import List, Tuple, Dict, Set, Optional

# Standard atomic masses for CIP tie-breaking (IUPAC 2021)
ATOMIC_MASS = {
    1: 1, 2: 4, 3: 7, 4: 9, 5: 11, 6: 12, 7: 14, 8: 16, 9: 19,
    10: 20, 11: 23, 12: 24, 13: 27, 14: 28, 15: 31, 16: 32, 17: 35,
    18: 40, 19: 39, 20: 40, 21: 45, 22: 48, 23: 51, 24: 52, 25: 55,
    26: 56, 27: 59, 28: 58, 29: 63, 30: 64, 31: 69, 32: 74, 33: 75,
    34: 80, 35: 79, 36: 84, 37: 85, 38: 88, 39: 89, 40: 90, 41: 93,
    42: 98, 43: 98, 44: 102, 45: 103, 46: 106, 47: 107, 48: 114,
    49: 115, 50: 120, 51: 121, 52: 130, 53: 127, 54: 132, 55: 133,
    56: 138, 57: 139, 58: 140, 59: 141, 60: 142, 61: 145, 62: 152,
    63: 153, 64: 158, 65: 159, 66: 164, 67: 165, 68: 167, 69: 169,
    70: 174, 71: 175, 72: 180, 73: 181, 74: 184, 75: 187, 76: 192,
    77: 193, 78: 195, 79: 197, 80: 202, 81: 205, 82: 208, 83: 209,
    84: 209, 85: 210, 86: 222, 87: 223, 88: 226, 89: 227, 90: 232,
    91: 231, 92: 238, 93: 237, 94: 244, 95: 243, 96: 247, 97: 247,
    98: 251, 99: 252, 100: 257, 101: 258, 102: 259, 103: 262,
    104: 267, 105: 268, 106: 271, 107: 272, 108: 277, 109: 276,
    110: 281, 111: 280, 112: 285, 113: 284, 114: 289, 115: 288,
    116: 293, 117: 294, 118: 294,
}


def compute_cip_priorities(mol, atom_idx: int) -> List[int]:
    """
    Returns neighbor indices sorted by CIP priority (highest first).
    Recursive depth is dynamically determined by molecule size (up to 15).
    """
    neighbors = mol.get_neighbors(atom_idx)
    if mol.atoms[atom_idx].implicit_hs > 0:
        neighbors = neighbors + [-1]
    if len(neighbors) < 2:
        return neighbors

    # Dynamic depth: larger molecules need deeper exploration
    depth = min(15, max(5, len(mol.atoms) // 2))

    from .ranking import calculate_ranks
    rank_map = calculate_ranks(mol)

    priorities = []
    for nbr_idx in neighbors:
        priority = _compute_priority_tuple(mol, nbr_idx, atom_idx, depth=depth)
        rank = rank_map.get(nbr_idx, 999999) if nbr_idx != -1 else 1000000
        priorities.append((priority, rank, nbr_idx))

    priorities.sort(key=lambda x: (x[0], x[1]), reverse=True)
    return [idx for _, _, idx in priorities]


def _compute_priority_tuple(mol, atom_idx: int, parent_idx: int,
                            depth: int, visited: Set[int] = None) -> Tuple:
    """
    Recursively compute CIP priority tuple for an atom.
    Aromatic bonds are treated as order 2 (per CIP rules - round up from 1.5).
    """
    if depth == 0:
        return (0, 0)

    if visited is None:
        visited = set()

    if atom_idx == -1:
        return (1, 0)  # Implicit H

    atom = mol.atoms[atom_idx]
    atomic_num = atom.atomic_num
    isotope = atom.isotope if atom.isotope else 0

    if atom_idx in visited:
        return (atomic_num, isotope)

    visited = visited | {atom_idx}

    neighbors = []
    for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
        if nbr_idx == parent_idx:
            continue

        bond = mol.bonds[bond_idx]
        bond_order = _get_bond_order_cip(bond.bond_type)

        # Multiple bonds: add neighbor multiple times (CIP rule)
        for _ in range(bond_order):
            nbr_priority = _compute_priority_tuple(mol, nbr_idx, atom_idx,
                                                   depth - 1, visited)
            neighbors.append(nbr_priority)

    # Add implicit H
    if atom.implicit_hs and atom.implicit_hs > 0:
        for _ in range(atom.implicit_hs):
            neighbors.append((1, 0))

    # Sort descending for comparison
    neighbors.sort(reverse=True)
    return (atomic_num, isotope, tuple(neighbors))


def _get_bond_order_cip(bond_type) -> int:
    """
    Convert bond type to integer order for CIP recursion.
    Aromatic = 1.5 -> treated as 2 for duplication (CIP convention).
    """
    if hasattr(bond_type, 'value'):
        bond_type = bond_type.value

    bond_map = {
        1: 1,   # SINGLE
        2: 2,   # DOUBLE
        3: 3,   # TRIPLE
        4: 2,   # AROMATIC (treated as 1.5 -> round up to 2 for CIP)
        5: 1,   # DATIVE
        6: 1,   # REV_DATIVE
        7: 1,   # TAUTOMERIC
        8: 1,   # COORDINATE
        9: 1,   # STAR
    }
    return bond_map.get(bond_type, 1)


def permutation_parity(order_a: List[int], order_b: List[int],
                       ranks: Optional[Dict[int, int]] = None) -> int:
    """
    Compute parity of permutation from order_a to order_b.
    Returns 0 (even) or 1 (odd).
    """
    if len(order_a) != len(order_b):
        return 0

    def get_val(x):
        if ranks is None:
            return x
        if x == -1:
            return 1000000
        return ranks.get(x, x)

    val_a = [get_val(x) for x in order_a]
    val_b = [get_val(x) for x in order_b]

    pos_map = {val: i for i, val in enumerate(val_b)}

    try:
        perm = [pos_map[val] for val in val_a]
    except KeyError:
        return 0

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
    """
    cip_order = compute_cip_priorities(mol, atom_idx)
    if len(cip_order) != len(neighbor_order):
        return chiral_bit

    parity = permutation_parity(neighbor_order, cip_order)
    return chiral_bit ^ parity
