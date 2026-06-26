"""
Standalone SCRIPT Ranking Engine
Implements a deterministic Morgan (Weisfeiler-Lehman) algorithm for atom invariants.
"""

from typing import Dict, List, Tuple
import hashlib

def calculate_ranks(mol) -> Dict[int, int]:
    """
    Assign a canonical rank to each atom in the molecule based on its environment.
    Ranks are deterministic and independent of the input atom order.
    Supports script.mol.CoreMolecule.

    For periodic molecules (mol.is_periodic = True), the lattice translation
    vector (tx, ty, tz) on each bond is included in the neighbor invariant so
    that atoms connected via different lattice directions receive distinct ranks.
    This satisfies the LQG canonicalization requirement from Krenn et al. (2022).
    """
    num_atoms = len(mol.atoms)
    if num_atoms == 0:
        return {}

    is_periodic = getattr(mol, 'is_periodic', False)

    # 1. Initial Invariants: (AtomicNum, Degree, TotalHs, Charge, Isotope, Radical)
    invariants = []
    for i in range(num_atoms):
        atom = mol.atoms[i]
        degree = len(mol.adj.get(i, []))
        inv = (
            atom.atomic_num,
            degree,
            atom.implicit_hs,
            atom.formal_charge,
            atom.isotope,
            atom.radical_electrons
        )
        invariants.append(inv)

    # 2. Iterative Refinement (Morgan/WL)
    for _ in range(num_atoms):
        new_invariants = []
        for i in range(num_atoms):
            nbr_info = []
            for nbr_idx, bond_idx in mol.adj.get(i, []):
                bond = mol.bonds[bond_idx]
                bond_val = _get_bond_order(bond.bond_type)

                if is_periodic:
                    # Include translation vector so atoms with different
                    # cross-cell connectivity get different ranks.
                    # Normalise: bond (u→v, t) and bond (v→u, -t) are the same
                    # undirected edge — always store the canonical orientation.
                    t = getattr(bond, 'translation', (0, 0, 0))
                    # Canonical form: the translation for the direction u→v
                    # where u < v by atom index, negated otherwise.
                    if bond.begin_atom_idx == i:
                        canon_t = t
                    else:
                        canon_t = (-t[0], -t[1], -t[2])
                    nbr_info.append((invariants[nbr_idx], bond_val, canon_t))
                else:
                    nbr_info.append((invariants[nbr_idx], bond_val))

            nbr_info.sort()
            combined = (invariants[i], tuple(nbr_info))
            new_invariants.append(_stable_hash(combined))

        if _ > 0 and _get_rank_order(new_invariants) == _get_rank_order(invariants):
            break
        invariants = new_invariants

    # 3. Handle Ties / Symmetry Breaking
    return _get_rank_order(invariants)

def _get_bond_order(bt) -> int:
    """Map bond types to stable integers."""
    # bt could be RDKit BondType or just an identifier
    s = str(bt)
    if 'SINGLE' in s: return 1
    if 'DOUBLE' in s: return 2
    if 'TRIPLE' in s: return 3
    if 'AROMATIC' in s: return 4
    return 0

def _stable_hash(obj) -> int:
    """Produces a stable 64-bit integer hash from a python object."""
    s = str(obj).encode('utf-8')
    h = hashlib.sha256(s).hexdigest()
    return int(h[:16], 16) 

def _get_rank_order(values: List) -> Dict[int, int]:
    """Converts a list of invariant hashes into a 0-indexed rank map."""
    indexed = [(val, i) for i, val in enumerate(values)]
    indexed.sort()
    
    ranks = {}
    current_rank = 0
    for i in range(len(indexed)):
        if i > 0 and indexed[i][0] != indexed[i-1][0]:
            current_rank = i
        ranks[indexed[i][1]] = current_rank
    return ranks
