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
    """
    num_atoms = len(mol.atoms)
    if num_atoms == 0:
        return {}
    
    # 1. Initial Invariants: (AtomicNum, Degree, TotalHs, Charge, Isotope, Radical)
    invariants = []
    for i in range(num_atoms):
        atom = mol.atoms[i]
        # Calculate degree from adj list
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
    # We run enough iterations to span the diameter of the molecule
    for _ in range(num_atoms):
        new_invariants = []
        for i in range(num_atoms):
            # Collect neighbor invariants with bond types for hashing
            nbr_info = []
            for nbr_idx, bond_idx in mol.adj.get(i, []):
                bond = mol.bonds[bond_idx]
                
                # We need a stable bond value. 
                # For RDKit-compatible types: SINGLE=1, DOUBLE=2, TRIPLE=3, AROMATIC=12
                # We'll use a simple integer mapping.
                bond_val = _get_bond_order(bond.bond_type)
                nbr_info.append((invariants[nbr_idx], bond_val))
            
            # Sort neighbors to ensure determinism
            nbr_info.sort()
            
            # Combine current invariant with surroundings
            combined = (invariants[i], tuple(nbr_info))
            new_invariants.append(_stable_hash(combined))
        
        # Check if refinement stabilized
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
