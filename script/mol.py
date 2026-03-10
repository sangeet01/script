"""
Core SCRIPT Graph Representation
Lightweight, RDKit-independent molecule data structure.
"""

from typing import List, Dict, Tuple, Optional, Any

class CoreAtom:
    def __init__(self, atomic_num: int, formal_charge: int = 0, 
                 isotope: int = 0, radical_electrons: int = 0,
                 symbol: str = "", is_aromatic: bool = False):
        self.atomic_num = atomic_num
        self.formal_charge = formal_charge
        self.isotope = isotope
        self.radical_electrons = radical_electrons
        self.symbol = symbol
        self.is_aromatic = is_aromatic
        self.rank = -1  # To be filled by ranking engine
        self.coords: Optional[Tuple[float, float, float]] = None
        self.implicit_hs: Optional[int] = None

        self._initial_tag = 0
        self._initial_nbrs: List[int] = []

class CoreBond:
    def __init__(self, begin_atom_idx: int, end_atom_idx: int, 
                 bond_type: Any, bond_dir: int = 0):
        self.begin_atom_idx = begin_atom_idx
        self.end_atom_idx = end_atom_idx
        self.bond_type = bond_type
        self.bond_dir = bond_dir # 0: None, 1: Up, 2: Down (relative to begin_atom)

class CoreMolecule:
    """
    RDKit-independent molecular graph.
    Used for canonicalization and parsing logic.
    """
    def __init__(self):
        self.atoms: List[CoreAtom] = []
        self.bonds: List[CoreBond] = []
        self.adj: Dict[int, List[int]] = {}  # idx -> list of (neighbor_idx, bond_idx)
        self.chiral_centers: Dict[int, int] = {} # atom_idx -> chirality_bit (0:CCW, 1:CW)

    def add_atom(self, atom: CoreAtom) -> int:
        idx = len(self.atoms)
        self.atoms.append(atom)
        self.adj[idx] = []
        return idx

    def add_bond(self, begin_idx: int, end_idx: int, bond_type: Any, bond_dir: int = 0):
        bond_idx = len(self.bonds)
        bond = CoreBond(begin_idx, end_idx, bond_type, bond_dir)
        self.bonds.append(bond)
        self.adj[begin_idx].append((end_idx, bond_idx))
        self.adj[end_idx].append((begin_idx, bond_idx))

    def get_neighbors(self, atom_idx: int) -> List[int]:
        return [nbr_idx for nbr_idx, _ in self.adj.get(atom_idx, [])]

    def get_bond(self, idx1: int, idx2: int) -> Optional[CoreBond]:
        for nbr_idx, bond_idx in self.adj.get(idx1, []):
            if nbr_idx == idx2:
                return self.bonds[bond_idx]
        return None


