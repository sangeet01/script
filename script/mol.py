"""
Core SCRIPT Graph Representation
Lightweight, RDKit-independent molecule data structure.
"""

from typing import List, Dict, Tuple, Optional, Any

class CoreAtom:
    def __init__(self, atomic_num: int, formal_charge: int = 0, 
                 isotope: int = 0, radical_electrons: int = 0,
                 symbol: str = "", is_aromatic: bool = False,
                 mapping: int = 0, occupancy: float = 1.0,
                 spin: int = 0, is_excited: bool = False):
        self.atomic_num = atomic_num
        self.formal_charge = formal_charge
        self.isotope = isotope
        self.radical_electrons = radical_electrons
        self.symbol = symbol
        self.is_aromatic = is_aromatic
        self.mapping = mapping
        self.occupancy = occupancy
        self.spin = spin
        self.is_excited = is_excited
        self.rank = -1  # To be filled by ranking engine
        self.coords: Optional[Tuple[float, float, float]] = None
        self.implicit_hs: int = 0

        self._initial_tag = 0
        self._initial_nbrs: List[int] = []
        self.chirality = 0  # 0: None, 1: CW (@@), 2: CCW (@)

class CoreBond:
    def __init__(self, begin_atom_idx: int, end_atom_idx: int, 
                 bond_type: Any, bond_dir: int = 0,
                 hapticity: int = 0, bond_class: str = ""):
        self.begin_atom_idx = begin_atom_idx
        self.end_atom_idx = end_atom_idx
        self.bond_type = bond_type
        self.bond_dir = bond_dir  # 0: None, 3: Up, 4: Down
        self.hapticity = hapticity   # eta-n for haptic organometallics (0 = not haptic)
        self.bond_class = bond_class # "dative", "rev_dative", "coordinate", "star", ""
        self.is_rc = False           # ring closure bond
        self.is_aromatic = False     # part of aromatic/resonant system

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
        self.macroscopic_context: Optional[str] = None

    def add_atom(self, atom: CoreAtom) -> int:
        idx = len(self.atoms)
        self.atoms.append(atom)
        self.adj[idx] = []
        return idx

    def add_bond(self, begin_idx: int, end_idx: int, bond_type: Any,
                 bond_dir: int = 0, hapticity: int = 0, bond_class: str = ""):
        bond_idx = len(self.bonds)
        bond = CoreBond(begin_idx, end_idx, bond_type, bond_dir,
                        hapticity=hapticity, bond_class=bond_class)
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


