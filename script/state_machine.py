"""
Generative State Machine for SCRIPT
Implements "Sandhi" (valence guards) and "Lopa" (elision) rules
to ensure every SCRIPT string produces a physically valid molecule.
"""

from typing import List, Dict, Tuple, Optional, Any
from .mol import CoreMolecule, CoreAtom, CoreBond

# Standard valences for organic atoms (Kekule form)
DEFAULT_VALENCE = {
    "H": 1, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1,
    "P": 3, "S": 2, "Cl": 1, "Br": 1, "I": 1, "At": 1, "Ts": 1
}

# Maximum valences for hypervalent atoms (only allowed in brackets)
HYPERVALENT_MAX = {
    "P": 5, "S": 6, "Cl": 7, "Br": 7, "I": 7, "Xe": 8, "As": 5, "Se": 6
}

class GenerativeStateMachine:
    def __init__(self):
        self.mol = CoreMolecule()
        self.current_atom_idx: Optional[int] = None
        self.stack: List[int] = [] # Stack of atom indices for branching
        self.valence_used: Dict[int, int] = {} # idx -> used valence
        self.registers: Dict[str, int] = {} # Named ring registers [A] -> atom_idx
        self.is_bracket: Dict[int, bool] = {} # idx -> was in brackets

    def add_atom(self, symbol: str, charge: int = 0, isotope: int = 0, 
                 hcount: Optional[int] = None, chiral: Optional[str] = None,
                 bond_order: int = 1, bond_dir: int = 0,
                 is_bracket: bool = False) -> int:
        """Add an atom and move the state pointer to it."""
        atomic_num = self._get_atomic_num(symbol)
        atom = CoreAtom(atomic_num=atomic_num, formal_charge=charge, isotope=isotope)
        atom.symbol = symbol # Store symbol for hypervalency check
        
        atom.implicit_hs = hcount # None means automatic
        if chiral:
            # RDKit convention: CHI_TETRAHEDRAL_CW=1 (@@), CHI_TETRAHEDRAL_CCW=2 (@)
            tag = 2 if chiral == "@" else 1 if chiral == "@@" else 0
            atom._initial_tag = tag
            
        atom_idx = len(self.mol.atoms)
        self.mol.add_atom(atom)
        
        # Mark as bracket atom for hypervalency rules
        self.is_bracket[atom_idx] = is_bracket or (hcount is not None)

        
        # Initial valence used includes explicit Hydrogens from brackets
        self.valence_used[atom_idx] = hcount if hcount is not None else 0
        
        # If there's a current atom, we implicitly create a bond
        if self.current_atom_idx is not None:
            self.add_bond(self.current_atom_idx, atom_idx, bond_order, bond_dir=bond_dir)
            
        self.current_atom_idx = atom_idx
        return atom_idx

    def add_bond(self, u_idx: int, v_idx: int, order: int, bond_dir: int = 0) -> bool:
        """
        Add or upgrade a bond between u and v with 'Sandhi' valence guards.
        """
        if u_idx == v_idx: return False
        
        existing_bond = self.mol.get_bond(u_idx, v_idx)
        
        max_u = self._get_max_valence(u_idx)
        max_v = self._get_max_valence(v_idx)
        
        v_inc = order if order != 4 else 1.5
        
        if existing_bond:
            # Upgrade logic
            current_order = 1
            cur_v = 1.0
            bt = existing_bond.bond_type
            if bt == 2: current_order, cur_v = 2, 2.0
            elif bt == 3: current_order, cur_v = 3, 3.0
            elif bt == 4: current_order, cur_v = 4, 1.5
            
            # If requesting aromatic on aromatic, or lower order, just return
            if order <= current_order and order != 4: return True
            if order == 4 and current_order == 4: return True
            
            diff = v_inc - cur_v
            if diff <= 0: return True 
            
            avail_u = max_u - self.valence_used[u_idx]
            avail_v = max_v - self.valence_used[v_idx]
            
            extra = min(diff, avail_u, avail_v)
            if extra <= 0: return False 
            
            existing_bond.bond_type = order
            self.valence_used[u_idx] += extra
            self.valence_used[v_idx] += extra
            return True
            
        avail_u = max_u - self.valence_used[u_idx]
        avail_v = max_v - self.valence_used[v_idx]
        
        actual_inc = min(v_inc, avail_u, avail_v)
        if actual_inc <= 0: return False
            
        bt = 1
        if actual_inc >= 3: bt = 3
        elif actual_inc >= 2: bt = 2
        elif actual_inc >= 1.5 and order == 4: bt = 4

        self.mol.add_bond(u_idx, v_idx, bt, bond_dir=bond_dir)
        self.valence_used[u_idx] += actual_inc
        self.valence_used[v_idx] += actual_inc
        return True

    def open_branch(self):
        if self.current_atom_idx is not None:
            self.stack.append(self.current_atom_idx)

    def close_branch(self):
        if self.stack:
            self.current_atom_idx = self.stack.pop()

    def add_ring(self, identifier: Any, bond_order: int = 1):
        """Close a ring using back-counting (int) or named register (str)."""
        if self.current_atom_idx is None: return
        
        target_idx = -1
        if isinstance(identifier, int):
            # identifier is "how many atoms back"
            # 1 means this atom itself (invalid), 2 means the previous atom.
            # Index-wise: target = current - (identifier - 1)
            target_idx = self.current_atom_idx - (identifier - 1)
        elif isinstance(identifier, str):
            if identifier in self.registers:
                target_idx = self.registers[identifier]
            else:
                self.registers[identifier] = self.current_atom_idx
                return
                
        if 0 <= target_idx < len(self.mol.atoms):
            self.add_bond(self.current_atom_idx, target_idx, bond_order)

    def _get_max_valence(self, atom_idx: int) -> int:
        atom = self.mol.atoms[atom_idx]
        symbol = atom.symbol
        
        # If in brackets, allow hypervalency
        if self.is_bracket.get(atom_idx, False):
            return HYPERVALENT_MAX.get(symbol, DEFAULT_VALENCE.get(symbol, 4))
        
        return DEFAULT_VALENCE.get(symbol, 4)

    def _get_atomic_num(self, symbol: str) -> int:
        periodic_table = {
            "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
            "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
            "K": 19, "Ca": 20, "Br": 35, "I": 53, "As": 33, "Se": 34, "Xe": 54
        }
        return periodic_table.get(symbol, 6)
