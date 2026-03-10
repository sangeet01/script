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
        self.parents: Dict[int, int] = {} # idx -> parent_idx in DFS tree

    def add_atom(self, symbol: str, charge: int = 0, isotope: int = 0, 
                 hcount: Optional[int] = None, chiral: Optional[str] = None,
                 bond_order: int = 1, bond_dir: int = 0,
                 is_bracket: bool = False, is_aromatic: bool = False) -> int:
        """Add an atom and move the state pointer to it."""
        atomic_num = self._get_atomic_num(symbol)
        atom = CoreAtom(atomic_num=atomic_num, formal_charge=charge, isotope=isotope, symbol=symbol, is_aromatic=is_aromatic)
        
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
            self.parents[atom_idx] = self.current_atom_idx
            
        self.current_atom_idx = atom_idx
        return atom_idx

    def add_bond(self, u_idx: int, v_idx: int, order: int, bond_dir: int = 0) -> bool:
        """
        Add or upgrade a bond between u and v with 'Sandhi' valence guards.
        """
        if u_idx == v_idx: return False
        
        # Resolve implicit bond (-1)
        if order == -1:
            u_arom = getattr(self.mol.atoms[u_idx], 'is_aromatic', False)
            v_arom = getattr(self.mol.atoms[v_idx], 'is_aromatic', False)
            order = 4 if (u_arom and v_arom) else 1
            
        existing_bond = self.mol.get_bond(u_idx, v_idx)
        
        max_u = self._get_max_valence(u_idx)
        max_v = self._get_max_valence(v_idx)
        
        v_inc = order if order != 4 else 1.5
        
        # Aromatic bonds bypass strict fractional valence checks 
        # to allow fused rings (3 * 1.5 = 4.5 > 4 for Carbon).
        is_arom_request = (order == 4)
        
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
            
            if is_arom_request:
                extra = diff
            else:
                extra = min(diff, avail_u, avail_v)
                if extra <= 0: return False 
            
            existing_bond.bond_type = order
            self.valence_used[u_idx] += extra
            self.valence_used[v_idx] += extra
            return True
            
        avail_u = max_u - self.valence_used[u_idx]
        avail_v = max_v - self.valence_used[v_idx]
        
        if is_arom_request:
            actual_inc = v_inc
        else:
            actual_inc = min(v_inc, avail_u, avail_v)
            if actual_inc <= 0: return False
            
        bt = order
        if not is_arom_request:
            if actual_inc >= 3: bt = 3
            elif actual_inc >= 2: bt = 2
            else: bt = 1

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

    def add_ring(self, identifier: Any, bond_order: int = -1):
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
            # Mark the bond as a ring closure for the target
            bond = self.mol.get_bond(self.current_atom_idx, target_idx)
            if bond:
                bond.is_rc = True

    def add_v2_ring(self, ring_size: int, is_resonant: bool, bond_order: int = -1):
        """Close a V2 ring setting aromatic/resonant flags automatically over the topological cycle."""
        if self.current_atom_idx is None: return
        if ring_size < 3: return # Invalid ring size
        
        # Topological walk back to find the target atom
        curr_trace = self.current_atom_idx
        for _ in range(ring_size - 1):
            if curr_trace in self.parents:
                curr_trace = self.parents[curr_trace]
            else:
                # Path terminated early? Fallback to index-based if something failed, but should not happen in valid DFS
                return
        
        target_idx = curr_trace
        
        if 0 <= target_idx < len(self.mol.atoms):
            # The bond itself
            bo = 4 if (is_resonant or bond_order == 4) else (1 if bond_order == -1 else bond_order)
            self.add_bond(self.current_atom_idx, target_idx, bo)
            
            bond = self.mol.get_bond(self.current_atom_idx, target_idx)
            if bond:
                bond.is_rc = True
                bond.is_aromatic = is_resonant
                
            # If resonant, walk back on the DFS path and mark atoms and intermediate bonds as aromatic
            if is_resonant:
                # Mark current and target
                self.mol.atoms[self.current_atom_idx].is_aromatic = True
                self.mol.atoms[target_idx].is_aromatic = True
                
                # Trace back from current to target via parents
                curr = self.current_atom_idx
                while curr != target_idx and curr in self.parents:
                    p = self.parents[curr]
                    self.mol.atoms[p].is_aromatic = True
                    b = self.mol.get_bond(curr, p)
                    if b:
                        b.bond_type = 4
                        b.is_aromatic = True
                    curr = p
                    if curr == target_idx: break
                    # Safety break
                    if len(self.parents) < 1: break # Should not happen

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
