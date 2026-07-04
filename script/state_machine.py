"""
Generative State Machine for SCRIPT
Implements "Sandhi" (valence guards) and "Lopa" (elision) rules
to ensure every SCRIPT string produces a physically valid molecule.
"""

from typing import List, Dict, Tuple, Optional, Any
from .mol import CoreMolecule, CoreAtom, CoreBond, BondType, StereoType

# Standard valences for organic atoms (Kekule form)
DEFAULT_VALENCE = {
    "H": 1, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1,
    "P": 3, "S": 2, "Cl": 1, "Br": 1, "I": 1, "At": 1, "Ts": 1
}

# Maximum valences for hypervalent atoms (only allowed in brackets)
HYPERVALENT_MAX = {
    "P": 5, "S": 6, "Cl": 7, "Br": 7, "I": 7, "Xe": 8, "As": 5, "Se": 6
}

# Transition metals: variable oxidation states, use generous upper limit
# Valence guard is bypassed for these to allow complex organometallic notation
TRANSITION_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am"
}
# Generous upper valence for all transition metals
TRANSITION_METAL_MAX_VALENCE = 12

class GenerativeStateMachine:
    def __init__(self):
        self.mol = CoreMolecule()
        self.current_atom_idx: Optional[int] = None
        self.stack: List[int] = [] # Stack of atom indices for branching
        self.valence_used: Dict[int, int] = {} # idx -> used valence
        self.registers: Dict[str, int] = {} # Named ring registers [A] -> atom_idx
        self.is_bracket: Dict[int, bool] = {} # idx -> was in brackets
        self.parents: Dict[int, int] = {} # idx -> parent_idx in DFS tree
        self.component_starts: List[int] = [] # Indices where new components start
        self.valence_violations: List[Tuple[int, int, int]] = [] # (u_idx, v_idx, order) on Sandhi reject

    def new_component(self):
        """Reset the cursor to start a new disconnected component (Sandhi break)."""
        self.current_atom_idx = None
        if self.mol.atoms:
            self.component_starts.append(len(self.mol.atoms))

    def add_atom(self, symbol: str, charge: int = 0, isotope: int = 0, 
                 hcount: Optional[int] = None, chiral: Optional[str] = None,
                 bond_order: int = 1, bond_dir: int = 0,
                 is_bracket: bool = False, is_aromatic: bool = False,
                 mapping: int = 0, occupancy: float = 1.0, 
                 spin: int = 0, is_excited: bool = False,
                 bond_class: str = "", radical: int = 0,
                 translation: tuple = (0, 0, 0)) -> int:
        """Add an atom and move the state pointer to it."""
        is_wildcard = (symbol == '*')
        atomic_num = 0 if is_wildcard else self._get_atomic_num(symbol)
        atom = CoreAtom(atomic_num=atomic_num, formal_charge=charge,
                        isotope=isotope, symbol=symbol,
                        is_aromatic=is_aromatic, mapping=mapping,
                        occupancy=occupancy, spin=spin, is_excited=is_excited,
                        is_wildcard=is_wildcard)
        
        atom.implicit_hs = hcount if hcount is not None else 0
        atom.radical_electrons = radical
        if chiral:
            # Tetrahedral: CHI_TETRAHEDRAL_CW=1 (@@), CHI_TETRAHEDRAL_CCW=2 (@)
            if chiral == '@':
                atom._initial_tag = 2
                atom.stereo_type = StereoType.TETRAHEDRAL
            elif chiral == '@@':
                atom._initial_tag = 1
                atom.stereo_type = StereoType.TETRAHEDRAL
            elif chiral == '@SP':
                atom._initial_tag = 1
                atom.stereo_type = StereoType.SQUARE_PLANAR
            elif chiral == '@OH':
                atom._initial_tag = 1
                atom.stereo_type = StereoType.OCTAHEDRAL
            elif chiral == '@AX':
                atom._initial_tag = 1
                atom.stereo_type = StereoType.ATROPISOMER
            elif chiral == '@TB':
                atom._initial_tag = 1
                atom.stereo_type = StereoType.TRIG_BIPYRAMIDAL
            elif chiral == '@PY':
                atom._initial_tag = 1
                atom.stereo_type = StereoType.PYRAMIDAL
            
        atom_idx = len(self.mol.atoms)
        self.mol.add_atom(atom)
        
        # Mark as bracket atom for hypervalency rules
        self.is_bracket[atom_idx] = is_bracket or (hcount is not None)

        
        # Initial valence used includes explicit Hydrogens from brackets
        self.valence_used[atom_idx] = hcount if hcount is not None else 0
        
        # If there's a current atom, we implicitly create a bond
        if self.current_atom_idx is not None:
            ok = self.add_bond(self.current_atom_idx, atom_idx, bond_order,
                          bond_dir=bond_dir, bond_class=bond_class,
                          translation=translation)
            if not ok:
                self.valence_violations.append((self.current_atom_idx, atom_idx, bond_order))
            self.parents[atom_idx] = self.current_atom_idx
            
        self.current_atom_idx = atom_idx
        return atom_idx

    def add_bond(self, u_idx: int, v_idx: int, order: int, bond_dir: int = 0,
                 hapticity: int = 0, bond_class: str = "",
                 translation: tuple = (0, 0, 0)) -> bool:
        """
        Add or upgrade a bond between u and v with 'Sandhi' valence guards.
        bond_class maps to BondType for semantically-typed bonds.
        """
        if u_idx == v_idx: return False

        # Map bond_class to BondType — these carry semantic meaning beyond order
        _class_to_type = {
            'dative':     BondType.DATIVE,
            'rev_dative': BondType.REV_DATIVE,
            'coordinate': BondType.COORDINATE,
            'star':       BondType.STAR,
            'tautomeric': BondType.TAUTOMERIC,
        }
        # For dative/coordinate/special bonds: bypass Sandhi valence guard
        # (metal centres and donor-acceptor pairs have variable valence)
        if bond_class in _class_to_type:
            special_type = _class_to_type[bond_class]
            # Check not already bonded
            if self.mol.get_bond(u_idx, v_idx) is None:
                self.mol.add_bond(u_idx, v_idx, special_type, bond_dir=bond_dir,
                                  hapticity=hapticity, bond_class=bond_class,
                                  translation=translation)
                self.mol.atoms[u_idx]._initial_nbrs.append(v_idx)
                self.mol.atoms[v_idx]._initial_nbrs.append(u_idx)
                # Paribhasa: dative bonds consume valence only on the DONOR.
                # DATIVE (N->B): u=donor, v=acceptor. Only u's valence used.
                # REV_DATIVE (B<-N): v=donor. Only v's valence used.
                # TAUTOMERIC: shared mobile bond, count both atoms.
                # COORDINATE/STAR: metal centre, no valence increment.
                if bond_class == 'dative':
                    self.valence_used[u_idx] = self.valence_used.get(u_idx, 0) + 1
                elif bond_class == 'rev_dative':
                    self.valence_used[v_idx] = self.valence_used.get(v_idx, 0) + 1
                elif bond_class == 'tautomeric':
                    self.valence_used[u_idx] = self.valence_used.get(u_idx, 0) + 1
                    self.valence_used[v_idx] = self.valence_used.get(v_idx, 0) + 1
                # coordinate / star: no valence increment (metal centre)
            return True
        
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
        
        # When adding a non-aromatic bond to an aromatic atom, the 1.5-per-aromatic-bond
        # counting may over-restrict availability (e.g. N with 2 aromatic bonds = 3.0 used,
        # hitting its valence cap of 3, leaving no room for a substituent like cyclopropyl).
        # For non-aromatic bond requests, recalculate availability using integer aromatic counts.
        if not is_arom_request:
            def integer_valence_used(atom_idx):
                total = 0
                for bond in self.mol.bonds:
                    if bond.begin_atom_idx == atom_idx or bond.end_atom_idx == atom_idx:
                        bt = bond.bond_type
                        # Special bond types (DATIVE=5 through STAR=9) are handled
                        # by the special-bond path and must not be counted here.
                        # Aromatic bonds (bt==4) count as 1 for non-aromatic purposes.
                        if bt >= 5:
                            pass  # coordinate / dative / star: 0 covalent contribution
                        elif bt == 4:
                            total += 1
                        else:
                            total += bt
                return total
            avail_u = max_u - integer_valence_used(u_idx)
            avail_v = max_v - integer_valence_used(v_idx)
        
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

        self.mol.add_bond(u_idx, v_idx, bt, bond_dir=bond_dir,
                          hapticity=hapticity, bond_class=bond_class,
                          translation=translation)
        
        # Track "Vak Order" for chirality resolution
        self.mol.atoms[u_idx]._initial_nbrs.append(v_idx)
        self.mol.atoms[v_idx]._initial_nbrs.append(u_idx)
        
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
        """Close a ring using back-counting (int) or named register (str).
        
        For integer identifiers: identifier is "how many atoms back" in the linear string.
        - identifier=1: current atom itself (invalid/ignored)
        - identifier=2: previous atom
        - identifier=N: the atom at position (current_idx - N + 1)
        
        This matches the flat position-based encoding in canonical.py.
        """
        if self.current_atom_idx is None: return
        
        target_idx = -1
        if isinstance(identifier, int):
            # identifier is "how many atoms back" from the current position
            # Calculate target using flat linear position arithmetic
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

        # V2 &N uses the same parent-chain semantics as canonical.py:
        # &N means "connect the current atom to the atom N positions back in
        # the generated SCRIPT atom stream". Branch close/open must not shrink
        # this history, otherwise branched rings decode to the wrong atom.
        target_idx = self._get_v2_ring_target(self.current_atom_idx, ring_size)
        if target_idx is None or not (0 <= target_idx < len(self.mol.atoms)):
            return

        aromatic_path = None
        if is_resonant:
            aromatic_path = self._find_existing_path(target_idx, self.current_atom_idx)

        # The bond itself
        bo = 4 if (is_resonant or bond_order == 4) else (1 if bond_order == -1 else bond_order)
        self.add_bond(self.current_atom_idx, target_idx, bo)

        bond = self.mol.get_bond(self.current_atom_idx, target_idx)
        if bond:
            bond.is_rc = True
            bond.is_aromatic = is_resonant

        # If resonant, walk back on the DFS path and mark atoms and intermediate bonds as aromatic
        if is_resonant:
            # Mark the actual graph path that existed before the closure.
            # The emitted atom interval may contain branches that are not in
            # this aromatic ring.
            path = aromatic_path or [target_idx, self.current_atom_idx]
            for idx in path:
                self.mol.atoms[idx].is_aromatic = True

            for u, v in zip(path, path[1:]):
                b = self.mol.get_bond(u, v)
                if b:
                    b.bond_type = 4
                    b.is_aromatic = True

    def _get_v2_ring_target(self, atom_idx: int, ring_size: int) -> Optional[int]:
        """Find the V2 ring closure target using the flat emitted atom stream.

        &N ring closures refer to the N-atom cycle in the generated SCRIPT
        atom stream, including branched atoms. Using the parent chain can
        miss atoms that were emitted in branches, so we calculate the target
        by flat atom index arithmetic.
        """
        if atom_idx is None or ring_size < 2:
            return None

        # Prefer flat emitted-stream arithmetic (canonicalizer uses this),
        # but only accept it if the current graph already contains a path
        # of the requested topological length between the two atoms. This
        # lets canonical strings roundtrip while preserving parent-chain
        # semantics for user-written branched encodings.
        flat_target = atom_idx - (ring_size - 1)
        if 0 <= flat_target < len(self.mol.atoms):
            path = self._find_existing_path(flat_target, atom_idx)
            if path is not None:
                # If the flat target equals the parent-chain target, prefer
                # the flat target (common for canonical roundtrips). If it
                # differs, don't greedily choose flat here because some
                # user-written branched scripts expect parent-chain semantics.
                parent = atom_idx
                for _ in range(ring_size - 1):
                    parent = self.parents.get(parent)
                    if parent is None:
                        break
                # If the parent chain is absent (None) or matches the flat
                # target, prefer flat. Absence of a parent chain commonly
                # occurs when parsing canonical output where the emission
                # ordering differs slightly; accepting flat in that case
                # preserves canonical roundtrips.
                if parent is None or parent == flat_target:
                    return flat_target

        # Fallback: walk parent chain (DFS emission history) and accept it
        # only if it corresponds to a path of the requested size.
        target = atom_idx
        for _ in range(ring_size - 1):
            target = self.parents.get(target)
            if target is None:
                return None
        path = self._find_existing_path(target, atom_idx)
        if path is not None and len(path) == ring_size:
            # Debug hook:
            # print(f"_get_v2_ring_target: using parent_target {target} for atom {atom_idx} ring_size {ring_size}")
            return target

        return None

    def _find_existing_path(self, start_idx: int, end_idx: int) -> Optional[List[int]]:
        """Find a path in the current graph before a ring closure bond is added."""
        if start_idx == end_idx:
            return [start_idx]

        queue: List[Tuple[int, List[int]]] = [(start_idx, [start_idx])]
        seen = {start_idx}

        while queue:
            atom_idx, path = queue.pop(0)
            for nbr_idx, _ in self.mol.adj.get(atom_idx, []):
                if nbr_idx in seen:
                    continue
                next_path = path + [nbr_idx]
                if nbr_idx == end_idx:
                    return next_path
                seen.add(nbr_idx)
                queue.append((nbr_idx, next_path))

        return None

    def finalize_valences(self):
        """
        Sandhi Finalization: Calculates remaining valences and fills them with
        implicit Hydrogens (Lopa). This permits bare symbols in output.
        """
        for i, atom in enumerate(self.mol.atoms):
            # Only fill for atoms that don't have explicit H specifying
            # and are in the organic subset.
            if self.is_bracket.get(i, False):
                continue
            
            symbol = atom.symbol
            if symbol not in DEFAULT_VALENCE:
                continue
                
            max_v = DEFAULT_VALENCE[symbol]
            used = self.valence_used.get(i, 0)
            
            # If used < max, fill the rest with implicit H
            if used < max_v:
                atom.implicit_hs = int(max_v - used)
                self.valence_used[i] = max_v

    def _get_max_valence(self, atom_idx: int) -> int:
        atom = self.mol.atoms[atom_idx]
        symbol = atom.symbol
        
        # Transition metals bypass strict valence guards
        if symbol in TRANSITION_METALS:
            return TRANSITION_METAL_MAX_VALENCE
        
        base = DEFAULT_VALENCE.get(symbol, 4)
        
        # Formal charge shifts available valence:
        # [N+] has valence 4, [N-] has valence 2, [O+] has valence 3, etc.
        charge = getattr(atom, 'formal_charge', 0)
        if charge != 0:
            charged_val = base + charge
            # If in brackets, also allow hypervalency on top of charge
            if self.is_bracket.get(atom_idx, False):
                hyper = HYPERVALENT_MAX.get(symbol, base)
                return max(charged_val, hyper)
            return max(charged_val, 1)
        
        # If in brackets, allow hypervalency
        if self.is_bracket.get(atom_idx, False):
            return HYPERVALENT_MAX.get(symbol, base)
        
        return base

    def _get_atomic_num(self, symbol: str) -> int:
        periodic_table = {
            "H": 1, "He": 2,
            "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
            "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
            "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
            "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32,
            "As": 33, "Se": 34, "Br": 35, "Kr": 36,
            "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
            "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
            "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
            "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61,
            "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68,
            "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75,
            "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
            "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
            "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93,
            "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
            "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106,
            "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
            "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
        }
        return periodic_table.get(symbol, 0 if symbol == "*" else 6)
