"""
SCRIPT Chirality Resolver - Paninian Model
==========================================
Implements native, RDKit-free chirality resolution using the
Ashtadhyayi-inspired rule system:

  Pratyaya   = @/@@  marker stored on the atom
  Vak Order  = DFS neighbor order (the order neighbors were "spoken")
  Sandhi     = CIP-like priority assignment (Sandhi rules)
  Lopa       = implicit H (elided, always lowest priority)

The resolver runs as a post-pass after the molecule is fully built,
just as Sanskrit Sandhi rules apply after morphological construction.
"""

from typing import List, Dict, Optional, Tuple, Any
from .mol import CoreMolecule, CoreAtom, StereoType

# Chirality bit convention (matches RDKit encoding for bridge compatibility)
CHI_NONE = 0
CHI_CCW  = 2   # @ = counterclockwise
CHI_CW   = 1   # @@ = clockwise


# ---------------------------------------------------------------------------
# CIP Priority Rules - Sandhi (junction rules)
# Priority tuple: (atomic_num, mass, -charge, ...)
# Higher tuple = higher CIP priority
# ---------------------------------------------------------------------------

def _priority_key(atom: CoreAtom) -> Tuple:
    """
    Sandhi Rule: returns the primary CIP priority key for an atom.
    Apavada (exception): isotope mass overrides atomic number for tie-breaking.
    Utsarga (general): atomic number is the base rule.
    """
    # Higher atomic num = higher priority (Utsarga)
    # Higher isotope mass = higher priority (Apavada)
    # More positive charge = higher priority
    mass = atom.isotope if atom.isotope > 0 else _standard_mass(atom.atomic_num)
    return (atom.atomic_num, mass, -atom.formal_charge)


def _standard_mass(atomic_num: int) -> int:
    """Approximate standard atomic mass for CIP tie-breaking."""
    # Only need coarse values; exact isotopes use atom.isotope directly
    MASS = {
        1: 1, 5: 11, 6: 12, 7: 14, 8: 16, 9: 19, 11: 23, 12: 24,
        13: 27, 14: 28, 15: 31, 16: 32, 17: 35, 19: 39, 20: 40,
        26: 56, 27: 59, 28: 58, 29: 64, 30: 65, 35: 80, 53: 127,
    }
    return MASS.get(atomic_num, atomic_num * 2)


# ---------------------------------------------------------------------------
# ChiralResolver
# ---------------------------------------------------------------------------

class ChiralResolver:
    """
    Post-pass chirality resolver.
    
    Usage:
        resolver = ChiralResolver(mol)
        resolver.resolve()
        # mol.atoms[i].chirality is now CHI_CW, CHI_CCW, or CHI_NONE
    """

    def __init__(self, mol: CoreMolecule):
        self.mol = mol

    def resolve(self):
        """
        Pratyaya Rule: for each atom that carries an _initial_tag,
        resolve the @/@@  sense into a definitive chirality bit.
        """
        # Compute the Morgan-WL rank map once for the whole molecule and
        # reuse it for every chiral atom. Without this, calculate_ranks was
        # being invoked twice per chiral atom (once here, once inside
        # compute_cip_priorities), making glucose spend ~75% of parse time
        # in ranking.py. The rank map depends only on the graph, which is
        # frozen by the time resolve() runs.
        from .ranking import calculate_ranks
        rank_map = calculate_ranks(self.mol)

        for idx, atom in enumerate(self.mol.atoms):
            tag = getattr(atom, '_initial_tag', 0)
            if tag == 0:
                atom.chirality = CHI_NONE
                continue

            # Sthiti form (CIP-absolute): @R/@S/@r/@s markers store the bit
            # directly, bypassing the Vak→CIP parity transform. This makes
            # the canonical form idempotent — the stored bit is frame-
            # independent, so re-parsing with a different DFS order produces
            # the same bit. The canonicalizer emits @R/@S (not @/@@), so
            # the bit never goes through a DFS-relative transform again.
            if getattr(atom, '_cip_absolute', False):
                bit = getattr(atom, '_cip_bit', 0)
                atom.chirality = CHI_CW if bit == 1 else CHI_CCW
                if not hasattr(self.mol, 'chiral_centers'):
                    self.mol.chiral_centers = {}
                if not hasattr(self.mol, '_chiral_ref_nbrs'):
                    self.mol._chiral_ref_nbrs = {}
                self.mol.chiral_centers[idx] = bit
                self.mol._chiral_ref_nbrs[idx] = list(getattr(atom, '_initial_nbrs', []))
                self.mol._cip_based_stereo = True
                continue

            # Atropisomer / allene axial centre: parity is encoded directly
            # in the marker (@AX1 = CW, @AX2 = CCW). No CIP resolution needed
            # because axial stereo doesn't follow the tetrahedral 4-neighbour
            # model. Store the bit and skip the CIP parity computation.
            if getattr(atom, '_is_allene_centre', False):
                parity = getattr(atom, '_allene_parity', 0)
                if parity == 0:
                    atom.chirality = CHI_NONE
                    continue
                # _allene_parity: 1 = CW, 2 = CCW
                # SCRIPT bit convention: 0 = CCW (@), 1 = CW (@@)
                bit = 1 if parity == 1 else 0
                atom.chirality = CHI_CW if bit == 1 else CHI_CCW
                if not hasattr(self.mol, 'chiral_centers'):
                    self.mol.chiral_centers = {}
                if not hasattr(self.mol, '_chiral_ref_nbrs'):
                    self.mol._chiral_ref_nbrs = {}
                self.mol.chiral_centers[idx] = bit
                self.mol._chiral_ref_nbrs[idx] = list(getattr(atom, '_initial_nbrs', []))
                self.mol._cip_based_stereo = True
                continue

            # Vak Order: neighbors in the order they were spoken during DFS
            vak_order = list(getattr(atom, '_initial_nbrs', []))

            # Lopa Rule: if atom has implicit/explicit H not already in vak_order,
            # insert H as a ghost at position 0 (lowest priority, spoken first)
            hcount = getattr(atom, 'implicit_hs', None)
            has_explicit_h = any(
                n >= 0 and n < len(self.mol.atoms) and self.mol.atoms[n].atomic_num == 1
                for n in vak_order
            )
            if not has_explicit_h and hcount and hcount > 0:
                # H ghost: atomic_num=1
                if -1 not in vak_order:
                    vak_order = [-1] + vak_order # -1 = implicit H sentinel

            # Lopa Rule (extended): for 3-coordinate tetrahedral centres on
            # elements with a stereochemically-active lone pair (S, N, P, Se,
            # Te, As, Sb), the lone pair is the 4th ghost neighbour. Without
            # this, sulfoxides and similar 3-coordinate chiral centres would
            # never reach the 4-neighbour threshold and would be silently
            # dropped as "not enough neighbours to define chirality".
            LONE_PAIR_ELEMENTS = {'S', 'N', 'P', 'Se', 'Te', 'As', 'Sb'}
            if (atom.symbol in LONE_PAIR_ELEMENTS
                    and len(vak_order) == 3
                    and hcount is not None and hcount == 0):
                vak_order = [-2] + vak_order  # -2 = lone pair sentinel

            if len(vak_order) < 3:
                # Not enough neighbors to define chirality.  Note: with a
                # lone-pair ghost (-2), 2 real neighbours + ghost = 3 total
                # is still insufficient for tetrahedral stereo (need 4).
                atom.chirality = CHI_NONE
                continue

            if len(vak_order) < 4:
                # 3 neighbours (with or without a ghost H) is not enough for
                # a tetrahedral centre. Need exactly 4.
                atom.chirality = CHI_NONE
                continue

            # Sandhi: sort the 4 neighbors by recursive CIP priority (descending)
            from .cip import compute_cip_priorities
            # Pass the pre-computed rank_map to avoid re-running Morgan-WL
            # once per chiral atom.
            cip_sorted = compute_cip_priorities(self.mol, idx, rank_map=rank_map)

            # Permutation sign from Vak Order to CIP-Sorted Order (Rank-Stable)
            from .cip import permutation_parity
            perm_parity = permutation_parity(vak_order, cip_sorted, ranks=rank_map)

            # Pratyaya interpretation:
            # tag=2 (@)  = CCW in Vak Order
            # tag=1 (@@) = CW  in Vak Order
            if tag == 2:  # @
                base_parity = CHI_CCW
            else:          # @@
                base_parity = CHI_CW

            # If the permutation from Vak -> CIP is odd, the sense is inverted
            if perm_parity == 1:
                base_parity = CHI_CW if base_parity == CHI_CCW else CHI_CCW

            atom.chirality = base_parity
            
            # Unify with Standard Registers (for SCRIPTCanonicalizer)
            if not hasattr(self.mol, 'chiral_centers'):
                self.mol.chiral_centers = {}
            if not hasattr(self.mol, '_chiral_ref_nbrs'):
                self.mol._chiral_ref_nbrs = {}
            
            # SCRIPT bit convention: 0=@ (CCW), 1=@@ (CW)
            stored_bit = 0 if base_parity == CHI_CCW else 1
            self.mol.chiral_centers[idx] = stored_bit
            self.mol._chiral_ref_nbrs[idx] = cip_sorted
            # Mark that this mol has CIP-based stereochemistry
            self.mol._cip_based_stereo = True

            # Sthiti transform: mark the atom as CIP-absolute so the
            # canonicalizer emits @R/@S (frame-independent) instead of
            # @/@@ (DFS-relative). This makes canonical forms idempotent
            # for ALL chiral molecules, including ring-containing ones
            # where the DFS order changes between input and canonical form.
            atom._cip_absolute = True
            atom._cip_bit = stored_bit

        # Yatha-samkhya post-pass: detect pseudoasymmetric centers.
        # A tetrahedral center is pseudoasymmetric (r/s) if it has exactly
        # one pair of enantiomorphic substituents — same topology but
        # opposite embedded R/S configuration. CIP assigns lowercase r/s
        # for these centers, distinct from uppercase R/S.
        self._detect_pseudoasymmetric_centers()

        # Sthana post-pass: detect planar chirality in metallocenes.
        # A metal (Fe, Ru, Os, etc.) bonded to 2+ aromatic 5-rings (Cp)
        # with substituents that break mirror symmetry is planar chiral.
        # Note: 2D graph cannot distinguish Rp from Sp — both enantiomers
        # have identical connectivity. We mark the structure as PLANAR
        # with variant=0 (unspecified); the user must use @PL1/@PL2 to
        # distinguish enantiomers.
        self._detect_planar_chirality()

    def _detect_planar_chirality(self):
        """Sthana: detect planar chirality in metallocene structures.

        A metallocene (e.g., ferrocene) consists of a metal sandwiched
        between two aromatic 5-membered rings (cyclopentadienyl). When
        substituents on a ring break mirror symmetry, the complex is
        planar chiral.

        Detection:
        1. Find metal atoms (transition metals)
        2. For each metal, find aromatic 5-rings bonded to it
        3. If 2+ such rings exist, check for asymmetric substitution
        4. Mark the metal as stereo_type=PLANAR with variant=0

        Limitation: 2D graph cannot distinguish Rp from Sp (enantiomers
        have identical connectivity). The user must specify @PL1/@PL2
        explicitly to distinguish enantiomers.
        """
        METAL_ATOMIC_NUMS = {
            26,  # Fe
            27,  # Co
            28,  # Ni
            44,  # Ru
            45,  # Rh
            46,  # Pd
            76,  # Os
            77,  # Ir
            78,  # Pt
            24,  # Cr
            25,  # Mn
            74,  # W
            42,  # Mo
            23,  # V
            22,  # Ti
            40,  # Zr
            72,  # Hf
        }

        for idx, atom in enumerate(self.mol.atoms):
            if atom.atomic_num not in METAL_ATOMIC_NUMS:
                continue

            # Find aromatic 5-rings bonded to this metal
            cp_rings = self._find_cp_rings(idx)
            if len(cp_rings) < 2:
                continue

            # Check if any ring has asymmetric substitution
            # (2+ different substituents at non-equivalent positions)
            for ring in cp_rings:
                if self._ring_has_asymmetric_substitution(ring, exclude_idx=idx):
                    # Mark the metal as planar chiral (variant=0 = unspecified)
                    atom.stereo_type = StereoType.PLANAR
                    atom._polyhedral_variant = 0
                    if not hasattr(self.mol, '_planar_chiral_metals'):
                        self.mol._planar_chiral_metals = []
                    self.mol._planar_chiral_metals.append(idx)
                    break  # one asymmetric ring is enough

    def _find_cp_rings(self, metal_idx: int) -> list:
        """Find aromatic 5-membered rings bonded to the metal.

        A Cp ring is a 5-membered ring where all 5 atoms are bonded to
        the metal (eta-5 coordination) OR at least 3 are bonded (eta-3).
        For simplicity, we check for 5-membered rings where the metal
        is bonded to at least one ring atom and the ring is aromatic.
        """
        rings = []
        # Get all atoms bonded to the metal
        metal_nbrs = set(n for n, _ in self.mol.adj.get(metal_idx, []))
        if not metal_nbrs:
            return rings

        # Find all 5-membered rings in the molecule
        all_rings = self._find_all_rings(max_size=5)
        for ring in all_rings:
            if len(ring) != 5:
                continue
            # Check if ring is aromatic (all bonds aromatic)
            is_aromatic = True
            for i in range(len(ring)):
                a, b = ring[i], ring[(i + 1) % len(ring)]
                bond = self.mol.get_bond(a, b)
                if bond is None or bond.bond_type != 4:  # 4 = AROMATIC
                    is_aromatic = False
                    break
            if not is_aromatic:
                continue
            # Check if metal is bonded to at least one ring atom
            if metal_nbrs & set(ring):
                rings.append(ring)
        return rings

    def _find_all_rings(self, max_size: int = 6) -> list:
        """Find all rings up to max_size using DFS."""
        rings = []
        seen_rings = set()

        def dfs(start: int, current: int, path: list, visited: set):
            if len(path) > max_size:
                return
            for nbr, _ in self.mol.adj.get(current, []):
                if nbr == start and len(path) >= 3:
                    # Found a ring
                    ring = tuple(sorted(path))
                    if ring not in seen_rings:
                        seen_rings.add(ring)
                        rings.append(list(path))
                elif nbr not in visited and nbr > start:
                    # Only explore atoms with higher index to avoid duplicates
                    dfs(start, nbr, path + [nbr], visited | {nbr})

        for start in range(len(self.mol.atoms)):
            dfs(start, start, [start], {start})
        return rings

    def _ring_has_asymmetric_substitution(self, ring: list, exclude_idx: int) -> bool:
        """Check if a ring has 2+ different substituents at non-equivalent
        positions (i.e., the substitution pattern breaks mirror symmetry).

        For a Cp ring, if there are 2+ substituents (excluding H and the
        metal bond) at different ring positions, the ring is prochiral.
        """
        substituents = []
        for ring_atom in ring:
            for nbr, _ in self.mol.adj.get(ring_atom, []):
                if nbr == exclude_idx:
                    continue
                if nbr in ring:
                    continue
                # This is an external substituent
                nbr_atom = self.mol.atoms[nbr]
                # Skip implicit H (atomic_num=1 with no further substituents)
                if nbr_atom.atomic_num == 1:
                    continue
                substituents.append((ring_atom, nbr_atom.atomic_num))

        # Need 2+ substituents at DIFFERENT ring positions
        if len(substituents) < 2:
            return False

        # Check if substituents are at different positions with different types
        positions = set(s[0] for s in substituents)
        if len(positions) < 2:
            return False  # all substituents at same position (not a ring)

        # Check if at least 2 substituents have different atomic numbers
        atomic_nums = set(s[1] for s in substituents)
        if len(atomic_nums) < 2:
            # Same substituent type — could still be chiral if at non-symmetric
            # positions, but for simplicity we require different types
            return False

        return True

    def _detect_pseudoasymmetric_centers(self):
        """Yatha-samkhya: detect pseudoasymmetric (r/s) centers.

        For each resolved tetrahedral chiral center, examine its 4
        substituents pairwise. If exactly one pair is enantiomorphic
        (same topology, opposite embedded R/S), mark the center as
        pseudoasymmetric so the canonicalizer emits @r/@s instead of
        @R/@S.

        Two substituents are enantiomorphic if:
        - Their CIP priority tuples are EQUAL (same topology)
        - BUT the stereocenters embedded in them have OPPOSITE configs
          (one is R where the other is S, at corresponding positions)
        """
        if not hasattr(self.mol, 'chiral_centers'):
            return

        for idx, atom in enumerate(self.mol.atoms):
            # Only check resolved tetrahedral centers (not already marked
            # pseudoasymmetric by explicit @r/@s input).
            if getattr(atom, '_pseudoasymmetric', False):
                continue
            if idx not in self.mol.chiral_centers:
                continue
            if getattr(atom, 'stereo_type', StereoType.NONE) != StereoType.TETRAHEDRAL:
                continue

            # Get the 4 substituents (neighbors), handling ghost neighbors
            nbrs = list(getattr(atom, '_initial_nbrs', []))
            hcount = getattr(atom, 'implicit_hs', 0) or 0
            if hcount > 0:
                nbrs = nbrs + [-1]  # implicit H

            # Lone-pair ghost for 3-coordinate S/N/P
            LONE_PAIR_ELEMENTS = {'S', 'N', 'P', 'Se', 'Te', 'As', 'Sb'}
            if (atom.symbol in LONE_PAIR_ELEMENTS
                    and len(nbrs) == 3
                    and hcount == 0):
                nbrs = nbrs + [-2]

            if len(nbrs) != 4:
                continue

            # Check all pairs for enantiomorphism
            enantiomorphic_pairs = 0
            for i in range(len(nbrs)):
                for j in range(i + 1, len(nbrs)):
                    if self._are_enantiomorphic(idx, nbrs[i], nbrs[j]):
                        enantiomorphic_pairs += 1

            # Pseudoasymmetric iff exactly one enantiomorphic pair
            # (the other two substituents must be distinct, which is the
            # common case; if there were two pairs it'd be non-stereogenic)
            if enantiomorphic_pairs == 1:
                atom._pseudoasymmetric = True

    def _are_enantiomorphic(self, center_idx: int, nbr_a: int, nbr_b: int,
                             visited: frozenset = None) -> bool:
        """Check if two substituents (rooted at nbr_a and nbr_b, excluding
        center_idx) are enantiomorphic — same topology but mirror images.

        Two substituents are enantiomorphic iff:
        1. Their atomic compositions and connectivities are identical
           (same topology)
        2. At every corresponding stereocenter in the two trees, the
           configurations are OPPOSITE (one R, one S)
        3. There is at least one stereocenter (otherwise they're just
           identical/homomorphic, not enantiomorphic)
        """
        if visited is None:
            visited = frozenset({center_idx})

        # Ghost neighbors can't be enantiomorphic with anything
        if nbr_a < 0 or nbr_b < 0:
            return False

        if nbr_a >= len(self.mol.atoms) or nbr_b >= len(self.mol.atoms):
            return False

        atom_a = self.mol.atoms[nbr_a]
        atom_b = self.mol.atoms[nbr_b]

        # Must be same element
        if atom_a.atomic_num != atom_b.atomic_num:
            return False

        # Must have same formal charge
        if atom_a.formal_charge != atom_b.formal_charge:
            return False

        # Must have same isotope
        if atom_a.isotope != atom_b.isotope:
            return False

        # Must have same number of neighbors (excluding the path back)
        nbrs_a = [n for n, _ in self.mol.adj.get(nbr_a, []) if n not in visited]
        nbrs_b = [n for n, _ in self.mol.adj.get(nbr_b, []) if n not in visited]
        if len(nbrs_a) != len(nbrs_b):
            return False

        # Check stereo: if both are chiral centers, they must be OPPOSITE
        # configs to be enantiomorphic. If only one is chiral, they're not
        # enantiomorphic (different topology).
        a_chiral = nbr_a in self.mol.chiral_centers
        b_chiral = nbr_b in self.mol.chiral_centers
        if a_chiral != b_chiral:
            return False
        if a_chiral and b_chiral:
            if self.mol.chiral_centers[nbr_a] == self.mol.chiral_centers[nbr_b]:
                # Same config → identical, not enantiomorphic
                return False
            # Opposite config at this center — continue checking subtrees
        # If neither is chiral, continue (no stereo to compare yet)

        # Recursively check that subtrees match (enantiomorphically)
        # Sort neighbors by atomic number for canonical matching
        def sort_key(n):
            if n < 0 or n >= len(self.mol.atoms):
                return (0,)
            return (self.mol.atoms[n].atomic_num,)

        nbrs_a_sorted = sorted(nbrs_a, key=sort_key)
        nbrs_b_sorted = sorted(nbrs_b, key=sort_key)

        # For simple cases (no further stereo), check that paired subtrees
        # are identical. For stereo-bearing subtrees, they must be enantiomorphic.
        # We use a conservative check: pair by sorted atomic number and
        # require each pair to be at least identical (homomorphic). The
        # enantiomorphic check at the root level is what makes the whole
        # substituent enantiomorphic.
        for na, nb in zip(nbrs_a_sorted, nbrs_b_sorted):
            if not self._subtrees_identical(na, nb, visited | {nbr_a, nbr_b}):
                return False

        # At least one of the substituent roots must be a chiral center
        # with opposite config (otherwise they're just identical)
        return a_chiral and b_chiral

    def _subtrees_identical(self, root_a: int, root_b: int,
                             visited: frozenset) -> bool:
        """Check if two subtrees have identical topology (ignoring stereo)."""
        if root_a < 0 or root_b < 0:
            return root_a == root_b
        if root_a >= len(self.mol.atoms) or root_b >= len(self.mol.atoms):
            return False

        atom_a = self.mol.atoms[root_a]
        atom_b = self.mol.atoms[root_b]
        if atom_a.atomic_num != atom_b.atomic_num:
            return False
        if atom_a.formal_charge != atom_b.formal_charge:
            return False

        nbrs_a = [n for n, _ in self.mol.adj.get(root_a, []) if n not in visited]
        nbrs_b = [n for n, _ in self.mol.adj.get(root_b, []) if n not in visited]
        if len(nbrs_a) != len(nbrs_b):
            return False

        # Match by sorted atomic number
        def sort_key(n):
            if n < 0 or n >= len(self.mol.atoms):
                return (0,)
            return (self.mol.atoms[n].atomic_num,)
        nbrs_a_sorted = sorted(nbrs_a, key=sort_key)
        nbrs_b_sorted = sorted(nbrs_b, key=sort_key)

        new_visited = visited | {root_a, root_b}
        for na, nb in zip(nbrs_a_sorted, nbrs_b_sorted):
            if not self._subtrees_identical(na, nb, new_visited):
                return False
        return True

    def _get_cip_priority(self, n_idx: int, parent_idx: int, depth: int, visited: set = None) -> tuple:
        """
        Recursive CIP priority assignment.
        Reference: (atomic_num, mass, -charge, (sorted_neighbor_priorities...))
        """
        if n_idx == -1: # Implicit H
            return (1, 1, 0, ())
            
        if depth == 0:
            return (0, 0, 0, ())

        if visited is None:
            visited = {parent_idx}
        
        if n_idx in visited:
            return (0, 0, 0, ()) # Avoid cycles

        atom = self.mol.atoms[n_idx]
        new_visited = visited | {n_idx}
        
        # Collect priorities of neighbors
        nbr_priorities = []
        for nbr_idx, _ in self.mol.adj.get(n_idx, []):
            if nbr_idx == parent_idx:
                continue
            nbr_priorities.append(self._get_cip_priority(nbr_idx, n_idx, depth - 1, new_visited))
            
        # Add implicit H priorities
        hcount = getattr(atom, 'implicit_hs', 0) or 0
        for _ in range(hcount):
            nbr_priorities.append((1, 1, 0, ()))
            
        # Sort neighbor priorities descending for comparison
        nbr_priorities.sort(reverse=True)
        
        mass = atom.isotope if atom.isotope > 0 else _standard_mass(atom.atomic_num)
        return (atom.atomic_num, mass, -atom.formal_charge, tuple(nbr_priorities))

    def _permutation_sign(self, original: List, target: List) -> int:
        """
        Compute the sign of the permutation that transforms `original` into `target`.
        Returns +1 (even) or -1 (odd).
        """
        if len(original) != len(target):
            return 1
            
        # Build mapping from value to position in target
        pos = {v: i for i, v in enumerate(target)}
        
        # Permutation list
        try:
            perm = [pos[v] for v in original]
        except KeyError:
            return 1
            
        # Count inversions
        inversions = 0
        for i in range(len(perm)):
            for j in range(i + 1, len(perm)):
                if perm[i] > perm[j]:
                    inversions += 1
        
        return 1 if inversions % 2 == 0 else -1
