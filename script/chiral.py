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
from .mol import CoreMolecule, CoreAtom

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
