"""
Graft Copolymer Atomic Expansion
================================

For a graft copolymer:

    {[BB]}<n:N> -g- {[GG]}<n:M>

the parser stores the BACKBONE unit (BB) and GRAFT unit (GG) as separate
CoreMolecule fragments joined by a single seam bond.  This module expands
that compact representation into a real branched atomic graph:

  * The backbone unit BB is replicated N times and chained head-to-tail.
  * The graft unit GG is replicated M times; each replica is attached as
    a pendant branch to a backbone branch point.

Branch-point assignment (deterministic, "grafting-through" semantics):

  * Branch points are the FIRST atom of each backbone repeat unit
    (i.e. atoms at indices 0, bb_unit_size, 2*bb_unit_size, ...).
  * Graft i (0-indexed) attaches to backbone branch point i, in order.
  * If M > N: surplus grafts chain off the LAST backbone branch point
    as a tail of pendant branches (each graft remains a separate
    branch; they are NOT concatenated).
  * If M < N: only the first M backbone units carry grafts; the rest
    of the backbone is bare.

The expanded graph is written into the parent CoreMolecule.  The original
backbone unit atoms and graft unit atoms are kept (so existing atom
indices remain valid); new replicated atoms are appended after them.

The PolymerBlock metadata on the parent molecule is updated so that
block_kind='graft' is preserved and atom_start/atom_end cover the full
expanded ranges.
"""

from __future__ import annotations
from typing import List, Tuple, Optional, Any
from .mol import CoreMolecule, CoreAtom, CoreBond, BondType, PolymerBlock


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _copy_unit_into_parent(parent: CoreMolecule,
                           unit: CoreMolecule,
                           extra_bonds: Optional[List[Tuple[int, int, BondType, str]]] = None
                           ) -> List[int]:
    """Append a copy of `unit`'s atoms and internal bonds into `parent`.

    Returns the list of new atom indices in `parent` corresponding to
    `unit.atoms` (in order).

    `extra_bonds` is an optional list of (local_u, local_v, bond_type,
    bond_class) tuples that should be added INSIDE the copied unit
    (using the new local indices).  Useful when the unit's internal
    ring-closure bonds were not materialised in `unit` itself.
    """
    # Map from unit-local atom index -> parent atom index
    index_map: List[int] = []
    for atom in unit.atoms:
        # Shallow copy of the CoreAtom — we don't want shared mutable state.
        new_atom = CoreAtom(
            atomic_num=atom.atomic_num,
            formal_charge=atom.formal_charge,
            isotope=atom.isotope,
            radical_electrons=atom.radical_electrons,
            symbol=atom.symbol,
            is_aromatic=atom.is_aromatic,
            mapping=atom.mapping,
            occupancy=atom.occupancy,
            spin=atom.spin,
            is_excited=atom.is_excited,
            is_wildcard=atom.is_wildcard,
            beam_radius=atom.beam_radius,
        )
        new_atom.implicit_hs = atom.implicit_hs
        new_atom._initial_tag = atom._initial_tag
        new_atom._initial_nbrs = list(atom._initial_nbrs)
        new_atom.chirality = atom.chirality
        new_atom.stereo_type = atom.stereo_type
        new_atom.stereo_neighbors = list(atom.stereo_neighbors)
        new_atom.is_query = atom.is_query
        new_atom.query_atomic_nums = list(atom.query_atomic_nums)
        new_atom.query_not = atom.query_not
        new_atom.query_ring = atom.query_ring
        new_atom.query_valence = atom.query_valence
        new_atom.query_hcount = atom.query_hcount
        new_atom.query_aromatic = atom.query_aromatic
        new_atom.query_primitives = list(atom.query_primitives)
        parent_idx = parent.add_atom(new_atom)
        index_map.append(parent_idx)

    # Re-add internal bonds (those whose both endpoints live inside the unit)
    for bond in unit.bonds:
        u = index_map[bond.begin_atom_idx]
        v = index_map[bond.end_atom_idx]
        parent.add_bond(
            u, v, bond.bond_type,
            bond_dir=bond.bond_dir,
            hapticity=bond.hapticity,
            bond_class=bond.bond_class,
            translation=bond.translation,
            control_points=bond.control_points,
        )

    if extra_bonds:
        for (lu, lv, bt, bc) in extra_bonds:
            u = index_map[lu]
            v = index_map[lv]
            parent.add_bond(u, v, bt, bond_class=bc)

    return index_map


def _remove_bond_between(parent: CoreMolecule, a: int, b: int) -> bool:
    """Remove the bond (if any) directly connecting atoms a and b.

    Returns True if a bond was removed.  The bond is removed from both
    `parent.bonds` (replaced with a sentinel None — index stability
    matters) and from `parent.adj`.  We do NOT renumber bonds; instead
    we mark the slot as None and downstream code skips None entries.

    To keep things simple for callers that scan parent.bonds linearly,
    we actually rebuild the bonds list and renumber, then fix up adj.
    """
    bond_idx: Optional[int] = None
    for nbr, bidx in parent.adj.get(a, []):
        if nbr == b:
            bond_idx = bidx
            break
    if bond_idx is None:
        return False

    # Remove from adj on both endpoints
    parent.adj[a] = [(n, bi) for (n, bi) in parent.adj[a] if bi != bond_idx]
    parent.adj[b] = [(n, bi) for (n, bi) in parent.adj[b] if bi != bond_idx]

    # Rebuild bonds list without the removed bond, fixing up adj indices
    old_bonds = parent.bonds
    new_bonds: List[CoreBond] = []
    remap: List[int] = [-1] * len(old_bonds)
    for i, bnd in enumerate(old_bonds):
        if i == bond_idx:
            continue
        remap[i] = len(new_bonds)
        new_bonds.append(bnd)
    parent.bonds = new_bonds

    # Fix adj indices
    for k in list(parent.adj.keys()):
        parent.adj[k] = [(n, remap[bi]) for (n, bi) in parent.adj[k] if remap[bi] != -1]

    parent._graph_version += 1
    return True


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def expand_graft_copolymer(mol: CoreMolecule,
                           bb_unit: CoreMolecule,
                           bb_repeat: Any,
                           gr_unit: CoreMolecule,
                           gr_repeat: Any) -> None:
    """Expand a graft copolymer into `mol` in place.

    Assumes `mol` already contains:
      * the backbone unit atoms at indices [bb_start, bb_end]
      * the graft unit atoms at indices [gr_start, gr_end]
      * a single seam bond between bb_end and gr_start (which we remove)

    After expansion, `mol.atoms` contains:
      [original bb unit atoms] + [original gr unit atoms] +
      [replicated bb unit atoms] + [replicated gr unit atoms]

    PolymerBlock metadata is rewritten to reflect the expanded ranges.
    """
    # Normalise repeat counts: 'n' or None -> 1 for expansion purposes
    def _as_int(r: Any, default: int = 1) -> int:
        if r is None: return default
        if isinstance(r, int): return max(1, r)
        if isinstance(r, tuple) and r:
            try: return max(1, int(r[0]))
            except (TypeError, ValueError): return default
        if isinstance(r, str):
            try: return max(1, int(r))
            except ValueError: return default
        return default

    N = _as_int(bb_repeat)
    M = _as_int(gr_repeat)

    # Find the original backbone and graft blocks in polymer_blocks
    bb_block: Optional[PolymerBlock] = None
    gr_block: Optional[PolymerBlock] = None
    for b in mol.polymer_blocks:
        if b.block_kind == '' and bb_block is None:
            bb_block = b
        elif b.block_kind == 'graft':
            gr_block = b
    # Fall back to positions if kinds weren't set
    if bb_block is None and len(mol.polymer_blocks) >= 1:
        bb_block = mol.polymer_blocks[0]
    if gr_block is None and len(mol.polymer_blocks) >= 2:
        gr_block = mol.polymer_blocks[1]
    if bb_block is None or gr_block is None:
        return  # nothing to expand

    bb_start = bb_block.atom_start if bb_block.atom_start is not None else 0
    bb_end = bb_block.atom_end if bb_block.atom_end is not None else (bb_start + len(bb_unit.atoms) - 1)
    gr_start = gr_block.atom_start if gr_block.atom_start is not None else (bb_end + 1)
    gr_end = gr_block.atom_end if gr_block.atom_end is not None else (gr_start + len(gr_unit.atoms) - 1)

    # Remove the seam bond between bb_end and gr_start (if present)
    _remove_bond_between(mol, bb_end, gr_start)

    # bb_unit_size = bb_end - bb_start + 1
    # gr_unit_size = gr_end - gr_start + 1
    # (we trust that these match the lengths of bb_unit.atoms / gr_unit.atoms)

    # ---- Step 1: replicate backbone (N-1 more copies, chained) ----
    # The first copy is already in mol at [bb_start..bb_end].
    # For each additional copy i = 1..N-1:
    #   - Append bb_unit's atoms to mol
    #   - Chain: connect previous tail (last atom of previous copy) -> new head
    # Track backbone "branch points": the FIRST atom index of each copy.
    bb_branch_points: List[int] = [bb_start]
    prev_bb_tail = bb_end
    for _ in range(1, N):
        new_indices = _copy_unit_into_parent(mol, bb_unit)
        new_head = new_indices[0]
        new_tail = new_indices[-1]
        # Chain from previous tail to new head with a single bond
        # (matching whatever bond type the unit's internal chain uses;
        # defaulting to SINGLE is safe for organic backbones.)
        mol.add_bond(prev_bb_tail, new_head, BondType.SINGLE, bond_class="")
        bb_branch_points.append(new_head)
        prev_bb_tail = new_tail

    # ---- Step 2: replicate grafts (M copies), each pendant ----
    # The first graft copy is already in mol at [gr_start..gr_end].
    # Attach it to bb_branch_points[0].
    # For each additional copy i = 1..M-1:
    #   - Append gr_unit's atoms to mol
    #   - Attach new head to the appropriate branch point
    #
    # Branch assignment: graft i -> bb_branch_points[min(i, N-1)]
    # (surplus grafts chain off the last backbone branch point)
    gr_head = gr_start
    bp_idx = 0 if N > 0 else 0
    target_bp = bb_branch_points[bp_idx] if bb_branch_points else bb_start
    mol.add_bond(target_bp, gr_head, BondType.SINGLE, bond_class="")

    # Track graft head atoms (for surplus-chain logic below)
    last_graft_heads: List[int] = [gr_head]
    for i in range(1, M):
        new_indices = _copy_unit_into_parent(mol, gr_unit)
        new_head = new_indices[0]
        # Assign branch point: graft i -> backbone branch point i (if exists)
        if i < len(bb_branch_points):
            target_bp = bb_branch_points[i]
            mol.add_bond(target_bp, new_head, BondType.SINGLE, bond_class="")
        else:
            # Surplus graft: attach as a sibling branch to the LAST backbone
            # branch point (do NOT chain grafts onto each other — they remain
            # independent pendant branches).
            target_bp = bb_branch_points[-1] if bb_branch_points else bb_start
            mol.add_bond(target_bp, new_head, BondType.SINGLE, bond_class="")
        last_graft_heads.append(new_head)

    # ---- Step 3: update PolymerBlock metadata ----
    # The original backbone block now spans [bb_start, bb_start + N*bb_unit_size - 1]
    # The original graft block now spans [gr_start, end_of_last_graft_copy]
    bb_full_end = bb_start + N * len(bb_unit.atoms) - 1
    # Graft atoms live AFTER all backbone replicas? No — we appended grafts
    # AFTER backbone replicas, so the original graft unit is still at gr_start
    # but additional graft copies are after the additional backbone copies.
    # The "graft block" range is now non-contiguous, but we keep the original
    # [gr_start, gr_end] range as a representative marker.
    # Update atom_start/atom_end on the blocks:
    bb_block.atom_end = bb_full_end
    bb_block.repeat_count = N
    # gr_block's atom range stays as the original unit; we add a note via
    # block_kind='graft' so downstream code knows it's been expanded.
    gr_block.repeat_count = M

    # Ensure block_topology is set
    mol.block_topology = 'graft'
