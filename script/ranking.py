"""
Standalone SCRIPT Ranking Engine
Implements a deterministic Morgan (Weisfeiler-Lehman) algorithm for atom invariants.
"""

from typing import Dict, List, Tuple

def calculate_ranks(mol) -> Dict[int, int]:
    """
    Assign a canonical rank to each atom in the molecule based on its environment.
    Ranks are deterministic and independent of the input atom order.
    Supports script.mol.CoreMolecule.

    For periodic molecules (mol.is_periodic = True), the lattice translation
    vector (tx, ty, tz) on each bond is included in the neighbor invariant so
    that atoms connected via different lattice directions receive distinct ranks.
    This satisfies the LQG canonicalization requirement from Krenn et al. (2022).

    Caching: the rank map is memoised on the CoreMolecule instance keyed by
    ``mol._graph_version``, which is bumped on every add_atom/add_bond call.
    This makes repeated calls (e.g. once per chiral atom in
    ChiralResolver.resolve) effectively free. The cache is opt-in: callers
    that pass a non-CoreMolecule object fall back to the uncached path.
    """
    # Opt-in cache: only CoreMolecule exposes _graph_version / _rank_cache.
    cache_version = getattr(mol, "_graph_version", None)
    if cache_version is not None and \
       getattr(mol, "_rank_cache", None) is not None and \
       getattr(mol, "_rank_cache_version", None) == cache_version:
        return mol._rank_cache

    rank_map = _calculate_ranks_uncached(mol)

    # Store back on the instance if the cache slot is available.
    if cache_version is not None:
        try:
            mol._rank_cache = rank_map
            mol._rank_cache_version = cache_version
        except (AttributeError, TypeError):
            pass

    return rank_map


def _calculate_ranks_uncached(mol) -> Dict[int, int]:
    """Original Morgan/WL implementation. See calculate_ranks for the cache wrapper."""
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
    """Map bond types to stable integers.

    Fast path: BondType enum from script.mol exposes ``.value`` (an int).
    Slow path: fall back to substring matching for RDKit BondType or
    other duck-typed bond objects that don't expose a value.
    """
    v = getattr(bt, "value", None)
    if isinstance(v, int):
        # script.mol.BondType values: SINGLE=1, DOUBLE=2, TRIPLE=3, AROMATIC=4,
        # DATIVE=5, REV_DATIVE=6, TAUTOMERIC=7, COORDINATE=8, STAR=9.
        # All non-aromatic typed bonds collapse to "other" (0) for ranking
        # purposes, matching the original substring-based semantics.
        if v in (1, 2, 3, 4):
            return v
        return 0
    # Fallback for RDKit BondType (no .value) or bare ints.
    if isinstance(bt, int):
        if bt in (1, 2, 3, 4):
            return bt
        return 0
    s = str(bt)
    if 'SINGLE' in s: return 1
    if 'DOUBLE' in s: return 2
    if 'TRIPLE' in s: return 3
    if 'AROMATIC' in s: return 4
    return 0

def _stable_hash(obj) -> int:
    """Produces a stable 64-bit integer hash from a python object.

    Original implementation used SHA-256, which dominated the parser's
    runtime for chiral molecules (1152 SHA-256 calls per glucose parse).
    The rank map only needs to be stable *within a single process* — it
    is never serialised or compared across runs — so Python's built-in
    ``hash()`` is sufficient and ~50× faster. We mask to 63 bits to keep
    the value positive (avoids edge cases with negative hashes in tuple
    comparisons downstream).
    """
    return hash(obj) & 0x7FFFFFFFFFFFFFFF

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
