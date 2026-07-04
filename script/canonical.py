from typing import List, Dict, Optional, Tuple
from .mol import CoreMolecule, BondType, StereoType
from .ranking import calculate_ranks
from .stereo import get_chiral_symbol, perceive_chirality
import math


def _lattice_vectors_to_params(
    lv: Tuple[Tuple[float,float,float], ...]
) -> Optional[Tuple[float,float,float,float,float,float]]:
    """
    Convert a 3x3 row-vector lattice matrix back to
    conventional cell parameters (a, b, c, alpha, beta, gamma).
    Returns None if the matrix is degenerate.
    """
    try:
        a_vec = lv[0]; b_vec = lv[1]; c_vec = lv[2]
        a = math.sqrt(sum(x**2 for x in a_vec))
        b = math.sqrt(sum(x**2 for x in b_vec))
        c = math.sqrt(sum(x**2 for x in c_vec))
        if a < 1e-9 or b < 1e-9 or c < 1e-9:
            return None
        cos_alpha = sum(b_vec[i]*c_vec[i] for i in range(3)) / (b*c)
        cos_beta  = sum(a_vec[i]*c_vec[i] for i in range(3)) / (a*c)
        cos_gamma = sum(a_vec[i]*b_vec[i] for i in range(3)) / (a*b)
        alpha = round(math.degrees(math.acos(max(-1.0, min(1.0, cos_alpha)))), 4)
        beta  = round(math.degrees(math.acos(max(-1.0, min(1.0, cos_beta)))),  4)
        gamma = round(math.degrees(math.acos(max(-1.0, min(1.0, cos_gamma)))), 4)
        return round(a, 4), round(b, 4), round(c, 4), alpha, beta, gamma
    except Exception:
        return None

# Organic subset for bare symbols (no brackets)
ORGANIC_SUBSET = {'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}

# Default valences for bare symbols
DEFAULT_VALENCE = {
    'B': 3, 'C': 4, 'N': 3, 'O': 2, 'P': 3, 'S': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1
}

class SCRIPTCanonicalizer:
    """
    Produces canonical SCRIPT strings using RDKit-independent logic.
    """

    def canonicalize_mol(self, mol):
        """Deprecated: Use canonicalize_core or rdkit_bridge.SCRIPTFromMol."""
        return self.canonicalize_core(mol)

    def canonicalize_core(self, mol: CoreMolecule) -> Optional[str]:
        """Convert a CoreMolecule to its canonical SCRIPT string."""
        if not mol.atoms:
            return None

        # Find all connected components (handles salts, multi-fragment molecules)
        components = self._find_components(mol)

        if len(components) == 1:
            body = self._canonicalize_component(mol, list(range(len(mol.atoms))))
        else:
            fragment_strings = []
            for component in components:
                s = self._canonicalize_component(mol, component)
                if s:
                    fragment_strings.append(s)
            fragment_strings = [s.rstrip('.') for s in fragment_strings]
            fragment_strings.sort(key=lambda s: (-len(s), s))
            body = '.'.join(fragment_strings)

        # V3.6: prepend crystallographic context + lattice parameters if present
        prefix = ''
        ctx = getattr(mol, 'macroscopic_context', None)
        lv  = getattr(mol, 'lattice_vectors',     None)
        if ctx:
            if lv is not None:
                # Convert 3x3 lattice vectors back to a,b,c,alpha,beta,gamma
                params = _lattice_vectors_to_params(lv)
                if params:
                    a, b, c, alpha, beta, gamma = params
                    prefix = (f"[[{ctx};{a:.4f},{b:.4f},{c:.4f},"
                              f"{alpha:.2f},{beta:.2f},{gamma:.2f}]]")
                else:
                    prefix = f"[[{ctx}]]"
            else:
                prefix = f"[[{ctx}]]"

        # Phase boundary
        pb = getattr(mol, 'phase_boundary', None)
        if pb:
            prefix = prefix + ' | ' if prefix else '| '

        return prefix + body if prefix else body

    def _find_components(self, mol: CoreMolecule) -> List[List[int]]:
        """Find all connected components via BFS.

        For periodic molecules, bonds with non-zero translation vectors are
        treated as connecting (so the whole lattice is one component) but BFS
        does not follow them a second time to avoid infinite traversal.
        """
        num_atoms = len(mol.atoms)
        visited = set()
        components = []
        for start in range(num_atoms):
            if start in visited:
                continue
            component = []
            queue = [start]
            visited.add(start)
            while queue:
                node = queue.pop(0)
                component.append(node)
                for nbr, bond_idx in mol.adj.get(node, []):
                    # For periodic bonds, only cross the boundary once.
                    # The bond is still used for connectivity, but we do
                    # not enqueue the neighbour if it was already visited.
                    bond = mol.bonds[bond_idx]
                    is_periodic_bond = (getattr(bond, 'translation', (0,0,0)) != (0,0,0))
                    if nbr not in visited:
                        if not is_periodic_bond:
                            visited.add(nbr)
                            queue.append(nbr)
                        else:
                            # Periodic bond to an unvisited image — include as
                            # same component but do not recurse into the image.
                            visited.add(nbr)
                            component.append(nbr)
            components.append(component)
        return components

    def _canonicalize_component(self, mol: CoreMolecule, atom_indices: List[int]) -> Optional[str]:
        """Canonicalize a single connected component.

        For periodic molecules (mol.is_periodic = True), the Morgan-WL ranking
        in ranking.py includes lattice translation vectors so topologically
        distinct sites receive distinct ranks, and _bond_symbol emits @tx,ty,tz
        on every cross-cell bond.  The result is a canonical LQG string that
        uniquely identifies the periodic net up to the choice of unit cell.
        """
        # 1. Ranking — periodic-aware via calculate_ranks
        rank_map = calculate_ranks(mol)
        num_atoms = len(mol.atoms)
        ranks = [rank_map.get(i, 0) for i in range(num_atoms)]

        # 2. DFS from lowest rank atom(s) within this component
        component_set = set(atom_indices)
        min_rank = min(ranks[i] for i in atom_indices)
        start_candidates = [i for i in atom_indices if ranks[i] == min_rank]

        if len(start_candidates) == 1:
            return self._canonicalize_component_from_root(mol, atom_indices, start_candidates[0], ranks)

        canonical_forms = []
        for start_atom in start_candidates:
            canon = self._canonicalize_component_from_root(mol, atom_indices, start_atom, ranks)
            if canon is not None:
                canonical_forms.append(canon)
        return min(canonical_forms) if canonical_forms else None

    def _canonicalize_component_from_root(self, mol: CoreMolecule, atom_indices: List[int], start_atom: int, ranks: List[int]) -> Optional[str]:
        visited = set()
        ring_bonds_set = set()
        self._find_ring_bonds(mol, start_atom, visited, ranks, -1, ring_bonds_set)

        atom_to_id_collect = {}
        dfs_neighbor_orders = {}
        self._collect_dfs_neighbor_orders(mol, start_atom, atom_to_id_collect, ranks, -1, ring_bonds_set, dfs_neighbor_orders)

        has_chiral_data = hasattr(mol, 'chiral_centers') and mol.chiral_centers
        if not has_chiral_data:
            perceive_chirality(mol, ranks, dfs_neighbor_orders)

        atom_to_id = {}
        depths = {start_atom: 1}
        ring_counter = [0]
        return self._dfs(mol, start_atom, atom_to_id, ranks, -1, ring_counter, ring_bonds_set, depths)
    
    def _find_ring_bonds(self, mol: CoreMolecule, atom_idx, visited, ranks, from_bond_idx, ring_bonds):
        """First pass: identify which bonds are ring closures.

        Periodic bonds (translation != (0,0,0)) are lattice edges, not ring
        closures in the quotient graph.  They are excluded from ring_bonds so
        that the DFS string builder does not emit Anubandha ring notation for
        lattice-spanning bonds.
        """
        visited.add(atom_idx)

        neighbors = []
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            if bond_idx == from_bond_idx:
                continue
            bond = mol.bonds[bond_idx]
            if getattr(bond, 'translation', (0,0,0)) != (0,0,0):
                continue  # skip periodic boundary bonds in ring detection
            neighbors.append((ranks[nbr_idx], nbr_idx, bond_idx))

        neighbors.sort(key=lambda x: (x[0], x[1]))

        for _, nbr_idx, bond_idx in neighbors:
            if nbr_idx in visited:
                ring_bonds.add(bond_idx)
            else:
                self._find_ring_bonds(mol, nbr_idx, visited, ranks, bond_idx, ring_bonds)
    
    def _collect_dfs_neighbor_orders(self, mol: CoreMolecule, atom_idx, atom_to_id, ranks, from_bond_idx, ring_bonds_set, dfs_orders):
        """Collect DFS neighbor orders for stereochemistry reference."""
        atom_string_id = len(atom_to_id)
        atom_to_id[atom_idx] = atom_string_id
        atom = mol.atoms[atom_idx]
        
        parent_idx = -1
        if from_bond_idx >= 0:
            parent_idx = mol.bonds[from_bond_idx].begin_atom_idx
            if parent_idx == atom_idx:
                parent_idx = mol.bonds[from_bond_idx].end_atom_idx
        
        tree_edges = []
        ring_closures = []
        ring_openings = []
        
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            if bond_idx == from_bond_idx:
                continue
            
            if bond_idx in ring_bonds_set:
                if nbr_idx in atom_to_id:
                    if atom_to_id[nbr_idx] < atom_string_id:
                        ring_closures.append((ranks[nbr_idx], nbr_idx, bond_idx))
                else:
                    ring_openings.append((ranks[nbr_idx], nbr_idx, bond_idx))
            else:
                if nbr_idx not in atom_to_id:
                    tree_edges.append((ranks[nbr_idx], nbr_idx, bond_idx))
        
        ring_closures.sort()
        tree_edges.sort()
        ring_openings.sort()
        
        # Build DFS neighbor order (SCRIPT Priority): [H] < Parent < Rings < Branches
        ordered_neighbors = []
        if atom.implicit_hs > 0: ordered_neighbors.append(-1)
        if parent_idx != -1: ordered_neighbors.append(parent_idx)
        for _, nbr_idx, _ in ring_closures: ordered_neighbors.append(nbr_idx)
        for _, nbr_idx, _ in ring_openings: ordered_neighbors.append(nbr_idx)
        for _, nbr_idx, _ in tree_edges: ordered_neighbors.append(nbr_idx)
        
        if len(ordered_neighbors) == 4:
            dfs_orders[atom_idx] = ordered_neighbors
        
        # Recurse
        all_edges = tree_edges + ring_openings
        for _, nbr_idx, bond_idx in all_edges:
            if nbr_idx not in atom_to_id:
                self._collect_dfs_neighbor_orders(mol, nbr_idx, atom_to_id, ranks, bond_idx, ring_bonds_set, dfs_orders)

    def _dfs(self, mol: CoreMolecule, atom_idx, atom_to_id, ranks, from_bond_idx, ring_counter, ring_bonds_set, depths):
        """DFS traversal that builds the SCRIPT string."""
        atom_string_id = len(atom_to_id)
        atom_to_id[atom_idx] = atom_string_id
        atom = mol.atoms[atom_idx]
        curr_depth = depths.get(atom_idx, 1)

        parent_idx = -1
        if from_bond_idx >= 0:
            parent_idx = mol.bonds[from_bond_idx].begin_atom_idx
            if parent_idx == atom_idx:
                parent_idx = mol.bonds[from_bond_idx].end_atom_idx

        tree_edges = []
        ring_closures = []  # Rings closing at this atom (back-count notation)
        
        ring_openings = []
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            if bond_idx == from_bond_idx:
                continue
            
            if bond_idx in ring_bonds_set:
                if nbr_idx in atom_to_id:
                    if atom_to_id[nbr_idx] < atom_string_id:
                        # FIXED: Calculate lookback based on flat linear string position, not depth
                        # This is the critical fix for bridged rings
                        target_position = atom_to_id[nbr_idx]
                        current_position = atom_string_id
                        topo_size = current_position - target_position + 1
                        is_arom = mol.bonds[bond_idx].bond_type == BondType.AROMATIC
                        anubandha = ":" if is_arom else "-"
                        bond_sym = self._bond_symbol(mol.bonds[bond_idx], nbr_idx, mol)
                        if is_arom: bond_sym = "" # redundant for &6:
                        ring_closures.append((topo_size, f"{bond_sym}&{topo_size}{anubandha}", nbr_idx))
                    else:
                        # This shouldn't happen with sorted DFS
                        pass
                else:
                    # Ring opening (target not visited yet)
                    ring_openings.append((ranks[nbr_idx], nbr_idx, bond_idx))
            else:
                if nbr_idx not in atom_to_id:
                    tree_edges.append((ranks[nbr_idx], nbr_idx, bond_idx))

        ring_closures.sort()
        tree_edges.sort()
        ring_openings.sort()
        # Build ordered neighbors for stereo perception
        # SCRIPT Priority: [H] < Parent < Ring-Closures < Ring-Openings < Branches < Main-Chain
        ordered_neighbors = []
        ihs = getattr(atom, 'implicit_hs', 0) or 0
        if ihs > 0: ordered_neighbors.append(-1)
        if parent_idx != -1: ordered_neighbors.append(parent_idx)
        for _, _, nbr_idx in ring_closures: ordered_neighbors.append(nbr_idx)
        for _, nbr_idx, _ in ring_openings: ordered_neighbors.append(nbr_idx)
        for _, nbr_idx, _ in tree_edges: ordered_neighbors.append(nbr_idx)

        # Build string parts
        parts = []
        if from_bond_idx >= 0:
            parts.append(self._bond_symbol(mol.bonds[from_bond_idx], atom_idx, mol))

        parts.append(self._atom_string(atom, atom_idx, mol, ranks, ordered_neighbors))
        
        for back_count, bond_sym, _ in ring_closures:
            parts.append(bond_sym)

        # Recursively process tree edges
        if tree_edges:
            # First process all but the last edge as branches
            for i in range(len(tree_edges) - 1):
                _, nbr_idx, bond_idx = tree_edges[i]
                if nbr_idx not in atom_to_id:
                    depths[nbr_idx] = curr_depth + 1
                    branch_result = self._dfs(mol, nbr_idx, atom_to_id, ranks, bond_idx, ring_counter, ring_bonds_set, depths)
                    if branch_result:
                        parts.append('(' + branch_result + ')')
            
            # Process last edge as main chain
            _, nbr_idx, bond_idx = tree_edges[-1]
            if nbr_idx not in atom_to_id:
                depths[nbr_idx] = curr_depth + 1
                parts.append(self._dfs(mol, nbr_idx, atom_to_id, ranks, bond_idx, ring_counter, ring_bonds_set, depths))
            else:
                # If the "main chain" edge was already visited via a branch, 
                # we might need to convert one of the previous branches to the main chain.
                # But for canonicality, if the last edge is gone, we just end the string or find the last sibling that wasn't visited.
                # Actually, the simplest way is to ensure parts is joined correctly.
                pass
        
        return "".join(parts)

    def _atom_string(self, atom, atom_idx, mol, ranks, ordered_neighbors):
        symbol = atom.symbol

        # Wildcard atom — emit bare '*', no brackets needed
        if getattr(atom, 'is_wildcard', False) or atom.atomic_num == 0:
            return '*'

        # Determine chiral symbol based on stereo_type
        stereo_t = getattr(atom, 'stereo_type', StereoType.NONE)
        chiral_sym = ""
        if stereo_t in (StereoType.NONE, StereoType.TETRAHEDRAL):
            if len(ordered_neighbors) == 4:
                chiral_sym = get_chiral_symbol(atom_idx, ordered_neighbors, mol, ranks)
        elif stereo_t == StereoType.SQUARE_PLANAR:
            chiral_sym = "@SP"
        elif stereo_t == StereoType.OCTAHEDRAL:
            chiral_sym = "@OH"
        elif stereo_t == StereoType.ATROPISOMER:
            # Route through get_chiral_symbol so that the geometry-resolved
            # @AX1 / @AX2 parity is used.  Falls back to bare "@AX" if the
            # dihedral could not be computed (no 3D coords available).
            resolved = get_chiral_symbol(atom_idx, ordered_neighbors, mol, ranks)
            chiral_sym = resolved if resolved else "@AX"
        elif stereo_t == StereoType.TRIG_BIPYRAMIDAL:
            chiral_sym = "@TB"
        elif stereo_t == StereoType.PYRAMIDAL:
            chiral_sym = "@PY"

        if (symbol in ORGANIC_SUBSET
                and atom.formal_charge == 0
                and atom.isotope == 0
                and atom.radical_electrons == 0
                and chiral_sym == ""
                and (not hasattr(atom, 'mapping') or atom.mapping == 0)
                and getattr(atom, 'spin', 0) == 0
                and not getattr(atom, 'is_excited', False)
                and self._has_default_valence(atom, mol, atom_idx)):
            return symbol

        parts = ['[']
        if atom.isotope > 0:
            parts.append(str(atom.isotope))
        parts.append(symbol)

        if chiral_sym:
            parts.append(chiral_sym)

        ihs = getattr(atom, 'implicit_hs', 0) or 0
        if ihs > 0:
            parts.append('H')
            if ihs > 1:
                parts.append(str(ihs))

        if atom.formal_charge > 0:
            parts.append('+' + (str(atom.formal_charge) if atom.formal_charge > 1 else ''))
        elif atom.formal_charge < 0:
            parts.append('-' + (str(abs(atom.formal_charge)) if atom.formal_charge < -1 else ''))

        if getattr(atom, 'mapping', 0) > 0:
            parts.append(f':{atom.mapping}')

        parts.append(']')
        return "".join(parts)

    def _bond_symbol(self, bond, to_atom_idx, mol):
        bt = bond.bond_type
        # Typed bonds
        if bt == BondType.DOUBLE:     base = '='
        elif bt == BondType.TRIPLE:   base = '#'
        elif bt == BondType.AROMATIC: base = ':'
        elif bt == BondType.DATIVE:   base = '->'
        elif bt == BondType.REV_DATIVE: base = '<-'
        elif bt == BondType.TAUTOMERIC: base = '=:'
        elif bt == BondType.COORDINATE: base = '>'
        elif bt == BondType.STAR:     base = '*'
        # Legacy int fallback
        elif bt == 2: base = '='
        elif bt == 3: base = '#'
        elif bt == 4: base = ':'
        elif bond.bond_dir != 0:
            is_reverse = (bond.end_atom_idx == to_atom_idx)
            if bond.bond_dir in (1, 3):
                base = '\\' if is_reverse else '/'
            elif bond.bond_dir in (2, 4):
                base = '/' if is_reverse else '\\'
            else:
                base = ''
        else:
            base = ''

        # V3.6: append translation vector for periodic/cross-cell bonds
        t = getattr(bond, 'translation', (0, 0, 0))
        if t != (0, 0, 0):
            # Normalise direction: store as seen from begin_atom_idx → end_atom_idx
            if bond.end_atom_idx == to_atom_idx:
                tx, ty, tz = t
            else:
                tx, ty, tz = -t[0], -t[1], -t[2]
            # Single bonds must be explicit when carrying a translation vector
            # so the output is reparseable (bare '@' after atom is ambiguous)
            if base == '':
                base = '-'
            base = f"{base}@{tx},{ty},{tz}"

        return base

    def _has_default_valence(self, atom, mol, atom_idx):
        if atom.symbol not in DEFAULT_VALENCE: return False

        valence = getattr(atom, 'implicit_hs', 0) or 0
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            bond = mol.bonds[bond_idx]
            bt = bond.bond_type
            if bt in (BondType.DOUBLE, 2):     valence += 2
            elif bt in (BondType.TRIPLE, 3):   valence += 3
            elif bt in (BondType.AROMATIC, 4): valence += 1.5
            elif bt in (BondType.DATIVE, BondType.REV_DATIVE,
                        BondType.COORDINATE, BondType.TAUTOMERIC,
                        BondType.STAR, 5, 6, 7, 8, 9):
                valence += 1
            else:
                valence += 1

        return abs(valence - DEFAULT_VALENCE[atom.symbol]) < 0.1

def canonicalize_mol(mol):
    """Note: This takes a CoreMolecule, not RDKit mol."""
    return SCRIPTCanonicalizer().canonicalize_core(mol)

def canonicalize_SCRIPT(script_string):
    from .parser import SCRIPTParser
    p = SCRIPTParser()
    result = p.parse(script_string)
    if not result["success"]: return None
    return SCRIPTCanonicalizer().canonicalize_core(result["molecule"])
