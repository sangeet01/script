from typing import List, Dict, Optional, Tuple
from .mol import CoreMolecule, BondType, StereoType
from .ranking import calculate_ranks
from .stereo import get_chiral_symbol, perceive_chirality

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
            return self._canonicalize_component(mol, list(range(len(mol.atoms))))
        
        # Canonicalize each component independently, sort for determinism
        fragment_strings = []
        for component in components:
            s = self._canonicalize_component(mol, component)
            if s:
                fragment_strings.append(s)
        
        # Strip trailing aliphatic anubandha '.' from each fragment before joining.
        # The grammar now accepts &N without trailing '.', so this is safe.
        # Aromatic anubandha ':' is never at end of a fragment so no need to strip it.
        fragment_strings = [s.rstrip('.') for s in fragment_strings]

        # Sort fragments canonically (largest first, then lexicographic)
        fragment_strings.sort(key=lambda s: (-len(s), s))
        return '.'.join(fragment_strings)

    def _find_components(self, mol: CoreMolecule) -> List[List[int]]:
        """Find all connected components via BFS."""
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
                for nbr, _ in mol.adj.get(node, []):
                    if nbr not in visited:
                        visited.add(nbr)
                        queue.append(nbr)
            components.append(component)
        return components

    def _canonicalize_component(self, mol: CoreMolecule, atom_indices: List[int]) -> Optional[str]:
        """Canonicalize a single connected component."""
        # 1. Ranking (Universal)
        rank_map = calculate_ranks(mol)
        num_atoms = len(mol.atoms)
        ranks = [rank_map.get(i, 0) for i in range(num_atoms)]
        
        # 2. DFS from lowest rank atom within this component
        component_set = set(atom_indices)
        start_candidates = [(ranks[i], i) for i in atom_indices if ranks[i] == min(ranks[j] for j in atom_indices)]
        start_atom = start_candidates[0][1] if start_candidates else atom_indices[0]
        
        # First pass: identify ring bonds
        visited = set()
        ring_bonds_set = set()
        self._find_ring_bonds(mol, start_atom, visited, -1, ring_bonds_set)
        
        # Second pass: collect DFS neighbor orders for stereochemistry
        atom_to_id_collect = {}
        dfs_neighbor_orders = {}
        self._collect_dfs_neighbor_orders(mol, start_atom, atom_to_id_collect, ranks, -1, ring_bonds_set, dfs_neighbor_orders)
        
        # 3. Perceive chirality (skip if already resolved natively)
        has_chiral_data = hasattr(mol, 'chiral_centers') and mol.chiral_centers
        if not has_chiral_data:
            perceive_chirality(mol, ranks, dfs_neighbor_orders)
        
        # 4. Build canonical string
        atom_to_id = {}
        depths = {start_atom: 1}
        ring_counter = [0]
        
        result = self._dfs(mol, start_atom, atom_to_id, ranks, -1, ring_counter, ring_bonds_set, depths)
        
        return result
    
    def _find_ring_bonds(self, mol: CoreMolecule, atom_idx, visited, from_bond_idx, ring_bonds):
        """First pass: identify which bonds are ring closures."""
        visited.add(atom_idx)
        
        neighbors = []
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            if bond_idx != from_bond_idx:
                neighbors.append((mol.atoms[nbr_idx].rank, nbr_idx, bond_idx)) # Use atom rank for sorting
                
        neighbors.sort(key=lambda x: (x[0], x[1]))
        
        for _, nbr_idx, bond_idx in neighbors:
            if nbr_idx in visited:
                ring_bonds.add(bond_idx)
            else:
                self._find_ring_bonds(mol, nbr_idx, visited, bond_idx, ring_bonds)
    
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
                        # Topological distance is the difference in depths
                        topo_size = curr_depth - depths.get(nbr_idx, 1) + 1
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
            chiral_sym = "@AX"
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
        if bt == BondType.DOUBLE:     return '='
        if bt == BondType.TRIPLE:     return '#'
        if bt == BondType.AROMATIC:   return ':'
        if bt == BondType.DATIVE:     return '->'
        if bt == BondType.REV_DATIVE: return '<-'
        if bt == BondType.TAUTOMERIC: return '=:'
        if bt == BondType.COORDINATE: return '>'
        if bt == BondType.STAR:       return '*'
        # Legacy int fallback
        if bt == 2: return '='
        if bt == 3: return '#'
        if bt == 4: return ':'

        if bond.bond_dir != 0:
            is_reverse = (bond.end_atom_idx == to_atom_idx)
            if bond.bond_dir in (1, 3):
                return '\\' if is_reverse else '/'
            elif bond.bond_dir in (2, 4):
                return '/' if is_reverse else '\\'
        return ""

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
