from typing import List, Dict, Optional, Tuple
from .mol import CoreMolecule
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

        # 1. Ranking (Universal)
        rank_map = calculate_ranks(mol)
        num_atoms = len(mol.atoms)
        ranks = [rank_map.get(i, 0) for i in range(num_atoms)]
        
        # 2. DFS from lowest rank atom
        start_indices = [i for i, r in enumerate(ranks) if r == 0]
        start_atom = start_indices[0] if start_indices else 0
        
        # First pass: identify ring bonds and collect DFS neighbor orders
        visited = set()
        ring_bonds_set = set()
        self._identify_rings(mol, start_atom, visited, ring_bonds_set, -1, ranks)
        
        # Second pass: collect DFS neighbor orders for stereochemistry
        atom_to_id = {}
        dfs_neighbor_orders = {}
        self._collect_dfs_neighbor_orders(mol, start_atom, atom_to_id, ranks, -1, ring_bonds_set, dfs_neighbor_orders)
        
        # 3. Perceive chirality using DFS neighbor orders as reference
        perceive_chirality(mol, ranks, dfs_neighbor_orders)
        
        # 4. Build canonical string
        atom_to_id = {}
        ring_bonds = {}
        ring_counter = [1]
        
        result = self._dfs(mol, start_atom, atom_to_id, ranks, -1, ring_bonds, ring_counter, ring_bonds_set)
        
        return result
    
    def _identify_rings(self, mol: CoreMolecule, atom_idx, visited, ring_bonds, from_bond_idx, ranks):
        """First pass: identify which bonds are ring closures."""
        visited.add(atom_idx)
        
        neighbors = []
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            if bond_idx != from_bond_idx:
                neighbors.append((ranks[nbr_idx], nbr_idx, bond_idx))
                
        neighbors.sort(key=lambda x: (x[0], x[1]))
        
        for _, nbr_idx, bond_idx in neighbors:
            if nbr_idx in visited:
                ring_bonds.add(bond_idx)
            else:
                self._identify_rings(mol, nbr_idx, visited, ring_bonds, bond_idx, ranks)
    
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
        
        # Build DFS neighbor order (SCRIPT Priority)
        ordered_neighbors = []
        if parent_idx != -1: ordered_neighbors.append(parent_idx)
        if atom.implicit_hs > 0: ordered_neighbors.append(-1)
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

    def _dfs(self, mol: CoreMolecule, atom_idx, atom_to_id, ranks, from_bond_idx, ring_bonds, ring_counter, ring_bonds_set):
        """DFS traversal that builds the SCRIPT string."""
        atom_string_id = len(atom_to_id)
        atom_to_id[atom_idx] = atom_string_id
        atom = mol.atoms[atom_idx]

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
                        # Back-count closure
                        back_count = atom_string_id - atom_to_id[nbr_idx] + 1
                        bond_sym = self._bond_symbol(mol.bonds[bond_idx], nbr_idx)
                        ring_closures.append((back_count, bond_sym, nbr_idx))
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
        # SCRIPT Priority: Parent < [H] < Ring-Closures < Ring-Openings < Branches < Main-Chain
        ordered_neighbors = []
        
        if parent_idx != -1: ordered_neighbors.append(parent_idx)
        if atom.implicit_hs > 0: ordered_neighbors.append(-1)
        for _, _, nbr_idx in ring_closures: ordered_neighbors.append(nbr_idx)
        for _, nbr_idx, _ in ring_openings: ordered_neighbors.append(nbr_idx)
        for _, nbr_idx, _ in tree_edges: ordered_neighbors.append(nbr_idx)

        # Build string parts
        parts = []
        if from_bond_idx >= 0:
            parts.append(self._bond_symbol(mol.bonds[from_bond_idx], atom_idx))

        parts.append(self._atom_string(atom, atom_idx, mol, ranks, ordered_neighbors))
        
        for back_count, bond_sym, _ in ring_closures:
            if back_count <= 9: parts.append(bond_sym + str(back_count))
            else: parts.append(bond_sym + '%' + str(back_count).zfill(2))

        # Recursively process tree edges
        if tree_edges:
            # First process all but the last edge as branches
            for i in range(len(tree_edges) - 1):
                _, nbr_idx, bond_idx = tree_edges[i]
                if nbr_idx not in atom_to_id:
                    branch_str = self._dfs(mol, nbr_idx, atom_to_id, ranks, bond_idx, ring_bonds, ring_counter, ring_bonds_set)
                    parts.append('(' + branch_str + ')')
            
            # Process last edge as main chain
            _, nbr_idx, bond_idx = tree_edges[-1]
            if nbr_idx not in atom_to_id:
                parts.append(self._dfs(mol, nbr_idx, atom_to_id, ranks, bond_idx, ring_bonds, ring_counter, ring_bonds_set))
            else:
                # If the "main chain" edge was already visited via a branch, 
                # we might need to convert one of the previous branches to the main chain.
                # But for canonicality, if the last edge is gone, we just end the string or find the last sibling that wasn't visited.
                # Actually, the simplest way is to ensure parts is joined correctly.
                pass
        
        return "".join(parts)

    def _atom_string(self, atom, atom_idx, mol, ranks, ordered_neighbors):
        symbol = atom.symbol
        chiral_sym = ""
        if len(ordered_neighbors) == 4:
            chiral_sym = get_chiral_symbol(atom_idx, ordered_neighbors, mol, ranks)

        if (symbol in ORGANIC_SUBSET
                and atom.formal_charge == 0
                and atom.isotope == 0
                and atom.radical_electrons == 0
                and chiral_sym == ""
                and self._has_default_valence(atom, mol, atom_idx)):
            return symbol

        parts = ['[']
        if atom.isotope > 0:
            parts.append(str(atom.isotope))
        parts.append(symbol)
        
        if chiral_sym:
            parts.append(chiral_sym)

        if atom.implicit_hs > 0:
            parts.append('H')
            if atom.implicit_hs > 1:
                parts.append(str(atom.implicit_hs))
        
        if atom.formal_charge > 0:
            parts.append('+' + (str(atom.formal_charge) if atom.formal_charge > 1 else ''))
        elif atom.formal_charge < 0:
            parts.append('-' + (str(abs(atom.formal_charge)) if atom.formal_charge < -1 else ''))
            
        parts.append(']')
        return "".join(parts)

    def _bond_symbol(self, bond, to_atom_idx):
        bt = bond.bond_type
        if bt == 2: return '='
        if bt == 3: return '#'
        if bt == 4: return ':' # SCRIPT uses : for aromatic if not implied
        
        if bond.bond_dir != 0:
            is_reverse = (bond.end_atom_idx == to_atom_idx)
            # BEGINWEDGE (1) / BEGINDASH (2) / ENDDOWNRIGHT (3) / ENDUPRIGHT (4)
            if bond.bond_dir in (1, 3): 
                return '\\' if is_reverse else '/'
            elif bond.bond_dir in (2, 4):
                return '/' if is_reverse else '\\'
        return ""

    def _has_default_valence(self, atom, mol, atom_idx):
        if atom.symbol not in DEFAULT_VALENCE: return False
        valence = atom.implicit_hs
        for nbr_idx, bond_idx in mol.adj.get(atom_idx, []):
            bt = mol.bonds[bond_idx].bond_type
            if bt == 2: valence += 2
            elif bt == 3: valence += 3
            elif bt == 4: valence += 1.5
            else: valence += 1
        return valence == DEFAULT_VALENCE[atom.symbol]

def canonicalize_mol(mol):
    """Note: This takes a CoreMolecule, not RDKit mol."""
    return SCRIPTCanonicalizer().canonicalize_core(mol)

def canonicalize_SCRIPT(script_string):
    from .parser import SCRIPTParser
    p = SCRIPTParser()
    result = p.parse(script_string)
    if not result["success"]: return None
    return SCRIPTCanonicalizer().canonicalize_core(result["molecule"])
