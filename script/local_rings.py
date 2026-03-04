"""
SCRIPT Local Ring Handler - Core differentiating feature
Implements local ring encoding: digit after atom = atoms back to connect
"""

import re
from typing import Dict, List, Tuple

class LocalRingHandler:
    """Handle SCRIPT local ring closures"""
    
    def __init__(self):
        self.ring_connections = {}  # {ring_id: (start_atom, end_atom, bond)}
    
    def parse_local_rings(self, script_string: str) -> Tuple[str, Dict[int, Tuple[int, int, str]]]:
        """Parse local rings and return cleaned string + ring connections"""
        self.ring_connections = {}
        
        # 1. Identify all atom positions in the string
        atom_info = self._find_atom_positions(script_string)
        if not atom_info:
            return script_string, {}
        
        # 2. Extract ring closures for each atom
        # SCRIPT: atom [optional_bond] digit [optional_bond] digit ...
        connections = []
        
        # We'll use a regex that matches the ring closure part specifically
        # (optional bond + digit)
        ring_pattern = re.compile(r'([-=#:/\\.]?|->|<-)(\d+|%\d\d)')
        
        cleaned_chars = list(script_string)
        # Track characters to remove (indices)
        to_remove = set()
        
        for atom_idx, (start, end) in enumerate(atom_info):
            # Look ahead from 'end' for multiple ring closures
            current_pos = end
            while current_pos < len(script_string):
                match = ring_pattern.match(script_string, current_pos)
                if not match:
                    break
                
                bond_part = match.group(1)
                ring_part = match.group(2)
                
                # Parse back-count
                if ring_part.startswith('%'):
                    back_count = int(ring_part[1:])
                else:
                    back_count = int(ring_part)
                
                # Calculate target atom (back_count atoms back from current)
                target_idx = atom_idx - back_count
                if target_idx >= 0:
                    connections.append((target_idx, atom_idx, bond_part))
                
                # Mark for removal
                for i in range(match.start(), match.end()):
                    to_remove.add(i)
                
                current_pos = match.end()

        # 3. Build cleaned string (removing digits, keeping bond symbols if needed)
        # Wait, SCRIPT ring closures *include* the bond. 
        # C=6 means a double bond to the atom 6 steps back.
        # However, for SMILES conversion, we need to strip these and re-insert them 
        # as SMILES labels.
        
        # Re-build connections dictionary
        result_connections = {}
        for i, conn in enumerate(connections):
            result_connections[i] = conn
            
        final_cleaned = "".join([c for idx, c in enumerate(script_string) if idx not in to_remove])
        return final_cleaned, result_connections

    def convert_to_smiles_rings(self, script_string: str) -> str:
        """Convert SCRIPT local rings to SMILES ring notation"""
        cleaned, connections = self.parse_local_rings(script_string)
        if not connections:
            return script_string
            
        atom_info = self._find_atom_positions(cleaned)
        insertions = []
        ring_label_counter = 1
        
        for conn_id, (target_atom, current_atom, bond) in connections.items():
            label = str(ring_label_counter)
            if ring_label_counter >= 10:
                label = "%" + label
            ring_label_counter += 1
            
            # Start atom (target)
            if target_atom < len(atom_info):
                insertions.append((atom_info[target_atom][1], label))
            
            # End atom (current)
            if current_atom < len(atom_info):
                # SMILES: atom + bond + label
                # In cleaned string, the bond from the ring closure was removed.
                # We should re-add it if it was provided.
                insertions.append((atom_info[current_atom][1], bond + label))
                
        insertions.sort(key=lambda x: x[0], reverse=True)
        
        result = list(cleaned)
        for pos, text in insertions:
            result.insert(pos, text)
            
        return "".join(result)

    def _find_atom_positions(self, script_string: str) -> List[Tuple[int, int]]:
        """Find start and end positions of each atom"""
        positions = []
        i = 0
        while i < len(script_string):
            if script_string[i] == '[':
                start = i
                end = script_string.find(']', i)
                if end != -1:
                    positions.append((start, end + 1))
                    i = end + 1
                else:
                    i += 1
            elif script_string[i].isupper():
                start = i
                i += 1
                if i < len(script_string) and script_string[i].islower():
                    i += 1
                positions.append((start, i))
            else:
                i += 1
        return positions

def parse_local_rings(script_string: str) -> Tuple[str, Dict[int, Tuple[int, int, str]]]:
    handler = LocalRingHandler()
    return handler.parse_local_rings(script_string)

def convert_local_to_smiles_rings(script_string: str) -> str:
    handler = LocalRingHandler()
    return handler.convert_to_smiles_rings(script_string)