"""
SCRIPT Validator - State-machine checker that prevents invalid strings
"""

from typing import Optional, Dict, List, Tuple, Set

class SCRIPTValidator:
    """State-machine validator for SCRIPT strings"""
    
    def __init__(self):
        self.valid_elements = self._get_valid_elements()
        self.valid_bonds = {'-', '=', '#', ':', '/', '\\', '->'}
        self.valid_amino_acids = {
            'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
        }
    
    def is_valid(self, script_string: str) -> bool:
        """Check if SCRIPT string is valid"""
        if not script_string:
            return False
        
        # Basic syntax validation
        if not self._check_basic_syntax(script_string):
            return False
        
        # For now, just do basic syntax checks
        # Full parser validation would be done separately
        return True
    
    def validate_with_errors(self, script_string: str) -> Dict[str, any]:
        """Validate and return detailed error information"""
        if not script_string:
            return {"valid": False, "error": "Empty string"}
        
        # Basic syntax check
        syntax_error = self._check_basic_syntax_detailed(script_string)
        if syntax_error:
            return {"valid": False, "error": syntax_error}
        
        return {"valid": True, "error": None}
    
    def _check_basic_syntax(self, script_string: str) -> bool:
        """Basic syntax checks"""
        # Balanced parentheses
        if not self._check_balanced_parens(script_string):
            return False
        
        # Balanced brackets
        if not self._check_balanced_brackets(script_string):
            return False
        
        # Balanced braces (peptide blocks)
        if not self._check_balanced_braces(script_string):
            return False
        
        return True
    
    def _check_basic_syntax_detailed(self, script_string: str) -> Optional[str]:
        """Basic syntax checks with error messages"""
        if not self._check_balanced_parens(script_string):
            return "Unbalanced parentheses"
        
        if not self._check_balanced_brackets(script_string):
            return "Unbalanced brackets"
        
        if not self._check_balanced_braces(script_string):
            return "Unbalanced braces"
        
        return None
    
    def _check_balanced_parens(self, s: str) -> bool:
        """Check balanced parentheses"""
        count = 0
        for char in s:
            if char == '(':
                count += 1
            elif char == ')':
                count -= 1
                if count < 0:
                    return False
        return count == 0
    
    def _check_balanced_brackets(self, s: str) -> bool:
        """Check balanced brackets"""
        count = 0
        for char in s:
            if char == '[':
                count += 1
            elif char == ']':
                count -= 1
                if count < 0:
                    return False
        return count == 0
    
    def _check_balanced_braces(self, s: str) -> bool:
        """Check balanced braces"""
        count = 0
        for char in s:
            if char == '{':
                count += 1
            elif char == '}':
                count -= 1
                if count < 0:
                    return False
        return count == 0
    
    def _check_chemical_validity(self, parse_result: Dict) -> bool:
        """Check chemical validity of parsed structure"""
        # Basic checks - can be expanded
        atoms = parse_result.get("atoms", [])
        
        # Must have at least one atom
        if not atoms:
            return False
        
        # Check valid elements
        for atom in atoms:
            element = atom.get("element", "")
            if element and not self._is_valid_element(element):
                return False
        
        return True
    
    def _check_chemical_validity_detailed(self, parse_result: Dict) -> Optional[str]:
        """Check chemical validity with error messages"""
        atoms = parse_result.get("atoms", [])
        
        if not atoms:
            return "No atoms found"
        
        for atom in atoms:
            element = atom.get("element", "")
            if element and not self._is_valid_element(element):
                return f"Invalid element: {element}"
        
        return None
    
    def _is_valid_element(self, element: str) -> bool:
        """Check if element symbol is valid"""
        return element in self.valid_elements
    
    def _get_valid_elements(self) -> Set[str]:
        """Get set of valid element symbols"""
        # Periodic table elements (simplified)
        elements = {
            'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
            'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
        }
        return elements

# Convenience function
def is_valid_SCRIPT(script_string: str) -> bool:
    """Check if SCRIPT string is valid"""
    validator = SCRIPTValidator()
    return validator.is_valid(script_string)