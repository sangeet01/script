from typing import Optional, Union, List, Any
from .mol import CoreMolecule
from .canonical import SCRIPTCanonicalizer

class SCRIPTWriter:
    """
    High-Level SCRIPT Writer (RDKit-Free)
    Converts CoreMolecule objects into canonical SCRIPT strings.
    """
    
    def __init__(self):
        self.canonicalizer = SCRIPTCanonicalizer()
    
    def to_script(self, data: Union[CoreMolecule, List[CoreMolecule], List[List[CoreMolecule]]], 
                  separator: str = ".") -> str:
        """
        Convert molecular data to SCRIPT.
        - CoreMolecule: Single component
        - List[CoreMolecule]: Multi-component (salts/solvates)
        - List[List[CoreMolecule]]: Reaction (Reactants >> Products)
        """
        if isinstance(data, CoreMolecule):
            return self.canonicalizer.canonicalize_core(data)
            
        if isinstance(data, list):
            if not data: return ""
            
            # Check for reaction (list of lists)
            if isinstance(data[0], list):
                sides = []
                for side in data:
                    sides.append(self.to_script(side, separator))
                return " >> ".join(sides)
            
            # Multi-component
            parts = [self.canonicalizer.canonicalize_core(m) for m in data]
            return separator.join(parts)
            
        return str(data)

def SCRIPTFromMol(mol: Any) -> str:
    """
    Convenience function. 
    If mol is CoreMolecule, writes SCRIPT.
    If mol is RDKit mol (and RDKit installed), converts via bridge (TBD) or errors.
    """
    writer = SCRIPTWriter()
    return writer.to_script(mol)