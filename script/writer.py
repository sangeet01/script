"""
SCRIPT Writer - Convert RDKit molecules to canonical SCRIPT strings
"""

from typing import Optional, Dict, List, Set
from .canonical import SCRIPTCanonicalizer

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

class SCRIPTWriter:
    """Convert RDKit molecules to canonical SCRIPT strings"""
    
    def __init__(self):
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit is required")
        self.canonicalizer = SCRIPTCanonicalizer()
    
    def mol_to_script(self, mol: Chem.Mol) -> Optional[str]:
        """Convert RDKit molecule to canonical SCRIPT"""
        if mol is None:
            return None
        
        try:
            # Use canonicalizer for consistent output
            return self.canonicalizer.canonicalize_mol(mol)
        except Exception:
            return None

def SCRIPTFromMol(mol: Chem.Mol) -> Optional[str]:
    """Convert RDKit molecule to canonical SCRIPT string"""
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required")
    
    writer = SCRIPTWriter()
    return writer.mol_to_script(mol)