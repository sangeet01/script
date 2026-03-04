"""
SCRIPT - Structural Chemical Representation in Plain Text

A next-generation chemical line notation that addresses fundamental 
limitations in SMILES while maintaining human readability.
"""

from .parser import SCRIPTParser
from .canonical import canonicalize_SCRIPT

# Convenience functions
def parse_script(script_string: str):
    """Parse SCRIPT string to internal representation"""
    parser = SCRIPTParser()
    return parser.parse(script_string)

def is_valid_SCRIPT(script_string: str) -> bool:
    """Check if SCRIPT string is valid"""
    result = parse_script(script_string)
    return result["success"]

# RDKit integration (optional)
try:
    from .rdkit_bridge import (
        script_to_smiles, 
        smiles_to_script, 
        MolFromSCRIPT, 
        SCRIPTFromMol,
        from_rdkit,
        CoreToRDKit
    )
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

__version__ = "0.1.0"
__author__ = "SCRIPT Development Team"
