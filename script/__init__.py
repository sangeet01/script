"""
SCRIPT - Structural Chemical Representation in Plain Text

A next-generation chemical line notation that addresses fundamental 
limitations in SMILES while maintaining human readability.
"""

from .parser import SCRIPTParser
from .canonical import canonicalize_SCRIPT
from .validator import is_valid_SCRIPT
from .peptide import expand_peptide, is_valid_peptide

# Convenience function
def parse_script(script_string: str):
    """Parse SCRIPT string to internal representation"""
    parser = SCRIPTParser()
    return parser.parse(script_string)

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

__version__ = "2.0.0"
__author__ = "SCRIPT Development Team"