"""
SCRIPT - Structural Chemical Representation in Plain Text

A next-generation chemical line notation that addresses fundamental 
limitations in SMILES while maintaining human readability.
"""

from .parser import SCRIPTParser
from .canonical import canonicalize_SCRIPT
from .validator import is_valid_SCRIPT
from .peptide import expand_peptide_to_script, is_valid_monomer
from .mol import CoreMolecule, CoreBond, CoreAtom, BondType, StereoType, PolymerBlock

# Convenience function
def parse_script(script_string: str):
    """Parse SCRIPT string to internal representation"""
    parser = SCRIPTParser()
    return parser.parse(script_string)

# Optional RDKit bridge — only loaded on demand
def get_rdkit_bridge():
    """Lazy import of RDKit bridge. Only loads RDKit when called."""
    try:
        from .rdkit_bridge import from_rdkit, to_rdkit, script_to_smiles, smiles_to_script
        return from_rdkit, to_rdkit, script_to_smiles, smiles_to_script
    except ImportError:
        raise ImportError("RDKit bridge not available. Install: pip install rdkit")

def get_visualizer():
    """Lazy import of visualizer. Only loads RDKit when called."""
    try:
        from .visualizer import draw_molecule
        return draw_molecule
    except ImportError:
        raise ImportError("Visualizer not available. Install: pip install rdkit")

# Constrained decoder for ML generation (100% validity guarantee)
def get_constrained_decoder():
    """Lazy import of ConstrainedSCRIPTDecoder for ML generation.
    
    Returns the ConstrainedSCRIPTDecoder class, which enforces grammar-state-
    aware token masking during autoregressive generation. This gives ~100%
    validity for ML-generated SCRIPT strings without changing the notation.
    
    Usage:
        decoder = get_constrained_decoder()()
        decoder.reset()
        mask = decoder.get_valid_token_mask(vocab)
    """
    from .constrained_decoder import ConstrainedSCRIPTDecoder
    return ConstrainedSCRIPTDecoder

__version__ = "3.0.0"
__author__ = "SCRIPT Development Team"