"""
SCRIPT - Structural Chemical Representation in Plain Text

A next-generation chemical line notation that addresses fundamental 
limitations in SMILES while maintaining human readability.
"""

from .parser import SCRIPTParser, parse_script
from .canonical import canonicalize_SCRIPT, SCRIPTCanonicalizer
from .validator import is_valid_SCRIPT
from .peptide import expand_peptide_to_script, is_valid_monomer
from .mol import CoreMolecule, CoreBond, CoreAtom, BondType, StereoType, PolymerBlock

def _make_rdkit_stub(fn_name: str):
    """Return a stub that raises a clear ImportError naming the missing function."""
    def _stub(*args, **kwargs):
        raise ImportError(
            f"script.{fn_name}() requires RDKit. "
            "Install it with:  pip install rdkit"
        )
    _stub.__name__ = fn_name
    return _stub

# RDKit-dependent names — real callables when RDKit is present, clear stubs otherwise.
# All are hoisted to the top level so the quickstart works without sub-module imports.
try:
    from .rdkit_bridge import (
        smiles_to_script as smiles_to_script,
        script_to_smiles as script_to_smiles,
        from_rdkit       as from_rdkit,
        to_rdkit         as to_rdkit,
        SCRIPTFromMol    as SCRIPTFromMol,
        MolFromSCRIPT    as MolFromSCRIPT,
    )
    _RDKIT_AVAILABLE = True
except ImportError:
    smiles_to_script = _make_rdkit_stub("smiles_to_script")
    script_to_smiles = _make_rdkit_stub("script_to_smiles")
    from_rdkit       = _make_rdkit_stub("from_rdkit")
    to_rdkit         = _make_rdkit_stub("to_rdkit")
    SCRIPTFromMol    = _make_rdkit_stub("SCRIPTFromMol")
    MolFromSCRIPT    = _make_rdkit_stub("MolFromSCRIPT")
    _RDKIT_AVAILABLE = False

try:
    from .visualizer import draw_molecule as draw_molecule
except ImportError:
    draw_molecule = _make_rdkit_stub("draw_molecule")

# Legacy shims — kept for backward compatibility only.
def get_rdkit_bridge():
    """Deprecated helper. Use script.smiles_to_script / script.from_rdkit directly."""
    if not _RDKIT_AVAILABLE:
        raise ImportError("RDKit bridge not available. Install: pip install rdkit")
    return from_rdkit, to_rdkit, script_to_smiles, smiles_to_script

def get_visualizer():
    """Deprecated helper. Use script.draw_molecule directly."""
    if _RDKIT_AVAILABLE:
        return draw_molecule
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