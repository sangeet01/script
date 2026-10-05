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
            if not data:
                return ""

            # Check for reaction (list of lists) -- inspect actual type, not duck-type.
            first = data[0]
            if isinstance(first, list):
                sides = []
                for side in data:
                    sides.append(self.to_script(side, separator))
                return " >> ".join(sides)

            # Multi-component
            parts = [self.canonicalizer.canonicalize_core(m) for m in data]
            return separator.join(parts)

        return str(data)


def write_core_to_script(mol: Any) -> str:
    """
    Convenience function: convert a CoreMolecule to a canonical SCRIPT string.

    For RDKit Mol objects use script.rdkit_bridge.SCRIPTFromMol instead --
    that function handles kekulization, aromaticity perception, and stereo
    transfer from the RDKit graph before conversion.
    """
    writer = SCRIPTWriter()
    return writer.to_script(mol)