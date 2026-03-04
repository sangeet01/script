"""
SCRIPT Peptide Handler - Handle {A.G.S...} macro mode
"""

from typing import Dict, List, Optional, Tuple

# Standard amino acid codes to SMILES mapping
AMINO_ACID_SMILES = {
    'A': 'N[C@@H](C)C(=O)O',           # Alanine
    'R': 'N[C@@H](CCCNC(=N)N)C(=O)O', # Arginine
    'N': 'N[C@@H](CC(=O)N)C(=O)O',    # Asparagine
    'D': 'N[C@@H](CC(=O)O)C(=O)O',    # Aspartic acid
    'C': 'N[C@@H](CS)C(=O)O',          # Cysteine
    'E': 'N[C@@H](CCC(=O)O)C(=O)O',   # Glutamic acid
    'Q': 'N[C@@H](CCC(=O)N)C(=O)O',   # Glutamine
    'G': 'NCC(=O)O',                   # Glycine
    'H': 'N[C@@H](CC1=CNC=N1)C(=O)O', # Histidine
    'I': 'N[C@@H]([C@@H](C)CC)C(=O)O',# Isoleucine
    'L': 'N[C@@H](CC(C)C)C(=O)O',     # Leucine
    'K': 'N[C@@H](CCCCN)C(=O)O',      # Lysine
    'M': 'N[C@@H](CCSC)C(=O)O',       # Methionine
    'F': 'N[C@@H](CC1=CC=CC=C1)C(=O)O', # Phenylalanine
    'P': 'N1[C@@H](CCC1)C(=O)O',      # Proline
    'S': 'N[C@@H](CO)C(=O)O',         # Serine
    'T': 'N[C@@H]([C@H](C)O)C(=O)O',  # Threonine
    'W': 'N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O', # Tryptophan
    'Y': 'N[C@@H](CC1=CC=C(O)C=C1)C(=O)O',   # Tyrosine
    'V': 'N[C@@H](C(C)C)C(=O)O'       # Valine
}

# Common post-translational modifications
PTM_SMILES = {
    'pS': 'N[C@@H](COP(=O)(O)O)C(=O)O',     # Phosphoserine
    'pT': 'N[C@@H]([C@H](C)OP(=O)(O)O)C(=O)O', # Phosphothreonine
    'pY': 'N[C@@H](CC1=CC=C(OP(=O)(O)O)C=C1)C(=O)O', # Phosphotyrosine
    'mK': 'N[C@@H](CCCCNC)C(=O)O',          # Methylated lysine
}

class PeptideHandler:
    """Handle peptide macro mode {A.G.S...}"""
    
    def __init__(self):
        self.amino_acids = AMINO_ACID_SMILES
        self.ptms = PTM_SMILES
    
    def expand_peptide_block(self, peptide_sequence: List[str]) -> Optional[str]:
        """Convert peptide sequence to full atomic SMILES"""
        if not peptide_sequence:
            return None
        
        smiles_parts = []
        
        for i, residue in enumerate(peptide_sequence):
            # Get SMILES for this residue
            if residue in self.amino_acids:
                aa_smiles = self.amino_acids[residue]
            elif residue in self.ptms:
                aa_smiles = self.ptms[residue]
            else:
                return None  # Unknown residue
            
            # Modify for peptide bond formation
            if i == 0:
                # First residue - keep N-terminus
                smiles_parts.append(aa_smiles)
            else:
                # Remove N-terminus and form peptide bond
                # This is simplified - full implementation would need proper bond formation
                modified_smiles = self._remove_n_terminus(aa_smiles)
                smiles_parts.append(modified_smiles)
        
        # Join with peptide bonds (simplified)
        return self._join_with_peptide_bonds(smiles_parts)
    
    def _remove_n_terminus(self, smiles: str) -> str:
        """Remove N-terminus for peptide bond formation"""
        # Simplified - remove leading N
        if smiles.startswith('N[C@@H]'):
            return smiles[1:]  # Remove N
        elif smiles.startswith('NCC'):
            return smiles[1:]  # Remove N (glycine case)
        elif smiles.startswith('N1[C@@H]'):
            return smiles[2:]  # Remove N1 (proline case)
        return smiles
    
    def _join_with_peptide_bonds(self, smiles_parts: List[str]) -> str:
        """Join amino acid SMILES with peptide bonds"""
        # Simplified joining - full implementation would need proper chemistry
        return '.'.join(smiles_parts)
    
    def parse_peptide_notation(self, peptide_string: str) -> Optional[List[str]]:
        """Parse {A.G.S.C[A]K} notation into residue list"""
        # Remove braces
        if not (peptide_string.startswith('{') and peptide_string.endswith('}')):
            return None
        
        inner = peptide_string[1:-1]
        
        # Split by dots
        residues = []
        i = 0
        while i < len(inner):
            if inner[i] == '.':
                i += 1
                continue
            
            # Check for PTM (multi-character)
            if i + 1 < len(inner) and inner[i:i+2] in self.ptms:
                residues.append(inner[i:i+2])
                i += 2
            # Check for bridge notation [A]
            elif inner[i] == '[':
                end = inner.find(']', i)
                if end != -1:
                    residues.append(inner[i:end+1])
                    i = end + 1
                else:
                    return None  # Malformed
            # Single character amino acid
            elif inner[i] in self.amino_acids:
                residues.append(inner[i])
                i += 1
            else:
                return None  # Unknown residue
        
        return residues
    
    def is_valid_peptide_sequence(self, sequence: List[str]) -> bool:
        """Check if peptide sequence is valid"""
        for residue in sequence:
            if residue.startswith('[') and residue.endswith(']'):
                # Bridge notation - assume valid for now
                continue
            elif residue not in self.amino_acids and residue not in self.ptms:
                return False
        return True

# Convenience functions
def expand_peptide(peptide_string: str) -> Optional[str]:
    """Expand peptide notation to full SMILES"""
    handler = PeptideHandler()
    sequence = handler.parse_peptide_notation(peptide_string)
    if sequence is None:
        return None
    return handler.expand_peptide_block(sequence)

def is_valid_peptide(peptide_string: str) -> bool:
    """Check if peptide notation is valid"""
    handler = PeptideHandler()
    sequence = handler.parse_peptide_notation(peptide_string)
    if sequence is None:
        return False
    return handler.is_valid_peptide_sequence(sequence)