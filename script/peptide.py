"""
SCRIPT Biopolymer Handler - RDKit-Free Implementation
Handles peptide {A.G.S} and nucleic acid {dA.dG.dC.dT} macro notation.
Translates monomer codes into CoreMolecule atoms via native SCRIPT parsing.
"""

from typing import Dict, List, Optional, Tuple, Any
from lark import Tree, Token

# Amino acid residue SCRIPT notation (native SCRIPT V2 strings).
# Aromatic amino acids (F, W, Y, H and aliases) use V2 `:` aromatic bond
# notation instead of SMILES lowercase atoms (c1...).  Converted from
# canonical SMILES via the RDKit bridge + canonicalize_SCRIPT.
AMINO_ACID_SCRIPT: Dict[str, str] = {
    'A':  'N[C@@H](C)C(=O)O',
    'R':  'N[C@@H](CCCNC(=N)N)C(=O)O',
    'N':  'N[C@@H](CC(=O)N)C(=O)O',
    'D':  'N[C@@H](CC(=O)O)C(=O)O',
    'C':  'N[C@@H](CS)C(=O)O',
    'E':  'N[C@@H](CCC(=O)O)C(=O)O',
    'Q':  'N[C@@H](CCC(=O)N)C(=O)O',
    'G':  'NCC(=O)O',
    'H':  'N:C:[NH]:C:C&5:C[C@RH](N)C(O)=O',
    'I':  'N[C@@H]([C@@H](C)CC)C(=O)O',
    'L':  'N[C@@H](CC(C)C)C(=O)O',
    'K':  'N[C@@H](CCCCN)C(=O)O',
    'M':  'N[C@@H](CCSC)C(=O)O',
    'F':  'C:C(:C:C:C:C&6:)C[C@RH](C(=O)O)N',
    'P':  'N1[C@@H](CCC1)C(=O)O',
    'S':  'N[C@@H](CO)C(=O)O',
    'T':  'N[C@@H]([C@H](C)O)C(=O)O',
    'W':  'C:C:C:C:[C]:C(:C:[NH]:[C]&5:&9:)C[C@RH](C(O)=O)N',
    'Y':  'O=C([C@RH](CC:C:C:C(:C:C&6:)O)N)O',
    'V':  'N[C@@H](C(C)C)C(=O)O',
    # 3-letter codes (alias)
    'Ala': 'N[C@@H](C)C(=O)O',
    'Gly': 'NCC(=O)O',
    'Val': 'N[C@@H](C(C)C)C(=O)O',
    'Leu': 'N[C@@H](CC(C)C)C(=O)O',
    'Ile': 'N[C@@H]([C@@H](C)CC)C(=O)O',
    'Pro': 'N1[C@@H](CCC1)C(=O)O',
    'Phe': 'C:C(:C:C:C:C&6:)C[C@RH](C(=O)O)N',
    'Trp': 'C:C:C:C:[C]:C(:C:[NH]:[C]&5:&9:)C[C@RH](C(O)=O)N',
    'Met': 'N[C@@H](CCSC)C(=O)O',
    'Ser': 'N[C@@H](CO)C(=O)O',
    'Thr': 'N[C@@H]([C@H](C)O)C(=O)O',
    'Cys': 'N[C@@H](CS)C(=O)O',
    'Tyr': 'O=C([C@RH](CC:C:C:C(:C:C&6:)O)N)O',
    'Asn': 'N[C@@H](CC(=O)N)C(=O)O',
    'Gln': 'N[C@@H](CCC(=O)N)C(=O)O',
    'Asp': 'N[C@@H](CC(=O)O)C(=O)O',
    'Glu': 'N[C@@H](CCC(=O)O)C(=O)O',
    'Lys': 'N[C@@H](CCCCN)C(=O)O',
    'Arg': 'N[C@@H](CCCNC(=N)N)C(=O)O',
    'His': 'N:C:[NH]:C:C&5:C[C@RH](N)C(O)=O',
}

# Post-translational modifications (native SCRIPT V2).
# Aromatic PTMs (pY, nitY, Pyl) converted from SMILES via RDKit bridge.
PTM_SCRIPT: Dict[str, str] = {
    'pS':  'N[C@@H](COP(=O)(O)O)C(=O)O',
    'pT':  'N[C@@H]([C@H](C)OP(=O)(O)O)C(=O)O',
    'pY':  'N[C@RH](C(O)=O)CC:C:C:C(:C:C&6:)O[P](=O)(O)O',
    'mK':  'N[C@@H](CCCCNC)C(=O)O',
    'mR':  'N[C@@H](CCCNC(=N)NC)C(=O)O',
    'acK': 'N[C@@H](CCCCNC(=O)C)C(=O)O',
    'ubK': 'N[C@@H](CCCCN)C(=O)O',
    'oxM': 'N[C@@H](CCS(=O)C)C(=O)O',
    'Hyp': 'N[C@@H](CO)C(=O)O',
    'Hyl': 'N[C@@H](CCCO)C(=O)O',
    'Orn': 'N[C@@H](CCCN)C(=O)O',
    'Cit': 'N[C@@H](CCCNC(=O)N)C(=O)O',
    'Sec': 'N[C@@H](C[Se])C(=O)O',
    'Pyl': 'C(CCNC(=O)CCC=NC&5-)C[C@RH](N)C(=O)O',
    'nitY': 'C([C@RH](N)C(=O)O)C:C:C:C(:C:C&6:)[N+](=O)[O-]',
    'suK': 'N[C@@H](CCCCNS(=O)(=O)O)C(=O)O',
    'Dpr': 'N[C@@H](CN)C(=O)O',
    'Dab': 'N[C@@H](CCN)C(=O)O',
    'Sar': 'CNCC(=O)O',
    'Gla': 'N[C@@H](CCC(=O)O)C(=O)O',
}

# Nucleotide base SCRIPT V2 notation.
# All nucleotide bases have been converted from SMILES (c1... aromatic
# notation) to native SCRIPT V2 (`:` aromatic bonds, `&N:` ring closures).
# This ensures {dA}, {rG}, {m5C} etc. parse correctly without SMILES.
NUCLEOTIDE_SCRIPT: Dict[str, str] = {
    'A':  'N:C:N:C(:[C]:[C]&6::[NH]:C:N&5:)N',
    'G':  '[C]:N:C:[NH]:[C]&5::N:C(N):[NH]:[C]=O',
    'C':  'N:[C](=O):[NH]:C:C:C&7:N',
    'T':  'C:[NH]:[C](:[NH]:[C](=O):C&7:C)=O',
    'U':  '[NH]:[C]:[O]:C:C:[NH]:[C]&7:=O',
    'I':  '[NH]:C:N:[C]:[NH]:C:N:[C]&5::[C]&9:=O',
    'dA': 'N:C:N:C(:[C]:[C]&6::[NH]:C:N&5:)N',
    'dG': '[C]:N:C:[NH]:[C]&5::N:C(N):[NH]:[C]=O',
    'dC': 'N:[C](=O):[NH]:C:C:C&7:N',
    'dT': 'C:[NH]:[C](:[NH]:[C](=O):C&7:C)=O',
    'rA': 'N:C:N:C(:[C]:[C]&6::[NH]:C:N&5:)N',
    'rG': '[C]:N:C:[NH]:[C]&5::N:C(N):[NH]:[C]=O',
    'rC': 'N:[C](=O):[NH]:C:C:C&7:N',
    'rU': '[NH]:[C]:[O]:C:C:[NH]:[C]&7:=O',
}

# Nucleotide modification SCRIPT V2 notation.
# All epigenetic modification codes converted from SMILES to V2.
NUC_MOD_SCRIPT: Dict[str, str] = {
    'm5C':  'N:[C](:[NH]:C:C(C):C&7:N)=O',
    'm6A':  '[C]:[C](:C(:N:C:N&6:)NC):N:C:[NH]&11:',
    'hm5C': 'N:[C](=O):[NH]:C:C:C&7:N',
    'f5C':  'C:C(:C(N):N:[C](:[NH]&7:)=O)C=O',
    'ca5C': 'C(=O)(O)C:C:[NH]:[C](=O):N:C&7:N',
    'm1A':  'C[N]:[C]:N:C:N:C(:[C]&6::N:C&9:)N',
    'm1G':  'N:C:[N](C):[C]:[C]&6::[C](=O):[NH]:C(N):N&8:',
    'm2G':  'C:N:[C]:[C](:N:C(:[NH]:[C]&6:=O)NC):[NH]&12:',
    'm22G': '[C](:[C]:[C](:N:C(:[NH]&6:)N(C)C):[NH]:C:N&11:)=O',
    'm7G':  '[NH]:[C]:[C](:N:C(:[C]&6:=O)N)=NN(CN=&11-)C',
    'psU':  'C(=O)NC(=O)N=CC&8-',
    's4U':  'C:C:[NH]:[C](:[NH]:[C]&6:=O)=S',
    'diU':  'NCCC(=O)NC&7-=O',
}

# Combined lookup: all monomer codes.
# Vikarana (expansion) lookup — amino acids take precedence over nucleotides
# for ambiguous single-letter codes. In standard biochemistry, bare `A` in a
# peptide context means alanine; `dA`/`rA` (unambiguous) mean adenine.
# Previously, NUCLEOTIDE_SCRIPT overwrote AMINO_ACID_SCRIPT for 'A', 'G', 'C',
# etc., causing {A} to expand as adenine (which also used unparseable SMILES
# aromatic notation). Now AMINO_ACID_SCRIPT is merged last so it wins.
ALL_MONOMERS: Dict[str, str] = {
    **NUCLEOTIDE_SCRIPT,
    **NUC_MOD_SCRIPT,
    **PTM_SCRIPT,
    **AMINO_ACID_SCRIPT,  # takes precedence for ambiguous codes (A, G, C, ...)
}


class PeptideHandler:
    """
    Handle biopolymer macro notation: {A.G.S} peptides and {dA.dG.dC.dT} nucleic acids.
    RDKit-free: uses native SCRIPTParser to expand monomer SCRIPT strings.
    """

    def __init__(self, state=None):
        self.state = state
        self.monomers = ALL_MONOMERS
        # Lazy-import parser to avoid circular dependency
        from .parser import SCRIPTParser
        self._parser = SCRIPTParser()

    # Public interface

    def handle(self, tree: Tree) -> None:
        """
        Walk a peptide_chain parse tree and populate self.state with atoms/bonds.
        Each monomer is expanded by parsing its native SCRIPT string and merging
        the resulting CoreMolecule into the state machine.
        """
        if self.state is None:
            return

        monomers = self._collect_monomers(tree)
        for code in monomers:
            script_str = self.monomers.get(code)
            if script_str:
                self._add_monomer_script(script_str)
            # Unknown codes are silently skipped (forward-compatible)

    # Tree walking

    def _collect_monomers(self, tree: Tree) -> List[str]:
        """Extract ordered list of monomer codes from a peptide_chain tree."""
        codes = []
        for subtree in tree.iter_subtrees():
            data = subtree.data.lstrip('!')
            if data == 'monomer':
                for token in subtree.children:
                    if isinstance(token, Token):
                        codes.append(str(token))
                        break
        return codes

    def _add_monomer_script(self, script_str: str) -> None:
        """Parse a monomer SCRIPT and transfer atoms/bonds into self.state."""
        try:
            result = self._parser.parse(script_str)
            if not result["success"] or result["molecule"] is None:
                return

            core = result["molecule"]
            if hasattr(core, 'atoms') and hasattr(core, 'bonds'):
                offset = len(self.state.mol.atoms)
                for atom in core.atoms:
                    self.state.mol.add_atom(atom)
                for bond in core.bonds:
                    from .mol import CoreBond
                    new_bond = CoreBond(
                        bond.begin_atom_idx + offset,
                        bond.end_atom_idx + offset,
                        bond.bond_type,
                        bond_dir=bond.bond_dir,
                        hapticity=bond.hapticity,
                        bond_class=bond.bond_class,
                    )
                    self.state.mol.add_bond_obj(new_bond)
        except Exception:
            pass  # Degrade gracefully - unknown monomers skip

    # Convenience methods

    def expand_to_script(self, sequence: List[str]) -> Optional[str]:
        """Return dot-joined SCRIPT for a sequence of monomer codes."""
        parts = []
        for code in sequence:
            s = self.monomers.get(code)
            if s:
                parts.append(s)
        return '.'.join(parts) if parts else None

    def is_known(self, code: str) -> bool:
        return code in self.monomers

    def monomer_type(self, code: str) -> str:
        """Return 'amino_acid', 'ptm', 'nucleotide', 'nuc_mod', or 'unknown'."""
        if code in AMINO_ACID_SCRIPT:
            return 'amino_acid'
        if code in PTM_SCRIPT:
            return 'ptm'
        if code in NUCLEOTIDE_SCRIPT:
            return 'nucleotide'
        if code in NUC_MOD_SCRIPT:
            return 'nuc_mod'
        return 'unknown'


# Convenience functions

def expand_peptide_to_script(sequence_str: str) -> Optional[str]:
    """Expand a dot-separated monomer sequence string to SCRIPT."""
    codes = [c.strip() for c in sequence_str.strip('{}').split('.') if c.strip()]
    h = PeptideHandler()
    return h.expand_to_script(codes)


def is_valid_monomer(code: str) -> bool:
    return code in ALL_MONOMERS