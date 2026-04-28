"""
SCRIPT Biopolymer Handler
Handles peptide {A.G.S} and nucleic acid {A.G.C.T} macro notation.
Translates monomer codes into CoreMolecule atoms via the GenerativeStateMachine.
"""

from typing import Dict, List, Optional, Tuple, Any
from lark import Tree, Token


# ──────────────────────────────────────────────────────────────
# Amino acid residue SCRIPT notation (internal representations)
# ──────────────────────────────────────────────────────────────
AMINO_ACID_SCRIPT: Dict[str, str] = {
    'A': 'N[C@@H](C)C(=O)O',
    'R': 'N[C@@H](CCCNC(=N)N)C(=O)O',
    'N': 'N[C@@H](CC(=O)N)C(=O)O',
    'D': 'N[C@@H](CC(=O)O)C(=O)O',
    'C': 'N[C@@H](CS)C(=O)O',
    'E': 'N[C@@H](CCC(=O)O)C(=O)O',
    'Q': 'N[C@@H](CCC(=O)N)C(=O)O',
    'G': 'NCC(=O)O',
    'H': 'N[C@@H](CC1=CNC=N1)C(=O)O',
    'I': 'N[C@@H]([C@@H](C)CC)C(=O)O',
    'L': 'N[C@@H](CC(C)C)C(=O)O',
    'K': 'N[C@@H](CCCCN)C(=O)O',
    'M': 'N[C@@H](CCSC)C(=O)O',
    'F': 'N[C@@H](Cc1ccccc1)C(=O)O',
    'P': 'N1[C@@H](CCC1)C(=O)O',
    'S': 'N[C@@H](CO)C(=O)O',
    'T': 'N[C@@H]([C@H](C)O)C(=O)O',
    'W': 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
    'Y': 'N[C@@H](Cc1ccc(O)cc1)C(=O)O',
    'V': 'N[C@@H](C(C)C)C(=O)O',
    # 3-letter codes (alias)
    'Ala':'N[C@@H](C)C(=O)O',   'Gly':'NCC(=O)O',
    'Val':'N[C@@H](C(C)C)C(=O)O', 'Leu':'N[C@@H](CC(C)C)C(=O)O',
    'Ile':'N[C@@H]([C@@H](C)CC)C(=O)O', 'Pro':'N1[C@@H](CCC1)C(=O)O',
    'Phe':'N[C@@H](Cc1ccccc1)C(=O)O', 'Trp':'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
    'Met':'N[C@@H](CCSC)C(=O)O', 'Ser':'N[C@@H](CO)C(=O)O',
    'Thr':'N[C@@H]([C@H](C)O)C(=O)O', 'Cys':'N[C@@H](CS)C(=O)O',
    'Tyr':'N[C@@H](Cc1ccc(O)cc1)C(=O)O', 'Asn':'N[C@@H](CC(=O)N)C(=O)O',
    'Gln':'N[C@@H](CCC(=O)N)C(=O)O', 'Asp':'N[C@@H](CC(=O)O)C(=O)O',
    'Glu':'N[C@@H](CCC(=O)O)C(=O)O', 'Lys':'N[C@@H](CCCCN)C(=O)O',
    'Arg':'N[C@@H](CCCNC(=N)N)C(=O)O', 'His':'N[C@@H](CC1=CNC=N1)C(=O)O',
}

# Post-translational modifications
PTM_SCRIPT: Dict[str, str] = {
    'pS':  'N[C@@H](COP(=O)(O)O)C(=O)O',
    'pT':  'N[C@@H]([C@H](C)OP(=O)(O)O)C(=O)O',
    'pY':  'N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O',
    'mK':  'N[C@@H](CCCCNC)C(=O)O',
    'mR':  'N[C@@H](CCCNC(=N)NC)C(=O)O',
    'acK': 'N[C@@H](CCCCNC(=O)C)C(=O)O',
    'ubK': 'N[C@@H](CCCCN)C(=O)O',   # simplified
    'oxM': 'N[C@@H](CCS(=O)C)C(=O)O',
    'Hyp': 'N[C@@H](CO)C(=O)O',       # hydroxyproline (simplified)
    'Hyl': 'N[C@@H](CCCO)C(=O)O',     # hydroxylysine (simplified)
    'Orn': 'N[C@@H](CCCN)C(=O)O',
    'Cit': 'N[C@@H](CCCNC(=O)N)C(=O)O',
    'Sec': 'N[C@@H](C[Se])C(=O)O',
    'Pyl': 'N[C@@H](CCCCNC(=O)C1CC=NC1)C(=O)O',
    'nitY':'N[C@@H](Cc1ccc([N+](=O)[O-])cc1)C(=O)O',
    'suK': 'N[C@@H](CCCCNS(=O)(=O)O)C(=O)O',
    'Dpr': 'N[C@@H](CN)C(=O)O',
    'Dab': 'N[C@@H](CCN)C(=O)O',
    'Sar': 'CNCC(=O)O',
    'Gla': 'N[C@@H](CCC(=O)O)C(=O)O',  # gamma-carboxyglutamic (simplified)
}

# Nucleotide base SCRIPT notation
NUCLEOTIDE_SCRIPT: Dict[str, str] = {
    'A':  'Nc1ncnc2[nH]cnc12',               # adenine
    'G':  'Nc1nc2[nH]cnc2c(=O)[nH]1',        # guanine
    'C':  'Nc1cc[nH]c(=O)n1',               # cytosine — use non-aromatic tautomer
    'T':  'Cc1c[nH]c(=O)[nH]c1=O',          # thymine
    'U':  'O=C1C=CN([H])C(=O)N1',           # uracil
    'I':  'O=c1[nH]cnc2[nH]cnc12',          # hypoxanthine / inosine base
    'dA': 'Nc1ncnc2[nH]cnc12',
    'dG': 'Nc1nc2[nH]cnc2c(=O)[nH]1',
    'dC': 'Nc1ccn([H])c(=O)n1',
    'dT': 'Cc1c[nH]c(=O)[nH]c1=O',
    'rA': 'Nc1ncnc2[nH]cnc12',
    'rG': 'Nc1nc2[nH]cnc2c(=O)[nH]1',
    'rC': 'Nc1ccn([H])c(=O)n1',
    'rU': 'O=C1C=CN([H])C(=O)N1',
}

# Nucleotide modification SCRIPT notation
NUC_MOD_SCRIPT: Dict[str, str] = {
    'm5C':  'Cc1cn([H])c(=O)nc1N',           # 5-methylcytosine
    'm6A':  'CNC1=NC=NC2=C1N=CN2',           # N6-methyladenine
    'hm5C': 'NCO.Nc1cc[nH]c(=O)n1',         # 5-hydroxymethylcytosine (simplified)
    'f5C':  'O=Cc1cn([H])c(=O)nc1N',        # 5-formylcytosine
    'ca5C': 'OC(=O)c1cn([H])c(=O)nc1N',     # 5-carboxylcytosine
    'm1A':  'Cn1cnc2c1ncnc2N',              # 1-methyladenine
    'm1G':  'CN1C=NC2=C1N=C(N)NC2=O',       # 1-methylguanosine
    'm2G':  'CNC1=NC2=C(N=CN2)C(=O)N1',     # N2-methylguanosine (simplified)
    'm22G': 'CN(C)c1nc2[nH]cnc2c(=O)[nH]1', # N2,N2-dimethylguanosine (simplified)
    'm7G':  'CN1CN=C2NC(=O)C(=NC2=N1)N',    # 7-methylguanosine
    'psU':  'O=C1NC(=O)CC=N1',              # pseudouridine
    's4U':  'O=C1C=CN([H])C(=S)N1',         # 4-thiouridine
    'diU':  'O=C1NC(=O)CCN1',              # dihydrouridine
    'I':    'O=c1[nH]cnc2[nH]cnc12',        # inosine/hypoxanthine
}


# Combined lookup: all monomer codes
ALL_MONOMERS: Dict[str, str] = {
    **AMINO_ACID_SCRIPT,
    **PTM_SCRIPT,
    **NUCLEOTIDE_SCRIPT,
    **NUC_MOD_SCRIPT,
}


class PeptideHandler:
    """
    Handle biopolymer macro notation: {A.G.S} peptides and {A.G.C.T} nucleic acids.

    Instantiated with a GenerativeStateMachine; the handle(tree) method
    walks the peptide_chain parse tree and adds atoms/bonds to the state machine
    by parsing each monomer's SCRIPT through RDKit and transferring atoms.
    """

    def __init__(self, state=None):
        self.state = state
        self.monomers = ALL_MONOMERS

    # ── Public interface ────────────────────────────────────────

    def handle(self, tree: Tree) -> None:
        """
        Walk a peptide_chain parse tree and populate self.state with atoms/bonds.
        Each monomer is expanded to its atomic SCRIPT and added to the molecule.
        The current implementation stores monomers as disconnected fragments
        (full polymer connectivity requires a proper peptide bond builder).
        """
        if self.state is None:
            return

        monomers = self._collect_monomers(tree)
        for code in monomers:
            script_str = self.monomers.get(code)
            if script_str:
                self._add_monomer_script(script_str)
            # Unknown codes are silently skipped (forward-compatible)

    # ── Tree walking ────────────────────────────────────────────

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
            # Also catch bare nucleotide/PTM tokens outside monomer tree
            # (grammar may route them differently under monomer_element)
            elif data == 'monomer_element':
                pass  # handled via monomer subtree above
        return codes

    def _add_monomer_script(self, script_str: str) -> None:
        """Parse a monomer SCRIPT and transfer atoms/bonds into self.state."""
        try:
            from rdkit import Chem
            from script.rdkit_bridge import from_rdkit
            mol = Chem.MolFromSmiles(script_str)
            if mol is None:
                return
            core = from_rdkit(mol)
            # Offset indices by current atom count in state
            offset = len(self.state.mol.atoms)
            for atom in core.atoms:
                self.state.mol.add_atom(atom)
            for bond in core.bonds:
                from script.mol import CoreBond
                new_bond = CoreBond(
                    bond.begin_atom_idx + offset,
                    bond.end_atom_idx + offset,
                    bond.bond_type,
                )
                self.state.mol.add_bond_obj(new_bond)
        except Exception:
            pass  # Degrade gracefully — unknown monomers skip

    # ── Convenience methods ─────────────────────────────────────

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
        if code in AMINO_ACID_SCRIPT: return 'amino_acid'
        if code in PTM_SCRIPT:        return 'ptm'
        if code in NUCLEOTIDE_SCRIPT: return 'nucleotide'
        if code in NUC_MOD_SCRIPT:    return 'nuc_mod'
        return 'unknown'


# ── Convenience functions ────────────────────────────────────────

def expand_peptide_to_script(sequence_str: str) -> Optional[str]:
    """Expand a dot-separated monomer sequence string to SCRIPT."""
    codes = [c.strip() for c in sequence_str.strip('{}').split('.') if c.strip()]
    h = PeptideHandler()
    return h.expand_to_script(codes)

def is_valid_monomer(code: str) -> bool:
    return code in ALL_MONOMERS