"""
Core SCRIPT Graph Representation
Lightweight, RDKit-independent molecule data structure.
"""

from __future__ import annotations
from typing import List, Dict, Tuple, Optional, Any
from enum import IntEnum


class BondType(IntEnum):
    """
    Typed bond order enum for CoreBond.
    Replaces bare integers — makes the IR self-documenting and safe.
    """
    SINGLE     = 1   # —   standard covalent single bond
    DOUBLE     = 2   # =   double bond
    TRIPLE     = 3   # #   triple bond
    AROMATIC   = 4   # :   resonant / aromatic (Kekule-free)
    DATIVE     = 5   # ->  donor->acceptor (Lewis acid/base, BN adducts)
    REV_DATIVE = 6   # <-  reverse dative
    TAUTOMERIC = 7   # =:  mobile / tautomeric (keto-enol etc.)
    COORDINATE = 8   # >   coordinate / haptic (organometallics)
    STAR       = 9   # *   resonance / haptic wildcard


class StereoType(IntEnum):
    """
    Extended stereochemistry type for non-tetrahedral centres.
    Tetrahedral R/S is handled by CoreAtom._initial_tag (existing).
    """
    NONE             = 0
    TETRAHEDRAL      = 1   # @ / @@   (existing, kept for completeness)
    SQUARE_PLANAR    = 2   # @SP
    OCTAHEDRAL       = 3   # @OH
    ATROPISOMER      = 4   # @AX  axial chirality (biaryls, allenes)
    TRIG_BIPYRAMIDAL = 5   # @TB
    PYRAMIDAL        = 6   # @PY


class CoreAtom:
    def __init__(self, atomic_num: int, formal_charge: int = 0, 
                 isotope: int = 0, radical_electrons: int = 0,
                 symbol: str = "", is_aromatic: bool = False,
                 mapping: int = 0, occupancy: float = 1.0,
                 spin: int = 0, is_excited: bool = False,
                 is_wildcard: bool = False):
        self.atomic_num = atomic_num
        self.formal_charge = formal_charge
        self.isotope = isotope
        self.radical_electrons = radical_electrons
        self.symbol = symbol
        self.is_aromatic = is_aromatic
        self.mapping = mapping
        self.occupancy = occupancy
        self.spin = spin
        self.is_excited = is_excited
        self.is_wildcard = is_wildcard       # True for '*' query/scaffold atoms
        self.rank = -1
        self.coords: Optional[Tuple[float, float, float]] = None
        self.implicit_hs: int = 0

        self._initial_tag = 0
        self._initial_nbrs: List[int] = []
        self.chirality = 0  # 0: None, 1: CW (@@), 2: CCW (@)

        # Extended stereochemistry (non-tetrahedral)
        self.stereo_type: StereoType = StereoType.NONE
        # Ordered neighbor indices for coordination stereo (filled by bridge)
        self.stereo_neighbors: List[int] = []

        # Query atom fields (Praśna — pattern atoms, not canonical representation)
        self.is_query: bool = False          # True if this is a query/pattern atom
        self.query_atomic_nums: List[int] = []   # [#6] -> [6], [#6,#7] -> [6,7]
        self.query_not: bool = False         # [!N] -> negation
        self.query_ring: Optional[int] = None    # [R] -> 0 (any ring), [r5] -> 5
        self.query_valence: Optional[int] = None # [v3] -> 3
        self.query_hcount: Optional[int] = None  # [H2] -> 2
        self.query_aromatic: Optional[bool] = None  # [a] -> True, [A] -> False
        self.query_primitives: List[dict] = []   # raw primitives for complex queries

class CoreBond:
    def __init__(self, begin_atom_idx: int, end_atom_idx: int,
                 bond_type: Any, bond_dir: int = 0,
                 hapticity: int = 0, bond_class: str = "",
                 translation: Tuple[int, int, int] = (0, 0, 0)):
        self.begin_atom_idx = begin_atom_idx
        self.end_atom_idx = end_atom_idx
        # Normalise to BondType enum if an int is passed (backward compat)
        if isinstance(bond_type, int) and not isinstance(bond_type, BondType):
            try:
                bond_type = BondType(bond_type)
            except ValueError:
                bond_type = BondType.SINGLE
        self.bond_type: BondType = bond_type
        self.bond_dir = bond_dir      # 0: None, 3: Up(/), 4: Down(\)
        self.hapticity = hapticity    # eta-n for haptic organometallics
        self.bond_class = bond_class  # "dative","rev_dative","coordinate","star",""
        self.is_rc = False            # ring closure bond
        self.is_aromatic = False      # part of aromatic/resonant system
        # Periodic topology: integer lattice translation vector (tx, ty, tz).
        # (0,0,0) for all intracell bonds; non-zero for bonds that cross unit-cell
        # boundaries in MOF/zeolite frameworks.  Ignored for non-periodic molecules.
        self.translation: Tuple[int, int, int] = translation

class PolymerBlock:
    """
    One segment in a block copolymer.

    A CoreMolecule that represents a block copolymer (e.g. ABA triblock)
    stores a list of PolymerBlock objects in its polymer_blocks attribute.
    Each block has:
      - unit:         the CoreMolecule for the repeat unit of this block
      - repeat_count: int, (min,max) tuple, or 'n' (symbolic)
      - block_kind:   junction token from the grammar
                      ('diblock', 'alternating', 'random', 'graft', or '')
    """
    def __init__(self, unit: 'CoreMolecule',
                 repeat_count: Any = None,
                 block_kind: str = ''):
        self.unit: 'CoreMolecule' = unit
        self.repeat_count: Any = repeat_count
        self.block_kind: str = block_kind

    def __repr__(self) -> str:
        return ("PolymerBlock(kind=" + repr(self.block_kind) +
                ", repeat=" + repr(self.repeat_count) +
                ", atoms=" + str(len(self.unit.atoms)) + ")")


class CoreMolecule:
    """
    RDKit-independent molecular graph.
    Used for canonicalization and parsing logic.
    """
    def __init__(self):
        self.atoms: List[CoreAtom] = []
        self.bonds: List[CoreBond] = []
        self.adj: Dict[int, List[Tuple[int, int]]] = {}  # idx -> list of (neighbor_idx, bond_idx)
        self.chiral_centers: Dict[int, int] = {} # atom_idx -> chirality_bit (0:CCW, 1:CW)
        self.macroscopic_context: Optional[str] = None

        # Tier 2 fields — semantic metadata above the atom/bond graph
        self.fragment_separator: str = "."      # "." solvate/salt | "~" ionic pair
        self.repeat_count: Optional[Any] = None # polymer: int, (min,max) tuple, or None
        self.phase_boundary: Optional[str] = None  # phase label when "|" is used
        # Block copolymer support: list of PolymerBlock segments
        self.polymer_blocks: List[PolymerBlock] = []
        # Junction type: 'diblock', 'triblock', 'alternating', 'random', 'graft', or None
        self.block_topology: Optional[str] = None
        # 3-D periodic topology (MOFs, zeolites, coordination polymers).
        # lattice_vectors: 3x3 row-vector matrix ((a1,a2,a3),(b1,b2,b3),(c1,c2,c3))
        #   in Angstroms.  None for non-periodic molecules.
        self.lattice_vectors: Optional[Tuple[Tuple[float,float,float], ...]] = None
        # is_periodic: True when at least one bond has a non-zero translation vector.
        self.is_periodic: bool = False

    def add_atom(self, atom: CoreAtom) -> int:
        idx = len(self.atoms)
        self.atoms.append(atom)
        self.adj[idx] = []
        return idx

    def add_bond(self, begin_idx: int, end_idx: int, bond_type: Any,
                 bond_dir: int = 0, hapticity: int = 0, bond_class: str = "",
                 translation: Tuple[int, int, int] = (0, 0, 0)):
        bond_idx = len(self.bonds)
        bond = CoreBond(begin_idx, end_idx, bond_type, bond_dir,
                        hapticity=hapticity, bond_class=bond_class,
                        translation=translation)
        self.bonds.append(bond)
        self.adj[begin_idx].append((end_idx, bond_idx))
        self.adj[end_idx].append((begin_idx, bond_idx))
        if translation != (0, 0, 0):
            self.is_periodic = True

    def add_bond_obj(self, bond: CoreBond) -> None:
        """Add a pre-constructed CoreBond object. Used by PeptideHandler."""
        bond_idx = len(self.bonds)
        self.bonds.append(bond)
        # Ensure adj lists exist for both endpoints
        if bond.begin_atom_idx not in self.adj:
            self.adj[bond.begin_atom_idx] = []
        if bond.end_atom_idx not in self.adj:
            self.adj[bond.end_atom_idx] = []
        self.adj[bond.begin_atom_idx].append((bond.end_atom_idx, bond_idx))
        self.adj[bond.end_atom_idx].append((bond.begin_atom_idx, bond_idx))

    def get_neighbors(self, atom_idx: int) -> List[int]:
        return [nbr_idx for nbr_idx, _ in self.adj.get(atom_idx, [])]

    def get_bond(self, idx1: int, idx2: int) -> Optional[CoreBond]:
        for nbr_idx, bond_idx in self.adj.get(idx1, []):
            if nbr_idx == idx2:
                return self.bonds[bond_idx]
        return None

class Reaction:
    """
    Represents a SCRIPT reaction: reactants >> products (with optional agents).

    Returned by SCRIPTParser when the input contains a '>>' arrow.
    Each side is a list of CoreMolecule objects.

    Example:
        CC>>CCO  ->  Reaction(reactants=[CC], products=[CCO], agents=[])
        CC>[Pd]>CCO  ->  Reaction(reactants=[CC], agents=[Pd], products=[CCO])
    """
    def __init__(self,
                 reactants: List[CoreMolecule],
                 products: List[CoreMolecule],
                 agents: Optional[List[CoreMolecule]] = None):
        self.reactants: List[CoreMolecule] = reactants
        self.products:  List[CoreMolecule] = products
        self.agents:    List[CoreMolecule] = agents or []

    def __repr__(self):
        r = len(self.reactants)
        a = len(self.agents)
        p = len(self.products)
        return f"Reaction(reactants={r}, agents={a}, products={p})"


# Nucleotide modification registry
# Maps NUC_MOD_CODE string to (base, modification_name)
NUCLEOTIDE_MODIFICATIONS = {
    "m5C":  ("C", "5-methylcytosine"),
    "m6A":  ("A", "N6-methyladenine"),
    "hm5C": ("C", "5-hydroxymethylcytosine"),
    "f5C":  ("C", "5-formylcytosine"),
    "ca5C": ("C", "5-carboxylcytosine"),
    "m1A":  ("A", "1-methyladenine"),
    "m1G":  ("G", "1-methylguanosine"),
    "m2G":  ("G", "N2-methylguanosine"),
    "m22G": ("G", "N2,N2-dimethylguanosine"),
    "m7G":  ("G", "7-methylguanosine"),
    "psU":  ("U", "pseudouridine"),
    "s4U":  ("U", "4-thiouridine"),
    "I":    ("A", "inosine"),
    "diU":  ("U", "dihydrouridine"),
}

# Nucleotide base registry
# Maps single/prefix codes to (full_name, is_dna, rdkit_smiles)
NUCLEOTIDE_BASES = {
    "A": ("adenine",   False, "n1cnc2ncnc(N)c12"),
    "G": ("guanine",   False, "n1cnc2c(=O)[nH]c(N)nc2n1"),
    "C": ("cytosine",  False, "n1ccc(N)nc1=O"),
    "T": ("thymine",   True,  "n1cc(C)c(=O)[nH]c1=O"),
    "U": ("uracil",    False, "n1cccc(=O)[nH]1=O"),
    "I": ("inosine",   False, "n1cnc2c(=O)[nH]cnc12"),
    "dA": ("deoxyadenosine", True,  "n1cnc2ncnc(N)c12"),
    "dG": ("deoxyguanosine", True,  "n1cnc2c(=O)[nH]c(N)nc2n1"),
    "dC": ("deoxycytidine",  True,  "n1ccc(N)nc1=O"),
    "dT": ("deoxythymidine", True,  "n1cc(C)c(=O)[nH]c1=O"),
    "rA": ("adenosine",      False, "n1cnc2ncnc(N)c12"),
    "rG": ("guanosine",      False, "n1cnc2c(=O)[nH]c(N)nc2n1"),
    "rC": ("cytidine",       False, "n1ccc(N)nc1=O"),
    "rU": ("uridine",        False, "n1cccc(=O)[nH]1=O"),
}
