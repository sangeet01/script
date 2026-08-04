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
    SPLINE     = 10  # ~>  spline beam (V4 L1 — LEAP71 RandomSplineLattice)
    BRIDGE     = 11  # <>  3-center-2-electron bridge (V4.3 — diborane, carboranes)


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
    PLANAR           = 7   # @PL  planar chirality (metallocenes, ferrocenes)


ATOMIC_NUM_TO_SYMBOL = {
    1: 'H',   2: 'He',  3: 'Li',  4: 'Be',  5: 'B',   6: 'C',   7: 'N',   8: 'O',  9: 'F',  10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',  16: 'S',  17: 'Cl', 18: 'Ar', 19: 'K',  20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V',  24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',  40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I',  54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W',  75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U',  93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',100: 'Fm',
   101: 'Md',102: 'No',103: 'Lr',104: 'Rf',105: 'Db',106: 'Sg',107: 'Bh',108: 'Hs',109: 'Mt',110: 'Ds',
   111: 'Rg',112: 'Cn',113: 'Nh',114: 'Fl',115: 'Mc',116: 'Lv',117: 'Ts',118: 'Og',
}


class CoreAtom:
    def __init__(self, atomic_num: int, formal_charge: int = 0, 
                 isotope: int = 0, radical_electrons: int = 0,
                 symbol: str = "", is_aromatic: bool = False,
                 mapping: int = 0, occupancy: float = 1.0,
                 spin: int = 0, is_excited: bool = False,
                 is_wildcard: bool = False,
                 beam_radius: Optional[float] = None):
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
        # V4 L2: beam radius ratio (LEAP71 IBeamThickness).  None = no
        # explicit per-vertex radius; the bridge falls back to the
        # Thickness:class on the molecule.  When set, overrides the
        # class for beams touching this vertex.
        self.beam_radius: Optional[float] = beam_radius
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
                 translation: Tuple[int, int, int] = (0, 0, 0),
                 control_points: Optional[List[Tuple[float, float, float]]] = None):
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
        self.bond_class = bond_class  # "dative","rev_dative","coordinate","star","spline",""
        self.is_rc = False            # ring closure bond
        self.is_aromatic = False      # part of aromatic/resonant system
        # Periodic topology: integer lattice translation vector (tx, ty, tz).
        # (0,0,0) for all intracell bonds; non-zero for bonds that cross unit-cell
        # boundaries in MOF/zeolite frameworks.  Ignored for non-periodic molecules.
        # [V4.2 Q2] Now accepts floats for crystallographic fractional translations.
        self.translation = translation
        # [V4.2 Q5] Explicit spline control points (None for non-spline bonds).
        # Each point is a (x, y, z) tuple. Used by RandomSplineLattice and
        # hand-authored curves (protein loops, RNA pseudoknots).
        self.control_points: Optional[List[Tuple[float, float, float]]] = control_points

class PolymerBlock:
    """
    One segment in a block copolymer.

    A CoreMolecule that represents a block copolymer (e.g. ABA triblock)
    stores a list of PolymerBlock objects in its polymer_blocks attribute.
    Each block has:
      - unit:         the CoreMolecule for the repeat unit of this block
                       (kept for informational purposes; may be None when
                       the block was expanded directly into the shared
                       parent graph — see atom_start/atom_end below)
      - repeat_count: int, (min,max) tuple, or 'n' (symbolic)
      - block_kind:   junction token from the grammar
                      ('diblock', 'alternating', 'random', 'graft', or '')
      - atom_start:   index of this block's first atom in the PARENT
                      molecule's shared atom list (V3.7 full expansion)
      - atom_end:     index of this block's last atom (inclusive) in the
                      parent molecule's shared atom list
    """
    def __init__(self, unit: Optional['CoreMolecule'],
                 repeat_count: Any = None,
                 block_kind: str = '',
                 atom_start: Optional[int] = None,
                 atom_end: Optional[int] = None):
        self.unit: Optional['CoreMolecule'] = unit
        self.repeat_count: Any = repeat_count
        self.block_kind: str = block_kind
        self.atom_start: Optional[int] = atom_start
        self.atom_end: Optional[int] = atom_end

    def __repr__(self) -> str:
        atom_count = len(getattr(self.unit, "atoms", []) or [])
        return ("PolymerBlock(kind=" + repr(self.block_kind) +
                ", repeat=" + repr(self.repeat_count) +
                ", atoms=" + str(atom_count) + ")")


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

        # Graph version counter — bumped on every add_atom / add_bond /
        # add_bond_obj call. Used as a cache key by ranking.calculate_ranks
        # so that repeated calls (e.g. once per chiral atom in
        # ChiralResolver.resolve) do not re-run the full Morgan-WL loop.
        # Initialised to -1 so the first call always misses.
        self._graph_version: int = -1
        self._rank_cache: Optional[Dict[int, int]] = None

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

        # ---- V4 Lattice Extension (LEAP71 Bridge) ----
        # lattice_type: tagged on molecules parsed from
        #   [Lattice:BodyCentered] {...} etc.  Tells the bridge which
        #   ILatticeType subclass to construct.
        self.lattice_type: Optional[str] = None
        # thickness_class + thickness_args: tagged from
        #   [Thickness:Constant(2.0)] etc.  Tells the bridge which
        #   IBeamThickness subclass to construct.
        self.thickness_class: Optional[str] = None
        self.thickness_args: Optional[Tuple] = None
        # post_process_ops: ordered list of (op_name, args) tuples
        #   applied to the voxelized lattice.  Render-time only —
        #   does not affect canonicalization of the graph.
        #   Example: [("overoffset", (3.0, 0.0)), ("intersect", ())]
        self.post_process_ops: List[Tuple[str, Tuple]] = []
        # namespace: optional disambiguator parsed from [[geom:..]] /
        #   [[xtal:..]].  "geom" = bounding shape (LEAP71 BaseShape),
        #   "xtal" = crystallographic context (existing V3 semantics),
        #   None = bare label (auto-dispatch by Base* prefix).
        self.context_namespace: Optional[str] = None

        # [V4.1] Generalized typed tags: list of (namespace, value, args)
        # tuples for any [Namespace:Value(args)] tags that weren't
        # recognized as Lattice or Thickness.  Order preserved.
        # Example: [("Mesh", "Icosphere", (2,)), ("Material", "Steel", ())]
        self.typed_tags: List[Tuple[str, str, Tuple]] = []

    def add_atom(self, atom: CoreAtom) -> int:
        idx = len(self.atoms)
        self.atoms.append(atom)
        self.adj[idx] = []
        self._graph_version += 1
        return idx

    def add_bond(self, begin_idx: int, end_idx: int, bond_type: Any,
                 bond_dir: int = 0, hapticity: int = 0, bond_class: str = "",
                 translation: Tuple[int, int, int] = (0, 0, 0),
                 control_points: Optional[List[Tuple[float, float, float]]] = None):
        bond_idx = len(self.bonds)
        bond = CoreBond(begin_idx, end_idx, bond_type, bond_dir,
                        hapticity=hapticity, bond_class=bond_class,
                        translation=translation,
                        control_points=control_points)
        self.bonds.append(bond)
        self.adj[begin_idx].append((end_idx, bond_idx))
        self.adj[end_idx].append((begin_idx, bond_idx))
        if translation != (0, 0, 0):
            self.is_periodic = True
        self._graph_version += 1

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
        self._graph_version += 1

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


# ============================================================================
# Polyatomic ion registry
# ============================================================================
# Maps FORMULA (canonical uppercase, no charge) -> {charge_str: expansion_script}
#
# The expansion is a valid SCRIPT substring that produces the correct atomic
# graph when parsed.  The parser's _handle_bracket_atom method looks up this
# registry when it encounters a POLYATOMIC_FORMULA token in the grammar.
#
# This follows the same pattern as amino acid residues in peptide.py: the
# expansion data lives as reference data in mol.py, and the parser handles
# expansion through the normal state-machine path (Sandhi-validated).
#
# Charge strings use the canonical SCRIPT charge suffix:
#   '-'   for -1,   '--'  for -2,   '-3'  for -3
#   '+'   for +1,   '++'  for +2,   '+3'  for +3

POLYATOMIC_IONS: Dict[str, Dict[str, str]] = {
    # ---- Sulfur oxoanions ----
    'SO4': {
        '-2': '[S]([O-])([O-])(=O)(=O)',         # sulfate
        '-':  '[S]([O-])(=O)(=O)([OH])',         # bisulfate (HSO4-)
    },
    'SO3': {
        '-2': '[S]([O-])([O-])(=O)',             # sulfite
        '-':  '[S]([O-])(=O)([OH])',             # bisulfite (HSO3-)
    },
    'S2O3': {
        '-2': '[S]([O-])([O-])(=O)[S]',          # thiosulfate
    },
    'S2O8': {
        '-2': '[S]([O-])(=O)(=O)[O][O][S]([O-])(=O)(=O)',  # peroxodisulfate
    },
    # ---- Nitrogen oxoanions ----
    'NO3': {
        '-': '[N+]([O-])(=O)=O',                  # nitrate
    },
    'NO2': {
        '-': '[N+]([O-])=O',                      # nitrite
    },
    'N3': {
        '-': '[N-]=[N+]=[N-]',                    # azide
    },
    # ---- Carbon oxoanions ----
    'CO3': {
        '-2': '[C]([O-])([O-])(=O)',             # carbonate
        '-':  '[C](=O)([O-])([OH])',             # bicarbonate (HCO3-)
    },
    'C2O4': {
        '-2': '[C](=O)([O-])[C](=O)([O-])',      # oxalate
    },
    'CN': {
        '-': '[C]#[N-]',                         # cyanide
    },
    'OCN': {
        '-': '[O][C]#[N-]',                      # cyanate
    },
    'SCN': {
        '-': '[S][C]#[N-]',                      # thiocyanate
    },
    # ---- Phosphorus oxoanions ----
    'PO4': {
        '-3': '[P]([O-])([O-])([O-])(=O)',       # phosphate
        '-2': '[P](=O)([O-])([O-])([OH])',       # monohydrogen phosphate (HPO4-2)
        '-':  '[P](=O)([O-])([OH])([OH])',       # dihydrogen phosphate (H2PO4-)
    },
    'P2O7': {
        '-4': '[P](=O)([O-])([O-])[O][P](=O)([O-])([O-])',  # pyrophosphate
    },
    # ---- Boron ----
    'BO3': {
        '-3': '[B]([O-])([O-])([O-])',           # borate
    },
    'B4O7': {
        '-2': '[B](=O)([O-])[B]([O-])[B]([O-])[B](=O)([O-])',  # tetraborate (simplified)
    },
    # ---- Halogen oxoanions ----
    'ClO': {
        '-': '[Cl][O-]',                         # hypochlorite
    },
    'ClO2': {
        '-': '[Cl]([O-])=O',                     # chlorite
    },
    'ClO3': {
        '-': '[Cl]([O-])(=O)=O',                 # chlorate
    },
    'ClO4': {
        '-': '[Cl]([O-])(=O)(=O)=O',             # perchlorate
    },
    'BrO3': {
        '-': '[Br]([O-])(=O)=O',                 # bromate
    },
    'BrO4': {
        '-': '[Br]([O-])(=O)(=O)=O',             # perbromate
    },
    'IO3': {
        '-': '[I]([O-])(=O)=O',                  # iodate
    },
    'IO4': {
        '-': '[I]([O-])(=O)(=O)=O',              # periodate
    },
    # ---- Manganese / Chrome / permanganate ----
    'MnO4': {
        '-':  '[Mn]([O-])(=O)(=O)=O',            # permanganate (Mn +VII)
        '-2': '[Mn]([O-])(=O)(=O)([O-])',        # manganate (Mn +VI)
    },
    'CrO4': {
        '-2': '[Cr]([O-])([O-])(=O)=O',          # chromate
    },
    'Cr2O7': {
        '-2': '[Cr](=O)([O-])[O][Cr](=O)(=O)([O-])',  # dichromate
    },
    # ---- Organic acid anions (common names) ----
    'CH3COO': {
        '-': 'CC(=O)[O-]',                       # acetate
    },
    'C2H3O2': {
        '-': 'CC(=O)[O-]',                       # acetate (IUPAC formula)
    },
    'AcO': {
        '-': 'CC(=O)[O-]',                       # acetate (shorthand AcO-)
    },
    'HCOO': {
        '-': '[C](=O)[O-]',                      # formate (HCOO-)
    },
    # ---- Silicate / aluminate ----
    'SiO3': {
        '-2': '[Si]([O-])([O-])(=O)',            # metasilicate
    },
    'AlO2': {
        '-': '[Al]([O-])=O',                     # aluminate (meta)
    },
    # ---- Peroxides / superoxides ----
    'O2': {
        '-':  '[O-][O]',                         # superoxide (O2-.)
        '-2': '[O-][O-]',                        # peroxide (O2-2)
    },
}