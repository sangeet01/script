**SCRIPT** is a deterministic molecular notation system and RDKit-independent cheminformatics engine. By shifting from atom-based states to a grammar-level **Anubandha (State Marker)** system, SCRIPT provides a "one true string" for every molecule with a verified 99.0% round-trip fidelity.

---

## Why SCRIPT?

SMILES has served chemistry for 35 years, but its limitations are critical for modern AI/ML applications:

- **Non-canonical**: Same molecule = multiple valid SMILES strings
- **Ambiguous rings**: Global ring labels (C1...C1) create parsing complexity
- **Stereochemistry fragility**: Neighbor ordering affects chirality interpretation
- **No validation**: Invalid strings parse without error

SCRIPT addresses these systematically:

| Problem | SMILES | SCRIPT 2.0 |
|---------|--------|------------|
| Canonicalization | Multiple valid strings | Path-invariant DFS traversal |
| Ring notation | Global labels `C1...C1` | **Topological &N:** (Invariant size) |
| Aromaticity | `c1ccccc1` (lowercase) | **Anubandha :** (Grammar state) |
| Tautomers | Multiple forms | **Mobile =:** (Unified form) |
| Validation | Post-hoc | Generative state machine (Sandhi) |

---

## Core Innovations

### 1. Deterministic Canonicalization
Morgan-invariant ranking with DFS traversal ensures every molecule has exactly one canonical SCRIPT string. Tie-breaking uses lexicographical ordering of neighbor subgraphs.

```python
# Aspirin has one canonical form
SCRIPT: CC(=O)OC1=CC=CC=C1C(=O)O
```

### 2. Topological Back-counting (&N)
SCRIPT 2.0 introduces **Topological Back-counting**. The ring closure index `&6:` or `&5.` indicates the size of the cycle along the DFS path, completely ignoring branches. This ensures the ring notation remains invariant even if side-chains are modified.

```python
SMILES:  C1CCCCC1      # Global label - must track state
SCRIPT:  C1CCCCC&6.    # Topological - connect to parent 5 steps back
```

### 3. Anubandha (The State Marker System)
Aromaticity and Tautomerism are represented as **states** of the topological path using the Anubandha markers (`:` for resonant, `.` for fixed).

```python
SMILES:  c1ccccc1      # Mixed case hack
SCRIPT:  C~C~C~C~C~C&6: # Resonant path (Aromatic)
SCRIPT:  CN=C(O)C      # Keto-Enol Tautomer
SCRIPT:  CN-C(=:O)C    # Delocalized representation using =: marker
```

### 4. CIP-Based Stereochemistry
Stereochemistry is resolved through Cahn-Ingold-Prelog priority as a universal reference frame. SCRIPT's DFS neighbor order (Parent < [H] < Ring-Closures < Ring-Openings < Branches < Chain) is transformed to CIP space using permutation parity.

```python
# Chiral transformation: local_chirality XOR parity(local_order -> CIP_order)
# Handles complex cases: glucose, fused rings, symmetric centers
```

### 5. Sandhi Validation (Panini-Inspired)
Generative state machine tracks valence in real-time during parsing. Invalid structures are caught immediately, not after molecule construction.

```python
# C(C)(C)(C)(C)(C) -> Rejected: 6-valent carbon
# Parser enforces chemical physics at grammar level
```

### 5. RDKit-Independent Core
Zero dependencies for core operations. RDKit is optional for interop/testing only.

```
Core Engine (script/):  Lark parser only
Bridge Layer (rdkit_bridge.py): Optional RDKit integration
```

---

## Benchmark Results

Tested on 100-compound dataset (alkanes, rings, stereocenters, drugs, natural products):

- **99.0% round-trip success** (96/97 valid compounds)
- **Complex molecules**: Aspirin, Penicillin G, Testosterone, Glucose, Cortisol, Caffeine
- **Resonance Support**: Native handling of fused aromatic systems and mobile tautomers

```bash
python benchmark.py
# Passed: 96 | Failed: 1 (Caffeine InChI parity) | Success Rate: 99.0%
```

---

## Installation

```bash
# Core engine (RDKit-free)
pip install script-notation

# With RDKit bridge for interop
pip install script-notation[rdkit]

# Development
pip install script-notation[dev]
```

---

## Quick Start

### RDKit-Free Usage

```python
from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer

# Parse SCRIPT string to CoreMolecule (standalone graph)
parser = SCRIPTParser()
result = parser.parse("CC(=O)OC1=CC=CC=C1C(=O)O")
mol = result["molecule"]

print(f"Atoms: {len(mol.atoms)}")
print(f"Bonds: {len(mol.bonds)}")

# Canonicalize CoreMolecule to SCRIPT string
canonicalizer = SCRIPTCanonicalizer()
script_str = canonicalizer.canonicalize_core(mol)
print(f"Canonical: {script_str}")
```

### RDKit Interop

```python
from rdkit import Chem
from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT

# SMILES -> SCRIPT
mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
script_str = SCRIPTFromMol(mol)
print(f"SCRIPT: {script_str}")

# SCRIPT -> RDKit Mol
mol_back = MolFromSCRIPT(script_str)
inchi = Chem.MolToInchi(mol_back)
print(f"InChI: {inchi}")
```

### Advanced Features

```python
# Peptide mode
parser.parse("{A.G.S[A]K}")  # Ala-Gly-Ser-Lys with disulfide

# Polymer mode
parser.parse("{CC}n")  # Polyethylene
parser.parse("{CC(C1=CC=CC=C6)}100")  # 100-mer polystyrene

# Validation
parser.is_valid("CCO")  # True
parser.is_valid("C(")   # False - unclosed branch
```

---

## Project Structure

```
script-notation/
├── script/                    # Core engine (RDKit-free)
│   ├── mol.py                 # CoreMolecule graph representation
│   ├── parser.py              # Lark-based SCRIPT parser
│   ├── canonical.py           # DFS canonicalization engine
│   ├── stereo.py              # Stereochemistry perception
│   ├── cip.py                 # CIP priority calculator
│   ├── state_machine.py       # Sandhi validation rules
│   ├── grammar.lark           # SCRIPT v2.0 EBNF grammar
│   ├── ranking.py             # Morgan invariant ranking
│   ├── local_rings.py         # Topological ring resolution
│   └── rdkit_bridge.py        # Optional RDKit interop
├── docs/
│   ├── SPEC.md                # Complete SCRIPT specification
│   ├── CIP_STEREO_THEORY.md   # Stereochemistry reconciliation theory
│   └── STANDALONE_ARCHITECTURE.md
├── tests/
│   ├── test_parser.py
│   └── test_rdkit_integration.py
├── examples/
│   ├── basic_usage.py
│   └── rdkit_demo.py
├── benchmark.py               # 100-compound validation suite
└── LICENSE                    # MIT + Commons Clause
```

---

## Technical Highlights

### Grammar (EBNF)
```
molecule: chain
chain: (atom | branch | local_ring | bond)+
atom: element | "[" element charge? chiral? isotope? "]" | peptide | polymer
local_ring: "&" INT anubandha
anubandha: ":" | "."
branch: "(" chain ")"
chiral: "@" | "@@"
bond: "-" | "=" | "#" | "~" | "=:" | "/" | "\"
```

### Canonicalization Algorithm
1. Compute Morgan invariants (atomic number, degree, bond orders)
2. Rank atoms with lexicographical tie-breaking
3. DFS traversal from highest-ranked atom
4. Prioritize neighbors: Parent < [H] < Ring-Closures < Ring-Openings < Branches < Chain
5. Perceive stereochemistry using CIP-reconciled neighbor order
6. Generate canonical SCRIPT string

### Stereochemistry Reconciliation
```python
# Transform local chirality to universal CIP space
chiral_CIP = chiral_local XOR parity(local_order -> CIP_order)

# SCRIPT uses DFS order, RDKit uses atom indices
# CIP priorities provide common reference frame
# Permutation parity handles coordinate transformation
```

---

## Comparison with Existing Notations

| Feature | SMILES | SELFIES | InChI | SCRIPT |
|---------|--------|---------|-------|--------|
| Canonical | No* | No | Yes | Yes |
| Human-readable | Yes | No | No | Yes |
| Writable | Yes | No | No | Yes |
| Invalid-proof | No | Yes | N/A | Yes (Sandhi) |
| Stereochemistry | Fragile | Limited | Robust | Robust (CIP) |
| Ring notation | Global | Encoded | N/A | **Topological &N:** |
| ML-ready | Partial | Yes | No | Yes |

*SMILES has canonical forms (e.g., RDKit canonical SMILES) but multiple tools produce different canonicalizations.

---

## Limitations & Future Work

- **Organometallics**: Dative bonds supported (C->N), but complex coordination chemistry needs extension.
- **Macrocycles**: Large rings (>9 atoms) use extended lookback notation.
- **Performance**: Parser optimized for correctness; speed improvements planned for production.

---

## Contributing

SCRIPT is open for academic and non-commercial use. Contributions welcome:

1. Fork the repository
2. Create feature branch (`git checkout -b feature/improvement`)
3. Add tests for new functionality
4. Submit pull request

See `docs/SPEC.md` for grammar details and `docs/CIP_STEREO_THEORY.md` for stereochemistry theory.

---

## Citation

If you use SCRIPT in academic work, please cite:

```
Sharma, S. (2026). SCRIPT: Structural Chemical Representation in Plain Text.
A Deterministic Molecular Notation System with CIP-Based Stereochemistry.
https://github.com/script-notation/script
```

---

## License

**MIT License with Commons Clause**

Free for academic research, personal projects, and non-commercial open-source development.

Commercial use (sale of software/services whose value derives substantially from SCRIPT) requires separate licensing agreement.

See `LICENSE` for full terms.

---

## Contact

Developed by **Sangeet Sharma** and the SCRIPT team.

For questions, collaborations, or commercial licensing:
- GitHub Issues: [script-notation/script/issues](https://github.com/script-notation/script/issues)
- Documentation: See `docs/` directory

---

*"A linear script to unfold molecular complexity."*

---

PS: Sangeet's the name, a daft undergrad splashing through chemistry and code like a toddler; my titrations are a mess, and I've used my mouth to pipette.
