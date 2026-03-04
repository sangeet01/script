# SCRIPT: Structural Chemical Representation in Plain Text

**SCRIPT** is a next-generation molecular representation language and a high-performance, RDKit-independent cheminformatics engine designed for the AI era.

Developed by **Sangeet Sharma**, SCRIPT addresses the fundamental flaws in legacy notations like SMILES and SELFIES by embedding chemical physics directly into its grammar and generative state machine.

---

## 🚀 Key Innovations

### 1. "Sandhi" (Valence Guards)
Inspired by Paninian linguistic principles, **Sandhi** implements real-time valence monitoring during parsing. If a SCRIPT string attempts to violate chemical laws (e.g., a 5-valent carbon), the engine either flags the error or applies rule-based corrections, ensuring the output is always a physically valid `CoreMolecule`.

### 2. "Lopa" (Elision Rules)
Deterministic elision (**Lopa**) of implicit hydrogens and bond orders enables highly concise strings without sacrificing the ability to reconstruct the full, unambiguous molecular graph.

### 3. Integral Parity (Deterministic Stereo)
SCRIPT solves the "black-box" stereochemistry problem. Instead of relying on complex perception algorithms, SCRIPT links stereochemical parity directly to the **DFS traversal order** of the canonicalizer.
- **Priority**: Parent < [H] < Ring-Closures < Ring-Openings < Branches < Chain.
- **Result**: 100% predictable stereochemistry anchored to the string's own topology.

### 4. RDKit-Independent Core
The SCRIPT library is architected to be zero-dependency for its core operations (`Parser`, `Graph`, `Canonicalizer`, `Stereo Engine`).

---

## 📊 Performance & Robustness
Tested against the ChEMBL dataset and medicinal chemistry benchmarks:
- **Success Rate**: 97% round-trip parity with RDKit-generated InChIs.
- **Stereo Alignment**: Successfully resolves complex cases like **Glucose** and symmetric fused rings using the **CIP Priority Reconciliation** theory.
- **Speed**: Optimized DFS traversal with integer-based bond representation.

---

## 📂 Project Structure

```text
open/
├── script/           # Core RDKit-free engine
│   ├── mol.py        # Standalone graph representation
│   ├── parser.py     # Lark-based high-performance parser
│   ├── canonical.py  # DFS-based canonicalization
│   └── stereo.py     # Standalone parity & CIP engine
├── docs/             # Technical Specifications & Theory
│   ├── SPEC.md       # SCRIPT v1.0 Grammar & Rules
│   └── CIP_THEORY.md # Coordinate-free stereo reconciliation
├── tests/            # Validation & Round-trip suites
├── examples/         # Usage demonstrations
└── LICENSE           # MIT + Commons Clause
```

---

## 🛠️ Quick Start (RDKit-Free)

```python
from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer

# 1. Parse a SCRIPT string to a CoreMolecule
parser = SCRIPTParser()
result = parser.parse("O[C@H](C)[C@@H](O)C")
mol = result["molecule"] # A standalone graph object

# 2. Canonicalize the core molecule
canonicalizer = SCRIPTCanonicalizer()
script_string = canonicalizer.canonicalize_core(mol)
print(script_string) # Deterministic, canonical SCRIPT
```

---

## ⚖️ License
This project is licensed under the **MIT License with Commons Clause Amendment**. 
It is free for use in academic research, personal projects, and non-commercial open-source development. **Commercial use, including sale or use in paid products/services, is prohibited without a separate agreement.**

---
*Developed by Sangeet Sharma (2026). Designed for safe, concise, and intelligent chemical processing.*
