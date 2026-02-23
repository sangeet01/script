# SCRIPT: Structural Chemical Representation In Plain Text



**SCRIPT** is a next-generation molecular notation system designed to replace SMILES with a robust, deterministic, and machine-learning-ready "One True String." By combining the algorithmic rigor of Panini's Sanskrit grammar with state-of-the-art graph theory, SCRIPT eliminates the ambiguities that have plagued chemical informatics for 35 years.

---

## Key Advantages

*   **Robust by Design**: Every SCRIPT string is a valid molecule. Our generative state machine "Sandhi rules" prevent valence violations before they happen.
*   **Deterministic Canonicalization**: Built on RDKit's Morgan Invariants, SCRIPT ensures that every molecule has exactly one unique string representation—no more "SMILES mismatch" hell.
*   **Local Ring Encoding**: Replaces arbitrary ring labels with "Lookback Counting" (e.g., `C1CCCCC6`). Connections are intrinsic to the string's memory, eliminating unclosed ring errors.
*   **Dual-Mode Support**: Seamlessly transition between atomic detail and monomer-level shortcuts (Peptides, Polymers, PTMs) in a single readable line.
*   **Audit-Ready**: Advanced error-correction logs every structural "morph" to maintain scientific integrity while ensuring AI compatibility.

## Theoretical Foundations

SCRIPT is not just another notation; it is a **Generative State Machine** based on:
1.  **Panini’s Sanskrit Grammar**: Using context-aware rules (*Anuvritti*) and junction logic (*Sandhi*) to ensure internal consistency.
2.  **Graph Linearization**: A unique DFS-traversal that prioritizes longest-chain fidelity and lexical back-indexing.
3.  **SELFIES-stye Robustness**: Every coordinate in SCRIPT-space maps to a physical reality.

## Current Implementation Status

Our production-ready Python engine currently supports:
- [x] **Canonicalizer**: Morgan-based ranking with lexicographical tie-breaking.
- [x] **Parser**: High-performance Lark-based grammar with local ring resolution.
- [x] **Validation**: Verified against a **100-compound benchmark** (Aspirin, Penicillin, complex drugs) with **90% round-trip fidelity**.
- [ ] **State Machine Validator**: (In progress) For 100% "invalid-proof" generation.

## Features at a Glance

| Feature | SCRIPT Notation | Classic SMILES |
| :--- | :--- | :--- |
| **Benzene** | `C=CC=CC=C6` | `c1ccccc1` (Ambiguous labels) |
| **Local Rings** | `C1CCCCC6` | `C1CCCCC1` (Global labels) |
| **Peptides** | `{A.G.S[A]K}` | `N[C@@H](C)C(=O)...` (Verbose) |
| **Dative Bonds** | `C->N` | Not standard |

---

## Quick Start

```python
from script.canonical import SCRIPTCanonicalizer
from script.parser import SCRIPTParser

# To SCRIPT
c = SCRIPTCanonicalizer()
script_str = c.canonicalize_smiles("CC(=O)Oc1ccccc1C(=O)O")
# Output: CC(=O)OC1=CC=CC=C1C(=O)O

# Back to Mol
p = SCRIPTParser()
mol = p.parse_to_mol(script_str)
```

## License
Open Source. (MIT/Apache 2.0)

---
*Developed by Sangeet Sharma and the SCRIPT Team.*
