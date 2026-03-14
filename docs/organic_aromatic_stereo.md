# SCRIPT Documentation: Organic, Aromatic, and Stereo

This document outlines the core SCRIPT syntax for organic chemistry, topological aromaticity, and stereochemical representation.

## 1. Organic Core (The Kekule Basis)
Unlike SMILES, SCRIPT is **Kekule-only** in its base form. Aromaticity is not a property of the atom symbol (no lowercase `c`), but a property of the topological cycle (The Anubandha).

- **Symbols**: Standard IUPAC uppercase (e.g., `C`, `N`, `O`, `P`, `S`, `F`, `Cl`, `Br`, `I`).
- **Hydrogen Fill (Lopa)**: Implicit hydrogens are calculated automatically based on the standard valence of the symbol.
  - `C` -> $CH_4$
  - `CC` -> Ethane ($CH_3-CH_3$)
  - `C=C` -> Ethylene ($CH_2=CH_2$)

## 2. Aromaticity & Rings (The Anubandha System `&`)
SCRIPT solves the "ring-crossing" and "aromaticity ambiguity" of SMILES using Anubandha markers.

### Topological Ring Closures
- **`&N.`**: Aliphatic ring of size `N`.
- **`&N:`**: Aromatic/Resonant ring of size `N`.
  *Example*: `C1CCCCC1` in SMILES is `C&6.C5` or simply `C&6.` in SCRIPT.
  *Example (Benzene)*: `c1ccccc1` is `C&6:` in SCRIPT.

### Legacy Jumps (V1 Compatible)
- **`1-9`**: Jumps $N-1$ atoms back in the DFS path.
- **`%NN`**: Jumps $NN$ atoms back.

## 3. Stereochemistry (The Vak Order)
SCRIPT eliminates the need for complex CIP (Cahn-Ingold-Prelog) ranking by using the **Vak Order** (Sentence Order).

- **Rule**: Chirality is resolved based on the order in which neighbors appear in the string.
- **`@`**: Counter-Clockwise (CCW) relative to the incoming path.
- **`@@`**: Clockwise (CW) relative to the incoming path.

*Example*: `C[C@H](O)C(=O)O` (L-Lactic Acid). 
The ordering is: `[Parent, H, O, C(=O)O]`. SCRIPT resolves this sequence natively without checking atomic weights.

## 4. Double Bond Stereo (The Sandhi)
- **`/`**: Up bond.
- **`\`**: Down bond.
- *Example*: `C/C=C\C` (trans-2-butene).
