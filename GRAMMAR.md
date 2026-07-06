# SCRIPT Grammar Guide

> A linear script to unfold molecular complexity.

SCRIPT is a molecular notation built on Paninian Sanskrit grammar principles. This guide explains the complete grammar in plain language. The formal LALR(1) grammar is defined in `script/grammar.lark`. 

---

## Table of Contents

1. [Basic Structure](#1-basic-structure)
2. [Atoms](#2-atoms)
3. [State Blocks](#3-state-blocks)
4. [Bonds](#4-bonds)
5. [Branches](#5-branches)
6. [Rings](#6-rings)
7. [Stereochemistry](#7-stereochemistry)
8. [Multi-Component Systems](#8-multi-component-systems)
9. [Reactions](#9-reactions)
10. [Materials Science](#10-materials-science)
11. [Quantum States](#11-quantum-states)
12. [Polymers](#12-polymers)
13. [Biopolymers](#13-biopolymers)
14. [Query Atoms](#14-query-atoms)
15. [Periodic Topology](#15-periodic-topology)
16. [Grammar Philosophy](#16-grammar-philosophy)
17. [Comparison with SMILES and SELFIES](#17-comparison)

---

## 1. Basic Structure

Every SCRIPT string represents a molecule, reaction, material, or polymer. The parser reads left-to-right, building atoms and bonds as it goes.

```
CCO              -> Ethanol
C(=O)O           -> Carboxylic acid group
C:C:C:C:C:C&6:   -> Benzene (aromatic ring)
CC>>CCO          -> Reaction: ethylene to ethanol
```

A SCRIPT string can contain:
- A single molecule: `CCO`
- A multi-fragment salt: `CN(C)C(=N)N=C(N)N.[Cl-]`
- A reaction: `CC>>CCO`
- A material: `[[Rutile]] Ti(O)2`
- A polymer: `{[CC]}n`
- A peptide: `{A.G.S}`

---

## 2. Atoms

### Organic Subset (no brackets needed)

These elements can be written without brackets, just like SMILES:

```
B, C, N, O, P, S, F, Cl, Br, I
```

Examples:
```
C            -> Carbon (methane with implicit hydrogens)
N            -> Nitrogen (ammonia with implicit hydrogens)
O            -> Oxygen (water with implicit hydrogens)
Cl           -> Chlorine
```

### Bracket Atoms (full periodic table)

Use brackets `[...]` for any element not in the organic subset, or when you need to specify charge, isotope, stereochemistry, or hydrogen count. All 118 elements are supported.

```
[Fe]         -> Iron
[Pt]         -> Platinum
[Ir]         -> Iridium
[U]          -> Uranium
[Og]         -> Oganesson
[Na+]        -> Sodium cation
[Cl-]        -> Chloride anion
[13C]        -> Carbon-13 isotope
[2H]         -> Deuterium
[CH3]        -> Carbon with 3 explicit hydrogens
[C@H]        -> Chiral carbon with 1 hydrogen
[O.]         -> Oxygen radical
[C..]        -> Carbon diradical (two unpaired electrons)
```

### Wildcard

```
*            -> Wildcard atom (scaffold attachment point)
[*]          -> Wildcard in brackets
```

---

## 3. State Blocks

State blocks use angle brackets `<...>` to carry metadata on atoms. Multiple attributes can be comma-separated inside a single pair of brackets.

### Available state attributes

| Syntax | Meaning | Example |
|--------|---------|---------|
| `<13>` | Isotope mass number | `C<13>` = Carbon-13 |
| `<+>` | Positive charge | `N<+>` = N cation |
| `<->` | Negative charge | `O<->` = O anion |
| `<2+>` | Charge +2 | `Mg<2+>` = Mg²⁺ |
| `<h3>` | Explicit H count | `C<h3>` = 3 hydrogens |
| `<~0.9>` | Fractional occupancy (alloys) | `Ti<~0.9>` = 90% Ti |
| `<s:3>` | Spin multiplicity | `O<s:3>` = triplet |
| `<*>` | Excited state | `[C<*>]` = excited carbon |
| `<m>` | Tautomer/mobility flag | `C<m>` = mobile tautomer |
| `<sqp>` | Coordination geometry | `Pt<sqp>` = square planar |

### Combined state blocks

```
C<13,+>         -> Carbon-13 cation
C<13,+,s:3>     -> Carbon-13 cation in triplet state
Ti<~0.9>        -> Titanium with 90% occupancy
```

---

## 4. Bonds

SCRIPT uses 9 typed bond types. Single bonds are usually implicit (just write two atoms next to each other).

### Bond Reference

| Symbol | Name | Use case |
|--------|------|----------|
| `-` | Single | Standard covalent (usually implicit) |
| `=` | Double | C=O, C=C |
| `#` | Triple | C#C, C#N |
| `:` | Aromatic | Benzene ring bonds |
| `->` | Dative | Donor to acceptor (N->B) |
| `<-` | Reverse dative | Acceptor from donor (B<-N) |
| `=:` | Tautomeric | Mobile tautomeric bond |
| `>` | Coordinate | Metal coordination (Fe>C) |
| `*N` | Haptic (eta-N) | Ferrocene (*5 = eta-5) |

### E/Z Stereo Bonds

```
/            -> Up bond (E/Z stereochemistry)
\            -> Down bond (E/Z stereochemistry)
```

### Examples

```
CC           -> Ethane (implicit single bond)
C=C          -> Ethene (double bond)
C#C          -> Ethyne (triple bond)
C:C          -> Aromatic bond
N->B(F)(F)F  -> Dative bond: N donates to B
B<-N(F)(F)F  -> Reverse dative: B accepts from N
Fe*5C:C:C:C:C -> Ferrocene (eta-5 haptic bond)
C/C=C/C      -> E-2-butene (trans)
C/C=C\C      -> Z-2-butene (cis)
C=:C-C       -> Tautomeric bond
```

---

## 5. Branches

Use parentheses `(...)` for side chains. Branches can be nested.

```
CC(C)C       -> Isobutane
CC(O)C       -> Isopropanol
CC(C)(C)C    -> Neopentane
C(=O)O       -> Carboxylic acid group
C(C(=O)O)N   -> Amino acid backbone
```

A bond before the branch content specifies how the branch connects:

```
C(=O)O       -> Carbon double-bonded to O, single-bonded to O
C(-O)N       -> Carbon single-bonded to O (explicit), single-bonded to N
```

---

## 6. Rings

### Anubandha Ring System (V2 — recommended)

Use `&N` where N is the ring size (number of atoms in the ring). The suffix specifies the bond type:

| Suffix | Meaning |
|--------|---------|
| `:` | Aromatic ring |
| `-` | Aliphatic ring |

```
C:C:C:C:C:C&6:    -> Benzene (6-membered aromatic ring)
CCCCCC&6-         -> Cyclohexane (6-membered aliphatic ring)
CCCCC&5-          -> Cyclopentane (5-membered aliphatic ring)
C:C:C:C:C:N&6:    -> Pyridine (6-membered aromatic with N)
```

The number N is the **ring size** — the number of atoms in the ring path through the DFS tree. Both the encoder and parser compute N the same way (depth difference), so they always agree on which atom the ring closes to.

### Legacy Ring Closures (SMILES-compatible)

For compatibility with SMILES, bare digits and `%NN` are also supported:

```
C1CCCCC1     -> Cyclohexane (SMILES-style: digit 1 opens and closes)
C1CCC1       -> Cyclobutane
C%12CCC%12   -> 12-membered ring (for ring sizes > 9)
```

### Named Rings (polycyclic bridges)

```
C[A]CCC[A]       -> Ring closure using register A
C[A]C[B]C[A]C[B] -> Two rings sharing atoms
```

---

## 7. Stereochemistry

SCRIPT supports 7 stereochemistry types with 25+ markers. This is the most complete stereo system of any molecular notation.

### Tetrahedral

| Marker | Meaning |
|--------|---------|
| `@` | Tetrahedral R (DFS-relative — the "Sandhi" form) |
| `@@` | Tetrahedral S (DFS-relative) |
| `@R` | CIP-absolute R (frame-independent — the "Sthiti" form) |
| `@S` | CIP-absolute S (frame-independent) |
| `@r` | Pseudoasymmetric r (auto-detected) |
| `@s` | Pseudoasymmetric s (auto-detected) |

The `@R`/`@S` markers are CIP-absolute: they mean the same thing regardless of atom ordering. The canonical form always emits `@R`/`@S` (never `@`/`@@`), making it idempotent.

```
[C@H](F)(Cl)Br        -> Chiral bromochlorofluoromethane
C[C@H](O)C(=O)O       -> L-Lactic acid
C[C@@H](O)C(=O)O      -> D-Lactic acid
[C@R](F)(Cl)(Br)I     -> CIP-absolute R
[C@S](F)(Cl)(Br)I     -> CIP-absolute S
```

### Pseudoasymmetric (auto-detected)

A center with enantiomorphic substituents (same topology, opposite embedded R/S) is auto-detected and emitted as `@r` or `@s` (lowercase). No manual specification needed.

### Axial (allenes, biaryls)

| Marker | Meaning |
|--------|---------|
| `@AX1` | Axial CW (geometry-resolved) |
| `@AX2` | Axial CCW (geometry-resolved) |
| `@AX` | Axial unspecified (no 3D coords) |

```
C[C@AX1]=C=CC    -> Allenic axial chirality, CW
C[C@AX2]=C=CC    -> Allenic axial chirality, CCW
```

### Square Planar

| Marker | Meaning |
|--------|---------|
| `@SP1` | cis |
| `@SP2` | trans |
| `@SP` | Unspecified |

```
[Pt@SP1](Cl)([NH3])(Cl)([NH3])   -> cisplatin (cis)
[Pt@SP2](Cl)([NH3])([NH3])(Cl)   -> transplatin (trans)
```

### Octahedral

| Marker | Meaning |
|--------|---------|
| `@OH1` | fac (Delta) |
| `@OH2` | fac (Lambda) |
| `@OH3`–`@OH5` | mer variants |
| `@OH` | Unspecified |

```
[Co@OH1](Cl)(Cl)(Cl)(N)(N)N   -> fac-Co complex
```

### Trigonal Bipyramidal

| Marker | Meaning |
|--------|---------|
| `@TB1` | Apical substitution variant 1 |
| `@TB2` | Apical substitution variant 2 |
| `@TB` | Unspecified |

```
[P@TB1](Cl)(Cl)(Cl)(Cl)Cl   -> Trigonal bipyramidal P
```

### Pyramidal

| Marker | Meaning |
|--------|---------|
| `@PY1` | CW |
| `@PY2` | CCW |
| `@PY` | Unspecified |

```
[N@PY1](C)(C)C   -> Pyramidal nitrogen
```

### Planar (metallocenes)

| Marker | Meaning |
|--------|---------|
| `@PL1` | Rp |
| `@PL2` | Sp |
| `@PL` | Unspecified (auto-detected) |

Planar chirality in ferrocenes is auto-detected when a transition metal is bonded to two aromatic 5-rings with asymmetric substitution. The 2D graph cannot distinguish Rp from Sp, so users must specify `@PL1`/`@PL2` to distinguish enantiomers.

### E/Z Double Bond Stereochemistry

```
C/C=C/C      -> E-2-butene (trans)
C/C=C\C      -> Z-2-butene (cis)
Cl/C=C/Cl    -> E-1,2-dichloroethene
Cl/C=C\Cl    -> Z-1,2-dichloroethene
```

### Sulfoxide Stereochemistry (Lopa rule)

Sulfoxides and other 3-coordinate chiral centers (S, N, P with lone pairs) automatically get a lone-pair "ghost" neighbor as the 4th coordination position:

```
C[S@](=O)C    -> Methyl sulfoxide (R)
C[S@@](=O)C   -> Methyl sulfoxide (S)
```

---

## 8. Multi-Component Systems

### Salts and Solvates

Use `.` to separate disconnected components (salts, solvates, mixtures):

```
[Na+].[Cl-]              -> Sodium chloride
CN(C)C(=N)N=C(N)N.[Cl-]  -> Metformin HCl
CCO.O                     -> Ethanol in water
```

### Ionic Pairs

Use `~` for ionic associations (the components are paired but not covalently bonded):

```
[Na+]~[Cl-]              -> Sodium chloride (ionic pair)
[Fe+3]~[Cl-]~[Cl-]~[Cl-] -> Iron(III) chloride complex
```

---

## 9. Reactions

Reactions use arrows to separate reactants, agents, and products.

### Reaction Arrows

| Arrow | Meaning |
|-------|---------|
| `>>` | Irreversible |
| `=>` | Irreversible (alternative) |
| `<=>` | Equilibrium |

### Examples

```
CC>>CCO                        -> Simple reaction
CC>[Pd]>CCO                    -> Reaction with Pd catalyst
[C:1]=O>>[C:1]O                -> Atom-mapped reaction (C atom 1)
C=O<=>CO                       -> Equilibrium reaction
CC>>C.CO                       -> Multi-product reaction
```

Agent notation: `reactants > agents > products` (the middle section between single `>` arrows contains catalysts/solvents).

---

## 10. Materials Science

### Crystallographic Context

Use `[[...]]` to specify crystal structure or space group:

```
[[Rutile]] Ti(O)2                          -> TiO2 in rutile phase
[[Anatase]] Ti(O)2                         -> TiO2 in anatase phase
[[bcc]] Fe                                 -> Body-centered cubic iron
[[fcc]] Fe                                 -> Face-centered cubic iron
[[Rutile;4.593,4.593,2.959,90,90,90]] Ti(O)2  -> With lattice parameters (a,b,c,alpha,beta,gamma)
```

### Surfaces and Interfaces

Use `|` to separate phases (surface adsorption, battery interfaces):

```
[[Pt_111]] | >C=O       -> CO adsorbed on Pt(111) surface
[[LiCoO2]] | Li<+>      -> Li+ in LiCoO2 battery lattice
```

### Alloys

Use `<~FLOAT>` for fractional occupancy:

```
Ti<~0.9>N<~0.1>        -> Doped titanium nitride (90% Ti, 10% N)
Fe<~0.5>Ni<~0.5>       -> 50-50 iron-nickel alloy
```

---

## 11. Quantum States

### Spin Multiplicity

```
O=O<s:3>       -> Triplet oxygen (ground state, spin=3)
O=O<s:1>       -> Singlet oxygen
O=O<s:1,*>     -> Singlet oxygen (excited)
```

### Excited States

```
[C<*>]         -> Excited carbon
```

---

## 12. Polymers

### Simple Polymers

Use `{[...]}` with repeat specifications:

```
{[CC]}n                  -> Polyethylene (unspecified length)
{[CC]}<n:50>             -> Exact 50-mer
{[CC]}<n:50-100>         -> Stochastic range (50-100 units)
{[CC(C)]}n               -> Polypropylene
{[CC(=O)O]}n             -> Polyester
```

### Block Copolymers

Use junction keywords to connect polymer blocks:

| Junction | Meaning |
|----------|---------|
| `-b-` | Diblock |
| `-alt-` | Alternating |
| `-a-` | Alternating (alias) |
| `-ran-` | Random |
| `-r-` | Random (alias) |
| `-stat-` | Statistical |
| `-g-` | Graft |

```
{[CC]}<n:50> -b- {[CCCO]}<n:100>              -> Diblock copolymer
{[CC]}<n:50> -b- {[CCCO]}<n:30> -b- {[CC]}<n:20>  -> Triblock
{[CC]} -alt- {[CC(N)C(=O)O]}                    -> Alternating copolymer
{[CC]} -ran- {[N]}                              -> Random copolymer
{[CC]} -g- {[CCCO]}                             -> Graft copolymer
{[CC]} -stat- {[N]}                             -> Statistical copolymer
```

---

## 13. Biopolymers

### Peptides

Use curly braces `{...}` with dot-separated monomer codes:

```
{A}              -> Alanine
{A.G.S}          -> Ala-Gly-Ser tripeptide
{A.G.S.K}        -> Tetrapeptide
{pS.G.acK.V}     -> PhosphoSer-Gly-acetylLys-Val
```

### Amino Acid Codes

Single-letter: A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V

Three-letter: Ala, Arg, Asn, Asp, Cys, Glu, Gln, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val

### Post-Translational Modifications

| Code | Modification |
|------|-------------|
| `pS`, `pT`, `pY` | Phosphorylation |
| `acK` | Acetylation |
| `mK`, `mR` | Methylation |
| `ubK` | Ubiquitination |
| `suK` | Sumoylation |
| `oxM` | Methionine oxidation |
| `nitY` | Tyrosine nitration |
| `Hyp`, `Hyl` | Hydroxylation |
| `Sec` | Selenocysteine |
| `Pyl` | Pyrrolysine |
| `Orn`, `Cit` | Non-standard residues |
| `Gla`, `Dpr`, `Dab` | Specialized residues |
| `Sar` | Sarcosine |

### Nucleic Acids

```
{dA.dG.dC.dT}     -> DNA strand
{rA.rG.rC.rU}     -> RNA strand
```

### Epigenetic Modification Codes

| Code | Modification |
|------|-------------|
| `m5C` | 5-methylcytosine (DNA methylation) |
| `m6A` | N6-methyladenine |
| `hm5C` | 5-hydroxymethylcytosine |
| `f5C` | 5-formylcytosine |
| `ca5C` | 5-carboxylcytosine |
| `m1A` | 1-methyladenine (RNA) |
| `m1G` | 1-methylguanosine |
| `m2G` | N2-methylguanosine |
| `m22G` | N2,N2-dimethylguanosine |
| `m7G` | 7-methylguanosine (RNA cap) |
| `psU` | Pseudouridine |
| `s4U` | 4-thiouridine |
| `I` | Inosine |
| `diU` | Dihydrouridine |

```
{m5C.dG.dA.dT}     -> Methylated DNA
{m6A.rG.rC.rU}     -> m6A RNA modification
{psU.rG.rC.rA}     -> Pseudouridine-containing RNA
```

---

## 14. Query Atoms

SMARTS-style query atoms for substructure search. These are pattern atoms (read-only), not canonical representations.

| Pattern | Meaning |
|---------|---------|
| `[#6]` | Atomic number 6 (carbon) |
| `[R]` | Any ring atom |
| `[r5]` | Atom in 5-membered ring |
| `[v3]` | Valence = 3 |
| `[X2]` | Connectivity (total degree) = 2 |
| `[!N]` | Not nitrogen |
| `[a]` | Aromatic atom |
| `[A]` | Aliphatic atom |
| `[#6,#7]` | Carbon OR nitrogen (comma = OR) |
| `[#6;R]` | Carbon AND ring atom (semicolon = AND) |

```
[#6]CC        -> Any carbon followed by CC
[R]CC         -> Any ring atom followed by CC
[!N]CC        -> Any atom except nitrogen
[#6,#7]C      -> Carbon or nitrogen
```

---

## 15. Periodic Topology

For periodic structures (MOFs, zeolites, coordination polymers), bonds can carry a lattice translation vector `@tx,ty,tz`:

```
Ti->@0,0,1O    -> Dative bond crossing +c unit cell (MOF pillar)
Fe-@1,0,0Fe   -> Single bond crossing +a unit cell boundary
C-@0,0,0C     -> Intracell bond (same as plain C-C)
```

Lattice parameters are specified in the macroscopic context:

```
[[Rutile;4.593,4.593,2.959,90,90,90]] Ti(O)2
     a=4.593Å  b=4.593Å  c=2.959Å  alpha=90° beta=90° gamma=90°
```

---

## 16. Grammar Philosophy

SCRIPT is built on Paninian Sanskrit grammar principles. The Ashtadhyayi (4th century BCE) is the most rigorous formal grammar system ever devised for any language. Nine concepts map directly to molecular notation:

| Sanskrit | SCRIPT meaning |
|----------|---------------|
| **Dhatu** (root) | Atomic element — the irreducible root |
| **Pratyaya** (marker) | Stereo suffix (@R, @SP1, @AX2) |
| **Vak** (speech order) | DFS neighbor ordering |
| **Sandhi** (junction) | Valence reconciliation — bonds must satisfy both atoms |
| **Lopa** (elision) | Ghost neighbors — implicit H, lone pairs |
| **Sthiti** (stable form) | CIP-absolute @R/@S — frame-independent |
| **Yatha-sankhya** (order) | Pseudoasymmetric @r/@s — secondary markers |
| **Sthana** (place) | Geometry classes — tetrahedral, planar, etc. |
| **Vikarana** (expansion) | Peptide/nucleic acid monomer expansion |
| **Ekadesa** (unification) | Aromatic `:` as unified bond |

### Three design principles

1. **Paninian Grammar** — generative validation through Sandhi (valence rules enforced at parse time, not after construction)

2. **Graph Linearization** — DFS traversal for canonicalization, local ring encoding (no global state), one molecule = one canonical string

3. **Invalid-Proof Design** — generative state machine prevents valence violations, every valid SCRIPT string maps to a valid molecule, parse errors caught immediately

---

## 17. Comparison

### SMILES vs SELFIES vs SCRIPT

| Feature | SMILES | SELFIES | SCRIPT |
|---------|--------|---------|--------|
| Canonical form | No (toolkit-dependent) | No | **Yes (deterministic)** |
| Human-readable | Yes | No | **Yes** |
| Hand-writable | Yes | No | **Yes** |
| Invalid-proof | No | Yes (construction) | Yes (constraint) |
| Grammar type | Ad hoc | CFG (Type 2) | **LALR(1)** |
| Stereo: tetrahedral | Fragile | Limited | **CIP-absolute (@R/@S)** |
| Stereo: non-tetrahedral | Partial | No | **7 types** |
| Ring closures | Global labels | Derivation state | **Ring size (&N)** |
| Bond types | Integer | Integer | **Typed enum (9 types)** |
| All 118 elements | Partial | Partial | **Yes** |
| Materials/alloys | No | No | **Yes** |
| Crystallography | No | No | **Yes** |
| Surfaces | No | No | **Yes** |
| Quantum states | No | No | **Yes** |
| Polymers | No | No | **Yes** |
| Biopolymers + PTMs | No | No | **Yes** |
| Reactions | Partial | No | **Yes** |
| Query atoms | SMARTS only | No | **Yes** |
| RDKit-free core | No | No | **Yes** |
| ML robustness | Low | **100% (construction)** | **100% (constraint)** |

### Ring notation comparison

```
SMILES:   C1CCCCC1     (global label — must track open ring 1)
SELFIES:  [C][C][C][C][C][Ring1][Branch1]  (derivation state)
SCRIPT:   CCCCCC&6-    (ring size = 6, local, no global state)
```

### Aromaticity comparison

```
SMILES:   c1ccccc1     (lowercase = aromatic by convention)
SCRIPT:   C:C:C:C:C:C&6:  (explicit : bonds + &6: ring closure)
```

---

## Further Reading

- `script/grammar.lark` — The complete LALR(1) grammar 
- `docs/SPEC.md` — Complete technical specification
- `docs/CIP_STEREO_THEORY.md` — Stereochemistry reconciliation theory
- `docs/STANDALONE_ARCHITECTURE.md` — Implementation architecture

---

*"A linear script to unfold molecular complexity."*
