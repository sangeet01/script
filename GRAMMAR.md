# SCRIPT Grammar Guide

A human-readable explanation of the SCRIPT notation system.

---

## Basic Structure

Every SCRIPT string represents a molecule, reaction, or material. The parser reads left-to-right, building atoms and bonds as it goes.

```
CCO          -> Ethanol (carbon-carbon-oxygen)
CC(O)C       -> Isopropanol (carbon-carbon with oxygen branch)
C=O          -> Formaldehyde (carbon double-bonded to oxygen)
```

---

## Atoms

### Simple Atoms (Organic Subset)

These don't need brackets:

```
B, C, N, O, P, S, F, Cl, Br, I
```

Examples:
```
C            -> Carbon (methane with implicit hydrogens)
N            -> Nitrogen (ammonia with implicit hydrogens)
O            -> Oxygen (water with implicit hydrogens)
```

### Bracket Atoms (Full Periodic Table)

Use brackets `[...]` for:
- Charges: `[Na+]`, `[O-]`, `[NH4+]`
- Isotopes: `[13C]`, `[2H]`, `[18O]`
- Explicit hydrogens: `[CH3]`, `[NH2]`
- Transition metals: `[Fe]`, `[Pt]`, `[Cu]`
- Stereochemistry: `[C@H]`, `[C@@H]`

Examples:
```
[Na+].[Cl-]  -> Sodium chloride (salt)
[13C]        -> Carbon-13 isotope
[C@H](O)C    -> Chiral carbon (lactic acid)
```

### State Blocks (V3 Materials)

Use angle brackets `<...>` for advanced properties:

```
<13>         -> Isotope mass 13
<+>          -> Positive charge
<->          -> Negative charge
<h3>         -> 3 explicit hydrogens
<~0.9>       -> 90% occupancy (alloys)
<s:3>        -> Spin triplet state
<*>          -> Excited state
```

Examples:
```
Ti<~0.9>N<~0.1>    -> Titanium nitride alloy (90% Ti, 10% N)
O<s:3>             -> Triplet oxygen (ground state)
C<13,+>            -> Carbon-13 cation
```

---

## Bonds

### Basic Bonds

```
-            -> Single bond (usually implicit)
=            -> Double bond
#            -> Triple bond
:            -> Aromatic bond
```

Examples:
```
C-C          -> Ethane (single bond, usually written CC)
C=C          -> Ethene (double bond)
C#C          -> Ethyne (triple bond)
C:C          -> Aromatic bond (benzene)
```

### Special Bonds

```
->           -> Dative bond (coordinate)
<-           -> Reverse dative
>            -> Coordinate bond (surface adsorption)
=:           -> Tautomeric/mobile bond
*            -> Resonant bond
*5           -> Haptic eta-5 bond (ferrocene)
/            -> Stereochemistry up
\            -> Stereochemistry down
```

Examples:
```
C->N         -> Dative bond (borane-amine complex)
[Fe]*5CCCCC  -> Ferrocene (eta-5 haptic bond)
C/C=C/C      -> E-alkene (trans)
C/C=C\C      -> Z-alkene (cis)
```

---

## Branches

Use parentheses `(...)` for side chains:

```
CC(C)C       -> Isobutane (two carbons, branch with carbon, one carbon)
CC(O)C       -> Isopropanol (carbon-carbon with OH branch)
CC(C)(C)C    -> Neopentane (central carbon with three methyl branches)
```

---

## Rings

### Local Ring Notation (Lookback Counting)

The digit tells you how many atoms back to connect:

```
C1CCCCC6     -> Cyclohexane (connect back 6 atoms)
C1CCC4       -> Cyclobutane (connect back 4 atoms)
```

### Anubandha Ring System (V2)

Use `&N:` for aromatic or `&N.` for aliphatic:

```
C1=CC=CC=C&6:    -> Benzene (6-membered aromatic ring)
C1CCCCC&6.       -> Cyclohexane (6-membered aliphatic ring)
```

The number after `&` is the ring size, `:` means aromatic, `.` means aliphatic.

### Named Rings (Polycycles)

Use `[A]`, `[B]`, etc. for complex ring systems:

```
C[A]CCC[A]       -> Ring closure using register A
C[A]C[B]C[A]C[B] -> Two rings sharing atoms
```

---

## Stereochemistry

### Tetrahedral Centers

```
@            -> Counter-clockwise (R in most cases)
@@           -> Clockwise (S in most cases)
```

Examples:
```
C[C@H](O)C(=O)O     -> L-Lactic acid
C[C@@H](O)C(=O)O    -> D-Lactic acid
```

### Double Bond Geometry

```
/            -> Up bond
\            -> Down bond
```

Examples:
```
C/C=C/C      -> E-2-butene (trans)
C/C=C\C      -> Z-2-butene (cis)
```

### Other Geometries

```
@SP          -> Square planar
@OH          -> Octahedral
@TB          -> Trigonal bipyramidal
@PY          -> Pyramidal
```

---

## Multi-Component Systems

### Salts and Mixtures

Use `.` to separate disconnected components:

```
[Na+].[Cl-]  -> Sodium chloride
CCO.O        -> Ethanol in water
```

### Ionic Complexes

Use `~` for ionic associations:

```
[Fe+3]~[Cl-]~[Cl-]~[Cl-]  -> Iron(III) chloride complex
```

---

## Reactions

Use `>>` or `=>` to show transformations:

```
CC(=O)O >> CCO + C=O     -> Acetic acid decomposition
C=C + [H-H] => CC        -> Ethene hydrogenation
```

---

## Peptides

Use curly braces `{...}` with single-letter amino acid codes:

```
{A}          -> Alanine
{A.G.S}      -> Ala-Gly-Ser tripeptide
{C[A].G.C[A]} -> Disulfide bridge between two cysteines
```

### Amino Acid Codes

Single letter: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y

Three letter: Ala, Arg, Asn, Asp, Cys, Glu, Gln, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val

### Post-Translational Modifications

```
pS, pT, pY   -> Phosphorylation
mK, mR       -> Methylation
acK          -> Acetylation
ubK          -> Ubiquitination
```

Example:
```
{A.pS.G}     -> Ala-phosphoSer-Gly
```

---

## Polymers

Use `{[...]}` with repeat specifications:

```
{[CC]}n              -> Polyethylene (unspecified length)
{[CC]}100            -> 100-mer polyethylene
{[CC]}<n:50-100>     -> Stochastic 50-100 unit polymer
```

### Block Copolymers

```
-b-          -> Block junction
-r-          -> Random junction
-a-          -> Alternating junction
```

Example:
```
{[CC]}<n:50> -b- {[CC(C)]}}<n:50>  -> Block copolymer
```

---

## Materials Science (V3)

### Crystallographic Context

Use `[[...]]` to specify crystal structure:

```
[[Rutile]] Ti(O)2    -> TiO2 in rutile phase
[[Anatase]] Ti(O)2   -> TiO2 in anatase phase
[[bcc]] Fe           -> Body-centered cubic iron
[[fcc]] Fe           -> Face-centered cubic iron
```

### Surface Chemistry

Use `|` to separate phases:

```
[[Pt_111]] | >C=O    -> CO adsorbed on Pt(111) surface
[[LiCoO2]] | Li<+>   -> Li+ in LiCoO2 lattice
```

### Alloys

Use `<~FLOAT>` for fractional occupancy:

```
Ti<~0.9>N<~0.1>      -> 90% Ti, 10% N alloy
Fe<~0.5>Ni<~0.5>     -> 50-50 iron-nickel alloy
```

### Electronic States

```
<s:3>        -> Spin triplet
<s:1>        -> Spin singlet
<*>          -> Excited state
```

Examples:
```
O=O<s:3>     -> Triplet oxygen (ground state)
O=O<s:1,*>   -> Singlet oxygen (excited state)
```

---

## Complete Examples

### Simple Molecules

```
C            -> Methane
CC           -> Ethane
CCO          -> Ethanol
CC(=O)O      -> Acetic acid
C#N          -> Hydrogen cyanide
```

### Rings

```
C1CCCCC6     -> Cyclohexane
C1=CC=CC=C6  -> Benzene (Kekule form)
C1CCC2CCCCC2CC1  -> Decalin (fused rings)
```

### Drugs

```
CC(=O)Oc1ccccc1C(=O)O           -> Aspirin
CN1C=NC2=C1C(=O)N(C(=O)N2C)C    -> Caffeine
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O   -> Ibuprofen
```

### Stereochemistry

```
C[C@H](O)C(=O)O      -> L-Lactic acid
C[C@@H](O)C(=O)O     -> D-Lactic acid
F[C@H](Cl)Br         -> Chiral halomethane
```

### Materials

```
[[Rutile]] Ti(O)2              -> Rutile TiO2
Ti<~0.9>N<~0.1>                -> TiN alloy
[[Pt_111]] | >C=O              -> CO on Pt surface
O=O<s:3>                       -> Triplet oxygen
```

### Peptides

```
{A.G.S}                        -> Ala-Gly-Ser
{C[A].G.S.C[A].K}              -> Peptide with disulfide bridge
{A.pS.G}                       -> Phosphorylated serine peptide
```

### Polymers

```
{[CC]}n                        -> Polyethylene
{[CC(C1=CC=CC=C6)]}n          -> Polystyrene
{[CC]}<n:50-100>               -> Stochastic PE (50-100 units)
```

---

## Grammar Philosophy

SCRIPT is built on three principles:

1. **Paninian Grammar** (Sanskrit linguistics)
   - Sandhi: Context-aware bond validation
   - Lopa: Automatic hydrogen handling
   - Anuvritti: State continuity through the string

2. **Graph Linearization**
   - DFS traversal for canonicalization
   - Local ring encoding (no global state)
   - One molecule = one canonical string

3. **Invalid-Proof Design**
   - Generative state machine prevents valence violations
   - Every valid SCRIPT string maps to a valid molecule
   - Parse errors caught immediately, not after construction

---

## Comparison with SMILES

| Feature | SMILES | SCRIPT |
|---------|--------|--------|
| Ring notation | `C1CCCCC1` (global) | `C1CCCCC6` (local lookback) |
| Aromaticity | `c1ccccc1` (lowercase) | `C=CC=CC=C6` (Kekule) |
| Canonicalization | Multiple algorithms | One deterministic algorithm |
| Stereochemistry | Neighbor-order dependent | CIP-priority reconciled |
| Materials | Not supported | Alloys, surfaces, crystals |
| Polymers | Not standard | Built-in `{[unit]}n` |
| Validation | Post-hoc | Real-time (Sandhi rules) |

---

## Further Reading

- `SPEC.md` - Complete technical specification
- `CIP_STEREO_THEORY.md` - Stereochemistry reconciliation theory
- `STANDALONE_ARCHITECTURE.md` - Implementation architecture
- `examples/` - Usage examples and demos
- `benchmark.py` - 100-compound validation suite

---

*"A linear script to unfold molecular complexity."*
