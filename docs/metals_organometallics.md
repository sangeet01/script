# SCRIPT Documentation: Metals & Organometallics

SCRIPT provides first-class support for coordination chemistry and organometallic complexes that are traditionally difficult to represent in linear notations.

## 1. Directional Coordination (Sandhi)
Traditional chemical formats treat all bonds as covalent shares. SCRIPT distinguishes between electron donation and sharing.

- **`->` (Dative Bond)**: Forward donation from L to M.
- **`<-` (Reverse Dative)**: Back-donation from M to L.
- **`>` (Coordinate Bond)**: Used specifically for surface adsorption or high-coordination centers where the bond is non-directional but distinct from a single covalent bond.

## 2. Haptic Bonds (The Star Bond `*n`)
Organometallics like Ferrocene involve a metal atom bonding to multiple atoms in a ring simultaneously (hapticity).

- **Syntax**: `*N` where `N` is the hapticity (eta-$n$).
- **Example (Ferrocene)**: `[Fe]*5&5:C5`
  *   `[Fe]` is the iron center.
  *   `*5` indicates an $\eta^5$ haptic bond.
  *   `&5:C5` is the cyclopentadienyl ring.

## 3. Bracket Atoms (Transition Metals)
Transitions metals and complex elements should be enclosed in brackets `[ ]` to bypass standard organic valence guards.

- **Example**: `[Pd]([P](Ph)3)2(Cl)2` (Bis(triphenylphosphine)palladium chloride).
- **Formal Charges**: `[Fe++]`, `[Pt+4]`.
- **Isotopes**: `[2H]`, `[13C]`.

## 4. Coordination Geometries
In the atom state block `< >`, you can specify the coordination geometry to disambiguate isomers that the topology alone cannot resolve.
- `<sqp>`: Square Planar
- `<oct>`: Octahedral
- `<tbp>`: Trigonal Bipyramidal
- `<tet>`: Tetrahedral
