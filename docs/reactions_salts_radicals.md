# SCRIPT Documentation: Reactions, Salts, and Radicals

## 1. Reaction Informatics (The Kriyā)
SCRIPT is designed for reaction handling as a first-class feature rather than a string modification.

- **Reaction Arrow (`>>`)**: Separator for reactants and products.
- **Atom Mapping (`:n`)**: Used to track atoms across the reaction arrow.
    - `[C:1]>>[C:1]`
    - Mapping is best placed inside brackets to avoid ambiguity with aromatic bonds.

## 2. Multi-Component Systems (Salts & Solvates)
Chemical structures often contain disconnected parts. SCRIPT uses two separators for isolation:

- **Dot Separator (`.`)**: Standard separation for disconnected components (e.g., salts).
    - `[Na+].[Cl-]`
- **Ionic/Salt Bridge (`~`)**: Represents an ionic interaction or a specific salt-bridge relationship.
    - `C~N` (Ionic interaction between components)

Each component in a salt or reaction is treated as an isolated **State Machine**, ensuring that valency calculations and ring-numberings do not leak across components.

## 3. Atomic States (State Blocks)
The `<...>` block within bracket atoms allows for high-density attribute tagging:

- **Isotopes**: `[C<13>]`
- **Charges**: `[N<+>]` or `[Ca<2+>]`
- **Tautomeric Mobility**: `[O<m>]` marks an atom as mobile in a transition state or resonant system.

## 4. Radicals
SCRIPT supports specific radical markers inside brackets:
- `.` : Monovalent radical (e.g., `[O.]`)
- `..` : Divalent diradical (e.g., `[C..]`)

## 5. Transition States
The **Star Bond (`*`)** represents resonant or partial bonds in a transition state.
Example: `C-H*O` indicates a hydrogen-bonding or transition state interaction between Hydrogen and Oxygen.
