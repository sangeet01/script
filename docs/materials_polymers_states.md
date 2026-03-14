# SCRIPT Documentation: Materials, Polymers, and States (V3)

V3 expands SCRIPT from discrete molecules into the realm of materials science, continuous phases, and quantum states.

## 1. Alloys & Occupancy (Vikalpa `~`)
For non-stoichiometric materials and doped semiconductors, SCRIPT uses the `~` symbol inside the atom state block to denote fractional occupancy.

- **Syntax**: `Atom<~FLOAT>`
- **Example**: `Ti<~0.9>N<~0.1>` (Titanium nitride doped with Nitrogen).

## 2. Macroscopic Context (Sthana `[[ ]]`)
Defines the bulk structure, space group, or crystalline phase in which the following entities exist.

- **Syntax**: `[[Label]]`
- **Example**: `[[Rutile]] Ti(O)2` vs `[[Anatase]] Ti(O)2`
- **Context Scope**: The context applies to all components following it until a phase boundary or string end.

## 3. Phase Boundaries (Adhikarana `|`)
Used for surface chemistry and interfaces where discrete molecules interact with a bulk phase.

- **Syntax**: `[[Surface]] | Molecule`
- **Example**: `[[Pt_111]] | >C=O` (Carbon Monoxide adsorbed onto a Platinum 111 surface).

## 4. Electronic & Excited States (Svara)
Quantum mechanical states are represented as "accents" on the atomic root.

- **Spin Multiplicity (`s:`)**: `<s:1>` (Singlet), `<s:3>` (Triplet).
- **Excited State (`*`)**: Marks a photo-excited or transition state.
- **Example**: `O=O<s:3>` (Ground state triplet oxygen).

## 5. Polymers & Repeating Units (Abhyasa)
Polymers use a distinct bracket notation to avoid ambiguity with peptide sequences.

- **Syntax**: `{[ Unit ]}Suffix`
- **Suffixes**:
  - `n`: Standard polymer repeat.
  - `<n:50-100>`: Stochastic range of units.
  - `-b-`: Block junction.
  - `-r-`: Random junction.
- **Example**: `{[CC]}n` (Polyethylene).

## 6. Reaction Mapping (`:n`)
Atom-to-atom mapping across reactions is preserved using the `:n` notation inside brackets.

- **Example**: `[C:1]>>[C:1]` (Preserves that the carbon in reactant is the same carbon in the product).
