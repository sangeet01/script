#!/usr/bin/env python3
"""
Comprehensive Chemistry Audit — test SCRIPT against every domain of chemistry.

Goal: find every molecule type that fails to parse, so we can fix the gaps.
"""

import sys
sys.path.insert(0, '/home/z/my-project/repos/script')

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer

p = SCRIPTParser()
c = SCRIPTCanonicalizer()

results = []

def test(category, name, script_str, expected_atoms=None):
    r = p.parse(script_str)
    ok = r.get('success', False)
    mol = r.get('molecule') if ok else None
    if isinstance(mol, list): mol = mol[0] if mol else None
    
    n_a = len(mol.atoms) if mol else 0
    n_b = len(mol.bonds) if mol else 0
    
    a_ok = expected_atoms is None or n_a == expected_atoms
    status = 'PASS' if (ok and a_ok) else 'FAIL'
    results.append((category, name, status))
    
    if status == 'PASS':
        print(f'  [PASS] {category:<25} {name:<35} atoms={n_a}')
    else:
        err = r.get('error','')[:70] if not ok else f'atoms={n_a} (expected {expected_atoms})'
        print(f'  [FAIL] {category:<25} {name:<35} {script_str[:40]}')
        print(f'         err: {err}')

print('='*80)
print('  COMPREHENSIVE CHEMISTRY AUDIT')
print('='*80)

# ============================================================
print('\n--- 1. Basic Organic Molecules ---')
test('organic', 'methane', 'C')
test('organic', 'ethane', 'CC')
test('organic', 'ethanol', 'CCO')
test('organic', 'acetic acid', 'CC(=O)O')
test('organic', 'aspirin (SCRIPT)', 'CC(=O)OC:C:C:C:C:C&6:C(=O)O')
test('organic', 'glucose', 'O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO')

# ============================================================
print('\n--- 2. Aromatic Compounds ---')
test('aromatic', 'benzene', 'C:C:C:C:C:C&6:')
test('aromatic', 'naphthalene', 'C:C:C:C:C&6:C:C:C:C&6:')
test('aromatic', 'pyridine', 'C:C:C:N:C:C&6:')
test('aromatic', 'pyrrole', 'C:C:C:C:N&5:')
test('aromatic', 'furan', 'C:C:C:C:O&5:')

# ============================================================
print('\n--- 3. Ionic Compounds ---')
test('ionic', 'NaCl', '[Na+].[Cl-]')
test('ionic', 'NaCl (ionic pair)', '[Na+]~[Cl-]')
test('ionic', 'MgSO4 (expanded)', '[Mg+2].[S]([O-])([O-])(=O)(=O)')
test('ionic', 'CaCO3 (expanded)', '[Ca+2].[C]([O-])([O-])(=O)')
test('ionic', 'ammonium chloride', '[NH4+].[Cl-]')

# ============================================================
print('\n--- 4. Coordination Complexes ---')
test('coord', 'hexaaquairon', '[Fe](O)(O)(O)(O)(O)(O)')
test('coord', 'cisplatin', '[Pt]([Cl])([Cl])([NH3])([NH3])')
test('coord', 'ferrocene', '[Fe]*5C:C:C:C:C&5:*5C:C:C:C:C&5:')
test('coord', 'EDTA-Mg (charged N)', '[Mg]([N+](C(=O)O)(C(=O)O)CC[N+](C(=O)O)C(=O)O)')
test('coord', 'heme (simplified)', '[Fe]*4N:C:C:C&6:C:C:C&6:*4N')

# ============================================================
print('\n--- 5. Radicals ---')
test('radical', 'methyl radical', '[C.]')
test('radical', 'hydroxyl radical', '[O.]')
test('radical', 'oxygen triplet', 'O=O<s:3>')
test('radical', 'nitric oxide', 'N=O<s:2>')
test('radical', 'divalent carbene', '[C..]')

# ============================================================
print('\n--- 6. Electron-Deficient (3c2e) ---')
test('electron-def', 'diborane', 'B(H)(H)<>H<>B(H)(H)H', expected_atoms=8)
test('electron-def', 'Al2Me6 dimer', '[Al](C)(C)<>C<>[Al](C)(C)', expected_atoms=7)
test('electron-def', 'BeH2 dimer', 'Be(H)<>H<>Be(H)', expected_atoms=5)

# ============================================================
print('\n--- 7. Hypercoordinate (hypervalent) ---')
test('hyper', 'SF6', '[S](F)(F)(F)(F)(F)(F)')
test('hyper', 'PF5', '[P](F)(F)(F)(F)(F)')
test('hyper', 'XeF2', '[Xe](F)(F)')
test('hyper', 'XeF4', '[Xe](F)(F)(F)(F)')
test('hyper', 'ClF3', '[Cl](F)(F)(F)')
test('hyper', 'IF7', '[I](F)(F)(F)(F)(F)(F)(F)')
test('hyper', 'sulfate SO4', '[S](=O)(=O)(O)(O)')

# ============================================================
print('\n--- 8. Isotopes ---')
test('isotope', 'deuterium', '[2H]')
test('isotope', 'tritium', '[3H]')
test('isotope', 'heavy water D2O', '[2H]O[2H]')
test('isotope', 'carbon-14', '[14C]')
test('isotope', 'carbon-13 methane', '[13C]')
test('isotope', 'iodine-131', '[131I]')
test('isotope', 'uranium-235', '[235U]')

# ============================================================
print('\n--- 9. Organometallics ---')
test('organometallic', 'Grignard', 'C[Mg]')
test('organometallic', 'MeLi tetramer (simplified)', 'C[Li]')
test('organometallic', 'Zeise salt', '[Pt]([Cl])([Cl])(C=C)')
test('organometallic', 'Wilkinson catalyst', '[Rh]([Cl])([P])([P])([P])')
test('organometallic', 'dative BN adduct', 'N->[B](F)(F)F')

# ============================================================
print('\n--- 10. Polymers ---')
test('polymer', 'polyethylene', '{[CC]}<n:50>')
test('polymer', 'polypropylene', '{[C(C)C]}<n:50>')
test('polymer', 'nylon-6,6', '{[CCCCCCC(=O)NCCCCCCN]}<n:50>')
test('polymer', 'KEVLAR (SCRIPT)', '{[C(=O)C:C:C:C(:C&6:)N]}<n:50>')
test('polymer', 'block copolymer', '{[CC]}<n:50> -b- {[CCCO]}<n:100>')

# ============================================================
print('\n--- 11. Biomolecules ---')
test('bio', 'glycine peptide', '{G}')
test('bio', 'alanine peptide', '{A}')
test('bio', 'tripeptide', '{A.G.S}')
test('bio', 'DNA strand', '{dA.dG.dC.dT}')
test('bio', 'RNA strand', '{rA.rG.rC.rU}')
test('bio', 'methylated DNA', '{m5C.dG}')
test('bio', 'phosphoserine', '{pS}')

# ============================================================
print('\n--- 12. Complex Stereochemistry ---')
test('stereo', 'R-lactic acid', '[C@H](O)(C)C(=O)O')
test('stereo', 'S-lactic acid', '[C@@H](O)(C)C(=O)O')
test('stereo', 'CIP-absolute R', '[C@RH](O)(C)C(=O)O')
test('stereo', 'allene', 'C=C=C')
test('stereo', 'biphenyl (SCRIPT)', 'C:C:C:C:C:C&6:-C:C:C:C:C:C&6:')
test('stereo', 'helicene (simplified)', 'C:C:C:C&6:C:C:C:C&6:C:C:C:C&6:')

# ============================================================
print('\n--- 13. Materials / Crystals ---')
test('materials', 'NaCl crystal', '[[xtal:Fm-3m]] Na-@0.5,0.5,0.5Cl')
test('materials', 'BCC iron', '[[xtal:Im-3m]] Fe-@0.5,0.5,0.5Fe')
test('materials', 'diamond', '[[xtal:Fd-3m]] C-@0.125,0.125,0.125C')
test('materials', 'TiO2 rutile', '[[xtal:Rutile]] Ti(O)2')
test('materials', 'steel alloy', 'Fe<~0.98>Cr<~0.02>')
test('materials', 'doped Si', '[[xtal:Fd-3m]] [Si<occ:0.999>][Df:P<occ:0.001>]')

# ============================================================
print('\n--- 14. Reactions ---')
test('reaction', 'combustion', 'CCO >> C(=O)O')
test('reaction', 'esterification', 'CC(=O)O.CO >> CC(=O)OC')
test('reaction', 'with catalyst', 'CC >[Pd]> CC=CC')
test('reaction', 'reversible', 'CC=O <=> CC(O)')

# ============================================================
print('\n--- 15. Edge Cases ---')
test('edge', 'single H atom', '[H]')
test('edge', 'single He', '[He]')
test('edge', 'hydrogen molecule', 'HH')
test('edge', 'noble gas', '[Ne]')
test('edge', 'transition metal', '[Fe]')
test('edge', 'lanthanide', '[Eu]')
test('edge', 'actinide', '[U]')
test('edge', 'transuranium', '[Pu]')

# ============================================================
print('\n--- 16. Surfaces and Interfaces ---')
test('surface', 'CO on Pt', '[[Pt_111]] | C=O')
test('surface', 'water on Au', '[[Au_111]] | O')

# ============================================================
print('\n--- 17. Mixed/Complex ---')
test('complex', 'methane hydrate (simplified)', 'C.C.C.C.C.C&6:')
test('complex', 'buckyball C60 (simplified)', 'C:C:C:C&6:C:C:C:C&6:C:C:C:C&6:C:C:C:C&6:C:C:C:C&6:')
test('complex', 'graphene fragment', 'C:C:C:C:C:C&6:C:C:C:C&6:')
test('complex', 'nanotube (simplified)', '{[C:C:C:C:C:C&6:]}<n:10>')

# ============================================================
print('\n--- 18. Additional V4.3 Features ---')
test('v4.3', 'defect vacancy', '[_:Fe]C')
test('v4.3', 'defect interstitial', '[Df:C]C')
test('v4.3', 'magnetic moment', 'Fe<m:0,0,2.2>')
test('v4.3', 'composition range', 'Fe<~0.5-0.7>Cr<~0.3-0.5>')
test('v4.3', 'phase transition', '[[xtal:BCC]] Fe >>(T=1185) [[xtal:FCC]] Fe')
test('v4.3', 'recursive SMARTS', '[?(C(=O)O)]C')
test('v4.3', 'bridge bond diborane', 'B(H)(H)<>H<>B(H)(H)H')
test('v4.3', 'spline bond', 'C~>C~>C')
test('v4.3', 'property tag', '[Prop:Bandgap(1.2,eV)] C')
test('v4.3', 'k-path', '[KPath:Gamma-X-M-Gamma] C')

# ============================================================
print('\n' + '='*80)
n_pass = sum(1 for _, _, s in results if s == 'PASS')
n_fail = sum(1 for _, _, s in results if s == 'FAIL')
print(f'  AUDIT RESULTS: {n_pass} PASS, {n_fail} FAIL out of {len(results)}')
print('='*80)

# List all failures
if n_fail > 0:
    print('\n  FAILURES BY CATEGORY:')
    by_cat = {}
    for cat, name, status in results:
        if status == 'FAIL':
            by_cat.setdefault(cat, []).append(name)
    for cat, names in by_cat.items():
        print(f'    {cat}: {", ".join(names)}')
