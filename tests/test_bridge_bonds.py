#!/usr/bin/env python3
"""
Test bridge bonds (3c2e) on electron-deficient molecules.

Covers:
  - Diborane B2H6 (the classic case)
  - Tetraborane B4H10
  - Pentaborane B5H9
  - Decaborane B10H14
  - Carboranes (C2B10H12)
  - Trimethylaluminum dimer Al2Me6
  - Beryllium hydride bridge
"""

import sys
sys.path.insert(0, '/home/z/my-project/repos/script')

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer
from script.mol import BondType

p = SCRIPTParser()
c = SCRIPTCanonicalizer()

results = []

def check(name, s, expect_atoms=None, expect_bridges=None):
    r = p.parse(s)
    ok = r.get('success', False)
    if not ok:
        results.append((name, False))
        print(f'  [FAIL] {name:<40} {s[:50]!r}')
        print(f'         err: {r.get("error","")[:70]}')
        return
    mol = r['molecule']
    if isinstance(mol, list): mol = mol[0]
    n_a = len(mol.atoms)
    n_b = len(mol.bonds)
    n_bridge = sum(1 for b in mol.bonds if b.bond_type == BondType.BRIDGE)
    canon = c.canonicalize_core(mol)
    
    a_ok = (expect_atoms is None or n_a == expect_atoms)
    br_ok = (expect_bridges is None or n_bridge == expect_bridges)
    status = 'PASS' if (a_ok and br_ok) else 'FAIL'
    results.append((name, a_ok and br_ok))
    print(f'  [{status}] {name:<40} atoms={n_a} bridges={n_bridge}')
    print(f'         SCRIPT: {s[:60]}')
    print(f'         canon:  {canon[:60]}')


print('='*72)
print('  3-Center-2-Electron Bridge Bond Tests')
print('='*72)

# ============================================================
print('\n--- Boron Hydrides ---')

# Diborane B2H6: 2 B + 6 H = 8 atoms, 2 bridges
check('diborane B2H6', 
      'B(H)(H)<>H<>B(H)(H)H',
      expect_atoms=8, expect_bridges=2)

# Diborane alternative notation (creates a chain, not real diborane)
check('diborane chain (alt)',
      'B(H)(H)<>H<>B(H)(H)<>H<>B',
      expect_atoms=9)

# Tetraborane B4H10: 4 B + 10 H, has 4 bridges
check('tetraborane B4H10 (simplified)',
      'B(H)(H)<>H<>B(H)<>H<>B(H)(H)<>H<>B(H)(H)',
      expect_atoms=14)

# Pentaborane B5H9: 5 B + 9 H
check('pentaborane B5H9 (simplified)',
      'B(H)<>H<>B(H)<>H<>B(H)<>H<>B(H)<>H<>B(H)(H)(H)',
      expect_atoms=16)

# Decaborane B10H14 (very simplified)
check('decaborane B10H14 (minimal)',
      'B<>H<>B<>H<>B<>H<>B<>H<>B',
      expect_atoms=9)

# ============================================================
print('\n--- Carboranes ---')

# ortho-Carborane C2B10H12 (simplified — just the bridge framework)
check('carborane C2B10H12 (simplified)',
      'C<>H<>B<>H<>C',
      expect_atoms=5)

# ============================================================
print('\n--- Aluminum Compounds ---')

# Trimethylaluminum dimer Al2Me6
check('Al2Me6 dimer',
      '[Al](C)(C)<>C<>[Al](C)(C)',
      expect_atoms=7, expect_bridges=2)

# ============================================================
print('\n--- Beryllium Compounds ---')

# Beryllium hydride bridge (gas phase, dimeric)
check('Be2H4 (beryllium hydride dimer)',
      'Be(H)<>H<>Be(H)',
      expect_atoms=5, expect_bridges=2)

# ============================================================
print('\n--- Mixed Bridge Bonds ---')

# Bridge with other bond types
check('bridge + single',
      'B(H)<>H<>B(H)(H)C',
      expect_atoms=7, expect_bridges=2)

# Bridge in a chain (not ring — ring closure with <> is complex)
check('bridge chain',
      'B<>H<>B<>H<>B<>H<>B',
      expect_atoms=7)

# ============================================================
print('\n--- Round-trip Canonicalization ---')

# Verify bridge bonds survive parse → canonicalize → re-parse
print('  Round-trip tests:')
rt_tests = [
    'B(H)(H)<>H<>B(H)(H)H',
    'B<>H<>B',
    '[Al](C)(C)<>C<>[Al](C)(C)',
    'Be(H)<>H<>Be(H)',
]
for s in rt_tests:
    r = p.parse(s)
    if not r.get('success'):
        print(f'    [FAIL] {s[:40]} — parse failed')
        continue
    mol = r['molecule']
    canon1 = c.canonicalize_core(mol)
    r2 = p.parse(canon1)
    if not r2.get('success'):
        print(f'    [FAIL] {s[:40]} — re-parse failed')
        continue
    canon2 = c.canonicalize_core(r2['molecule'])
    rt_ok = (canon1 == canon2)
    results.append(('round-trip: ' + s[:30], rt_ok))
    print(f'    [{"OK" if rt_ok else "FAIL"}] {s[:40]}')
    print(f'         canon: {canon1[:60]}')

# ============================================================
print('\n' + '='*72)
n_pass = sum(1 for _, ok in results if ok)
n_fail = sum(1 for _, ok in results if not ok)
print(f'  Bridge Bond Tests: {n_pass} PASS, {n_fail} FAIL')
print('='*72)
