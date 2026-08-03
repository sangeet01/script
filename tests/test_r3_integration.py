#!/usr/bin/env python3
"""
R3: End-to-end integration test for all V4/V4.1/V4.2 features.

Tests the complete pipeline:
  parse(S) → mol → canonicalize(mol) → S' → parse(S') → mol' → canonicalize(mol') → S''
  Verify S' == S'' (canonical form is stable under re-canonicalization)

Covers every feature added in Parts 1-4:
  - V4: spline ~>, beam radius <r:>, lattice tags, thickness tags, post-process |>
  - V4.1: generalized typed tags, generalized PIPE_OP, generalized namespaces
  - V4.2: float translations @0.5,0.5,0.5, explicit spline ~{x,y,z;...}
  - Sanskrit/Paninian: Sandhi validation, Lopa implicit Hs
  - CIP: @R/@S chirality, ring closures &N
  - Crystallography: [[xtal:SG]] namespace, periodic bonds
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer
from script.mol import BondType, CoreMolecule

p = SCRIPTParser()
c = SCRIPTCanonicalizer()

results = []

def check(name, s, expect_parse=True, expect_roundtrip=True):
    """End-to-end test: parse → canon → re-parse → re-canon → verify stable."""
    r = p.parse(s)
    if not r.get('success'):
        if expect_parse:
            results.append((name, False, f'parse failed: {r.get("error","")[:60]}'))
            print(f'  [FAIL] {name}: parse failed')
        else:
            results.append((name, True, 'expected parse failure'))
            print(f'  [PASS] {name}: correctly rejected')
        return

    if not expect_parse:
        results.append((name, False, 'expected parse failure but succeeded'))
        print(f'  [FAIL] {name}: should have been rejected')
        return

    mol = r['molecule']
    if isinstance(mol, list):
        mol = mol[0] if mol else None
    if mol is None:
        results.append((name, False, 'no molecule returned'))
        print(f'  [FAIL] {name}: no molecule')
        return

    canon1 = c.canonicalize_core(mol)
    if canon1 is None:
        results.append((name, False, 'canonicalize returned None'))
        print(f'  [FAIL] {name}: canon=None')
        return

    # Re-parse the canonical form
    r2 = p.parse(canon1)
    if not r2.get('success'):
        results.append((name, False, f're-parse failed: {r2.get("error","")[:60]}'))
        print(f'  [FAIL] {name}: re-parse failed')
        return

    mol2 = r2['molecule']
    if isinstance(mol2, list):
        mol2 = mol2[0] if mol2 else None
    canon2 = c.canonicalize_core(mol2)

    roundtrip_ok = (canon1 == canon2)
    if expect_roundtrip and roundtrip_ok:
        results.append((name, True, canon1))
        print(f'  [PASS] {name}')
        print(f'         canon: {canon1[:80]}')
    elif not expect_roundtrip:
        results.append((name, True, f'no roundtrip (expected): {canon1[:40]}'))
        print(f'  [PASS] {name}: parsed (no roundtrip expected)')
    else:
        results.append((name, False, f'{canon1!r} != {canon2!r}'))
        print(f'  [FAIL] {name}: roundtrip mismatch')
        print(f'         canon1: {canon1[:80]}')
        print(f'         canon2: {canon2[:80]}')


# ============================================================
print('='*72)
print('  R3: End-to-End Integration Test')
print('='*72)

# --- V1-V3 baseline (must still round-trip) ---
print('\n--- V1-V3 Baseline ---')
check('ethanol',        'CCO')
check('aspirin',        'C(C(=O)O):C(OC(=O)C):C:C:C:C&6:')
check('benzene',        'C:C:C:C:C:C&6:')
check('chiral glucose', 'O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO')

# --- V4 Lattice features ---
print('\n--- V4 Lattice ---')
check('spline bond',         'C~>C~>C')
check('spline + lattice',    '[Lattice:RandomSpline] C~>C~>C')
check('beam radius',         'C<r:0.5>C')
check('lattice tag',         '[Lattice:BodyCentered] C')
check('thickness tag',       '[Thickness:Constant(2.0)] C')
check('lattice + thickness', '[Thickness:Constant(2.0)] [Lattice:BodyCentered] C')
check('bounding + lattice',  '[[geom:BaseSphere_r50]] [Lattice:BodyCentered] C')
check('post-process',        'C |> smoothen')
check('post-process args',   'C |> overoffset(3,0) |> intersect')

# --- V4.1 Generalized ---
print('\n--- V4.1 Generalized ---')
check('typed tag Mesh',         '[Mesh:Icosphere(2)] C')
check('typed tag Material',     '[Material:Steel] C')
check('typed tag Solver kv',    '[Solver:FEM(order=2)] C')
check('typed tag Export str',   '[Export:STL(binary)] C')
check('multi typed tags',       '[Mesh:Icosphere(2)] [Material:Steel] C')
check('xtal namespace',         '[[xtal:Rutile]] Ti(O)2')
check('geom namespace',         '[[geom:BaseSphere_r50]] C')
check('custom namespace',       '[[bio:1UBQ]] C')

# --- V4.2 New ---
print('\n--- V4.2 New ---')
check('float translation',      'Na-@0.5,0.5,0.5Cl')
check('integer translation',    'C-@1,0,0C')
check('negative translation',   'C-@-0.5,0.5,-0.5C')
check('spline explicit',        'C~{0.0,0.0,0.0;1.0,1.0,1.0}C')
check('spline explicit neg',    'C~{-1.0,-2.0,-3.0;0.0,0.0,0.0}C')
check('crystal Fm-3m',          '[[xtal:Fm-3m]] Na-@0.5,0.5,0.5Cl')
check('crystal Im-3m',          '[[xtal:Im-3m]] Fe-@0.5,0.5,0.5Fe')
check('crystal Pm-3m',          '[[xtal:Pm-3m]] C-@1,0,0C')

# --- Combined V4 + V4.1 + V4.2 ---
print('\n--- Combined ---')
check('full lattice pipeline',
      '[[geom:BaseSphere_r50]] [Thickness:Constant(2.0)] [Lattice:BodyCentered] {[C]}<n:5> |> overoffset(3,0) |> intersect')
check('crystal + tags + pipeline',
      '[[xtal:Fm-3m]] [Wyckoff:Na_4a] [Wyckoff:Cl_4b] Na-@0.5,0.5,0.5Cl |> export(cif)')
check('gyroid lattice',
      '[[geom:Gyroid_k1.0]] [Lattice:Gyroid] C |> smoothen')

# --- Paninian features (Sandhi, Lopa, Sthiti) ---
print('\n--- Paninian Features ---')
check('Sandhi valid neopentane', 'C(C)(C)(C)C')        # 4-valent C (CH3)4 — valid
check('Sandhi valid 4-valent',  'C(C)(C)(C)')           # 4-valent C — valid
check('Sandhi charge adjust',   '[N+](=O)[O-]')         # charge adjusts valence
# Negative tests (should fail) — 6-valent carbon (5 substituents + would need 6th)
check('Sandhi invalid 6-valent', 'C(C)(C)(C)(C)C',      expect_parse=False)

# --- Summary ---
print('\n' + '='*72)
n_pass = sum(1 for _, ok, _ in results if ok)
n_fail = sum(1 for _, ok, _ in results if not ok)
print(f'  R3 Results: {n_pass} PASS, {n_fail} FAIL')
print('='*72)
