#!/usr/bin/env python3
"""
V4 Lattice Extension — test suite for the LEAP71 bridge grammar.

Tests that the new grammar + IR + parser all wire together correctly
for the four canonical LEAP71 worked examples:

  1. Regular body-centered lattice in a sphere  (RegularTask)
  2. Conformal body-centered lattice in a pipe segment
  3. Random-spline lattice on a noised regular grid
  4. Post-processed lattice infill with OverOffset + Intersect

Plus edge cases: namespaces, beam radius, no-thickness fallback, etc.

Run: python /home/z/my-project/scripts/test_lattice_v4.py
"""

import sys
sys.path.insert(0, '/home/z/my-project/repos/script')

from script.parser import SCRIPTParser
from script.mol import BondType, CoreMolecule

p = SCRIPTParser()

# ---------- helpers ----------
def parse(s):
    r = p.parse(s)
    if not r.get('success'):
        raise AssertionError(f"Parse failed: {s!r}\n  err: {r.get('error','')[:200]}")
    mol = r.get('molecule')
    if isinstance(mol, list):
        # Salts / multi-component — pick the first that is a CoreMolecule
        mol = next((m for m in mol if isinstance(m, CoreMolecule)), mol[0])
    return mol

def expect(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f'  [OK] {msg}')

# ---------- Test 1: spline bond ~>
print('\n=== Test 1: spline bond (~>) ===')
mol = parse('C~>C')
expect(len(mol.atoms) == 2, 'spline bond has 2 atoms')
expect(len(mol.bonds) == 1, 'spline bond has 1 bond')
expect(mol.bonds[0].bond_class == 'spline', f"bond_class is 'spline' (got {mol.bonds[0].bond_class!r})")

# Multi-spline
mol = parse('C~>C~>C~>C')
expect(len(mol.atoms) == 4, 'multi-spline has 4 atoms')
expect(len(mol.bonds) == 3, 'multi-spline has 3 bonds')
expect(all(b.bond_class == 'spline' for b in mol.bonds), 'all bonds are spline')

# ---------- Test 2: beam radius state attr <r:..> ----------
print('\n=== Test 2: beam radius <r:..> ===')
mol = parse('C<r:0.5>C')
expect(mol.atoms[0].beam_radius == 0.5, f"atom[0].beam_radius=0.5 (got {mol.atoms[0].beam_radius})")
expect(mol.atoms[1].beam_radius is None, f"atom[1].beam_radius=None (got {mol.atoms[1].beam_radius})")

# Multiple atoms with beam radius
mol = parse('C<r:1.0>C<r:2.0>C<r:3.0>')
expect(all(mol.atoms[i].beam_radius == float(i+1) for i in range(3)),
       f"all three atoms have beam_radius 1.0, 2.0, 3.0 (got {[a.beam_radius for a in mol.atoms]})")

# ---------- Test 3: lattice_type tag ----------
print('\n=== Test 3: lattice_type tag ===')
for ltype in ('BodyCentered', 'Octahedron', 'RandomSpline', 'Custom'):
    mol = parse(f'[Lattice:{ltype}] C')
    expect(mol.lattice_type == ltype,
           f"[Lattice:{ltype}] sets lattice_type={ltype!r} (got {mol.lattice_type!r})")

# Custom-named lattice type (regex fallback)
mol = parse('[Lattice:MyWeirdLattice] C')
expect(mol.lattice_type == 'MyWeirdLattice',
       f"custom lattice type 'MyWeirdLattice' (got {mol.lattice_type!r})")

# ---------- Test 4: thickness tag with args ----------
print('\n=== Test 4: thickness tag with args ===')
mol = parse('[Thickness:Constant(2.0)] C')
expect(mol.thickness_class == 'Constant', f"thickness_class=Constant (got {mol.thickness_class})")
expect(mol.thickness_args == (2.0,), f"thickness_args=(2.0,) (got {mol.thickness_args})")

mol = parse('[Thickness:CellBased(1,4)] C')
expect(mol.thickness_class == 'CellBased', f"thickness_class=CellBased")
expect(mol.thickness_args == (1, 4), f"thickness_args=(1,4) (got {mol.thickness_args})")

mol = parse('[Thickness:Boundary] C')
expect(mol.thickness_class == 'Boundary', f"thickness_class=Boundary")
expect(mol.thickness_args == (), f"thickness_args=() for no-arg class (got {mol.thickness_args})")

mol = parse('[Thickness:GlobalFunc(x)] C')
expect(mol.thickness_class == 'GlobalFunc', f"thickness_class=GlobalFunc")
expect(mol.thickness_args == ('x',), f"thickness_args=('x',) (got {mol.thickness_args})")

# ---------- Test 5: macroscopic context + namespace ----------
print('\n=== Test 5: macroscopic context + namespace ===')
mol = parse('[[BaseSphere_r50]] C')
expect(mol.macroscopic_context == 'BaseSphere_r50', f"context=BaseSphere_r50 (got {mol.macroscopic_context})")
expect(mol.context_namespace is None, f"namespace=None for bare (got {mol.context_namespace})")

mol = parse('[[geom:BaseSphere_r50]] C')
expect(mol.context_namespace == 'geom', f"namespace=geom (got {mol.context_namespace})")

mol = parse('[[xtal:Rutile]] Ti(O)2')
expect(mol.context_namespace == 'xtal', f"namespace=xtal (got {mol.context_namespace})")

# Backward compat: bare Rutile still works
mol = parse('[[Rutile]] Ti(O)2')
expect(mol.macroscopic_context == 'Rutile', f"bare Rutile backward-compat (got {mol.macroscopic_context})")

# ---------- Test 6: post-process pipeline ----------
print('\n=== Test 6: post-process pipeline ===')
mol = parse('C |> intersect')
expect(mol.post_process_ops == [('intersect', ())],
       f"single no-arg post-process (got {mol.post_process_ops})")

mol = parse('C |> overoffset(3,0)')
expect(mol.post_process_ops == [('overoffset', (3, 0))],
       f"single positional post-process (got {mol.post_process_ops})")

mol = parse('C |> overoffset(3,0) |> intersect')
expect(mol.post_process_ops == [('overoffset', (3, 0)), ('intersect', ())],
       f"chained post-process preserves order (got {mol.post_process_ops})")

mol = parse('C |> render(subsample=5)')
expect(mol.post_process_ops == [('render', (('subsample', 5),))],
       f"kv post-process (got {mol.post_process_ops})")

# ---------- Test 7: Worked Example 1 — RegularTask ----------
print('\n=== Test 7: Worked Example 1 — RegularTask ===')
s = '[[BaseSphere_r50]] [Thickness:Constant(2.0)] [Lattice:BodyCentered] {[C]}<n:5>'
mol = parse(s)
expect(mol.macroscopic_context == 'BaseSphere_r50', f"bounding=BaseSphere_r50")
expect(mol.thickness_class == 'Constant', f"thickness=Constant")
expect(mol.thickness_args == (2.0,), f"thickness args=(2.0,)")
expect(mol.lattice_type == 'BodyCentered', f"lattice=BodyCentered")
expect(mol.repeat_count == 5, f"repeat_count=5 (got {mol.repeat_count})")

# ---------- Test 8: Worked Example 2 — ConformalTask ----------
print('\n=== Test 8: Worked Example 2 — ConformalTask ===')
s = '[[BasePipeSegment_r10_h20_a45]] [Thickness:Boundary] [Lattice:BodyCentered] {[C]}<n:8> |> intersect'
mol = parse(s)
expect(mol.macroscopic_context == 'BasePipeSegment_r10_h20_a45', f"conformal bounding")
expect(mol.thickness_class == 'Boundary', f"boundary thickness")
expect(mol.lattice_type == 'BodyCentered', f"lattice=BodyCentered")
expect(mol.post_process_ops == [('intersect', ())], f"post-process=intersect")

# ---------- Test 9: Worked Example 3 — RandomSpline + noised ----------
print('\n=== Test 9: Worked Example 3 — RandomSpline + noised ===')
s = '[[BaseBox_100x100x100]] [Thickness:CellBased(1,4)] [Lattice:RandomSpline] {[C~>[*]~>[*]~>C]}<n:8>'
mol = parse(s)
expect(mol.lattice_type == 'RandomSpline', f"lattice=RandomSpline")
expect(mol.thickness_class == 'CellBased', f"thickness=CellBased")
expect(mol.thickness_args == (1, 4), f"thickness args=(1,4)")
expect(any(b.bond_class == 'spline' for b in mol.bonds), f"contains spline bonds")
expect(mol.macroscopic_context == 'BaseBox_100x100x100', f"box bounding")

# ---------- Test 10: Worked Example 4 — post-processed infill ----------
print('\n=== Test 10: Worked Example 4 — post-processed infill ===')
s = '[[BaseSphere_r50]] [Lattice:BodyCentered] [Thickness:Constant(2.0)] {[C]}<n:5> |> overoffset(3,0) |> intersect'
mol = parse(s)
expect(len(mol.post_process_ops) == 2, f"2 post-process ops (got {len(mol.post_process_ops)})")
expect(mol.post_process_ops[0] == ('overoffset', (3, 0)), f"first op=overoffset(3,0)")
expect(mol.post_process_ops[1] == ('intersect', ()), f"second op=intersect")

# ---------- Test 11: Periodic bonds still work ----------
print('\n=== Test 11: periodic bonds still work ===')
mol = parse('Ti->@0,0,1O')
expect(mol.bonds[0].translation == (0, 0, 1), f"translation=(0,0,1) (got {mol.bonds[0].translation})")
expect(mol.is_periodic == True, f"is_periodic=True")

# ---------- Test 12: Combined lattice + periodic bond ----------
print('\n=== Test 12: lattice with periodic inter-cell bond ===')
s = '[[BaseBox_20x20x20]] [Lattice:BodyCentered] {[C-@1,0,0C-@0,1,0C-@0,0,1C]}<n:5>'
mol = parse(s)
expect(mol.lattice_type == 'BodyCentered', f"lattice=BodyCentered")
expect(mol.is_periodic == True, f"is_periodic=True (got {mol.is_periodic})")
expect(any(b.translation != (0,0,0) for b in mol.bonds), f"has at least one cross-cell bond")

# ---------- Test 13: mixed lattice cell + beam radius ----------
print('\n=== Test 13: lattice cell with per-vertex beam radius ===')
s = '[Lattice:BodyCentered] C<r:1.0>~>C<r:2.0>'
mol = parse(s)
expect(mol.lattice_type == 'BodyCentered', f"lattice=BodyCentered")
expect(mol.atoms[0].beam_radius == 1.0, f"atom[0] radius=1.0")
expect(mol.atoms[1].beam_radius == 2.0, f"atom[1] radius=2.0")
expect(any(b.bond_class == 'spline' for b in mol.bonds), f"has spline bond")

print('\n' + '='*60)
print('  ALL V4 LATTICE TESTS PASSED')
print('='*60)
