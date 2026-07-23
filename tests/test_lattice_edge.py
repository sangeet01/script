#!/usr/bin/env python3
"""
V4 Lattice Extension — edge-case stress test.

Push the grammar/parser into weird corners and see what breaks. Useful for
finding the seams where the V4 extension interacts with V3 features.
"""

import sys
sys.path.insert(0, '/home/z/my-project/repos/script')

from script.parser import SCRIPTParser
from script.mol import BondType, CoreMolecule

p = SCRIPTParser()

results = []

def try_parse(name, s, expect_ok=True):
    r = p.parse(s)
    ok = r.get('success', False)
    mol = r.get('molecule') if ok else None
    if isinstance(mol, list):
        mol = next((m for m in mol if isinstance(m, CoreMolecule)), None)
    status = 'PASS' if (ok == expect_ok) else 'FAIL'
    results.append((status, name, s, ok, r.get('error', '') if not ok else ''))
    print(f'  [{status}] ok={ok}  {name:<40} {s!r}')
    if status == 'FAIL':
        err = r.get('error','') if not ok else '(parsed but expected to fail)'
        print(f'           err: {err[:120]}')
    return mol if ok else None

print('=== A. Spline bond combinations ===')
try_parse('spline + double bond',        'C~>C=C')
try_parse('spline + aromatic',           'C~>C:C')
try_parse('spline + dative',             'C~>N->B')
try_parse('spline ring closure',         'C~>C~>C&3-')
try_parse('spline branch',               'C~>(C~>C)C')
try_parse('spline in branch',            'C(C~>C)~>C')
try_parse('spline + periodic',           'C~>@1,0,0C')
try_parse('spline + beam radius',        'C<r:0.5>~>C<r:1.0>')
try_parse('multi-spline chain',          'C~>C~>C~>C~>C')
try_parse('spline + lattice tag',        '[Lattice:RandomSpline] C~>C~>C')

print('\n=== B. Beam radius edge cases ===')
try_parse('zero radius',                 'C<r:0>~>C')
try_parse('large radius',                'C<r:100.0>~>C')
try_parse('fractional radius',           'C<r:0.001>~>C')
try_parse('radius in bracket atom',      '[C<r:0.5>]~>C')   # may need state_block support inside bracket_atom
try_parse('radius with isotope',         'C<13,r:0.5>~>C')
try_parse('radius + charge',             'C<+1,r:0.5>~>C')
try_parse('radius + spin',               'C<s:3,r:0.5>~>C')

print('\n=== C. Lattice / thickness combinations ===')
try_parse('lattice only no thickness',   '[Lattice:BodyCentered] C')
try_parse('thickness only no lattice',   '[Thickness:Constant(2.0)] C')
try_parse('thickness then lattice',      '[Thickness:Constant(2.0)] [Lattice:BodyCentered] C')
try_parse('lattice then thickness',      '[Lattice:BodyCentered] [Thickness:Constant(2.0)] C')
try_parse('with polymer',                '[Lattice:BodyCentered] {[C]}<n:5>')
try_parse('with stochastic repeat',      '[Lattice:BodyCentered] {[C]}<n:50-100>')
try_parse('custom lattice type',         '[Lattice:MyCustomType] C')
try_parse('custom thickness class',      '[Thickness:MyCustom(1,2,3)] C')
try_parse('thickness no args',           '[Thickness:Boundary] C')

print('\n=== D. Namespace + bounding shape ===')
try_parse('geom namespace',              '[[geom:BaseSphere_r50]] C')
try_parse('xtal namespace',              '[[xtal:Rutile]] Ti(O)2')
try_parse('chem namespace',              '[[chem:Benzene]] C')
try_parse('bare BaseShape',              '[[BaseBox_20x20x20]] C')
try_parse('bare crystal',                '[[Rutile]] Ti(O)2')
try_parse('geom + lattice',              '[[geom:BaseSphere_r50]] [Lattice:BodyCentered] C')
try_parse('geom + lattice + thickness',  '[[geom:BaseSphere_r50]] [Thickness:Constant(2.0)] [Lattice:BodyCentered] C')
try_parse('geom + lattice + post',       '[[geom:BaseSphere_r50]] [Lattice:BodyCentered] C |> intersect')

print('\n=== E. Post-process pipeline ===')
try_parse('single no-arg',               'C |> smoothen')
try_parse('single with args',            'C |> overoffset(3,0)')
try_parse('chained 2 ops',               'C |> overoffset(3,0) |> intersect')
try_parse('chained 3 ops',               'C |> overoffset(3,0) |> smoothen |> intersect')
try_parse('kv args',                     'C |> render(subsample=5)')
try_parse('mixed kv + positional',       'C |> render(subsample=5,3)')  # may be ambiguous
try_parse('unknown op name',             'C |> custom_op(1,2)')
try_parse('no parens',                   'C |> intersect')
try_parse('empty parens',                'C |> intersect()')

print('\n=== F. Combined V3 + V4 (the killer tests) ===')
try_parse('lattice + periodic MOF bond', '[[geom:BaseBox_20x20x20]] [Lattice:BodyCentered] {[C-@1,0,0C-@0,1,0C-@0,0,1C]}<n:5>')
try_parse('lattice + alloy Vikalpa',     '[Lattice:BodyCentered] C<~0.5>~>C<~0.5>')
try_parse('lattice + spin state',        '[Lattice:BodyCentered] C<s:3>~>C<s:3>')
try_parse('lattice + stereochemistry',   '[Lattice:BodyCentered] [C@H](O)~>[C@@H](O)')
try_parse('lattice + charged',           '[Lattice:BodyCentered] [N+](=O)~>[O-]')
try_parse('lattice + peptide',           '[Lattice:Custom] {A.G.S.K}')  # may not be useful but should parse
try_parse('lattice + reaction',          '[Lattice:BodyCentered] C >> [Lattice:BodyCentered] CO')
try_parse('lattice + surface',           '[[geom:Pt_111]] | [Lattice:BodyCentered] >C=O')  # surface + lattice
try_parse('lattice + polymer block',     '[Lattice:BodyCentered] {[C]}<n:5> -b- {[C]}<n:5>')  # block copoly lattice

print('\n=== G. Should FAIL (negative tests) ===')
try_parse('unclosed lattice tag',        '[Lattice:BodyCentered C', expect_ok=False)
try_parse('lattice type lowercase',      '[lattice:bodycentered] C', expect_ok=False)
try_parse('thickness missing class',     '[Thickness:(2.0)] C', expect_ok=False)
try_parse('post-process no op',          'C |>', expect_ok=False)
try_parse('post-process bad op',         'C |> 123', expect_ok=False)

print('\n' + '='*60)
n_pass = sum(1 for r in results if r[0] == 'PASS')
n_fail = sum(1 for r in results if r[0] == 'FAIL')
print(f'  Edge case results: {n_pass} PASS, {n_fail} FAIL')
print('='*60)
