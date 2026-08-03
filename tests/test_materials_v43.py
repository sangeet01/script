#!/usr/bin/env python3
"""
Materials Science Test Suite — V4.3

Tests all materials science grammar extensions:
  T1.1: Defect atoms [_:Fe], [Df:C]C
  T1.2: Magnetic spin direction <m:0,0,2.2>
  T1.3: Property tags [Prop:Bandgap(1.2,eV)]
  T1.4: K-point paths [KPath:Gamma-X-M-Gamma]
  T1.5: Composition ranges <~0.5-0.7>
  T2.1: Phase transitions >>(T=1185)
  T2.3: Modulation q-vectors [[xtal:P1;q:0.31,0,0]]
  T3:   Topological, 2D, grain boundaries, etc.

Run: python /home/z/my-project/scripts/test_materials_v43.py
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer

p = SCRIPTParser()
c = SCRIPTCanonicalizer()

results = []

def check(name, s, expect_parse=True):
    r = p.parse(s)
    ok = r.get('success', False) == expect_parse
    results.append((name, ok))
    status = 'PASS' if ok else 'FAIL'
    print(f'  [{status}] {name:<45} {s[:55]!r}')
    if not ok and r.get('success'):
        print(f'         Expected fail but parsed')
    elif not ok and not r.get('success'):
        print(f'         err: {r.get("error","")[:60]}')


print('='*72)
print('  V4.3 Materials Science Test Suite')
print('='*72)

# ============================================================
print('\n--- T1.5: Composition Ranges (HEA, solid solutions) ---')
check('binary composition range',    'Fe<~0.5-0.7>Cr<~0.3-0.5>')
check('ternary HEA range',           'Fe<~0.4-0.6>Cr<~0.2-0.3>Ni<~0.2-0.3>')
check('Ti-6Al-4V range',             'Ti<~0.8-0.9>Al<~0.05-0.07>V<~0.03-0.05>')
check('point occupancy still works', 'Ti<~0.9>N<~0.1>')

# ============================================================
print('\n--- T1.1: Defect Atoms (vacancies, interstitials) ---')
check('Fe vacancy',                  '[_:Fe]')
check('Cl vacancy',                  '[_:Cl]')
check('C interstitial',              '[Df:C]C')
check('F-center (vacancy + e-)',     '[_:Cl<->]')
check('interstitial with occupancy', '[Df:C]C<occ:0.5>')
check('vacancy in lattice',          '[[xtal:Fm-3m]] Na-@0.5,0.5,0.5[_:Cl]')

# ============================================================
print('\n--- T1.2: Magnetic Ordering ---')
check('Fe magnetic moment',          'Fe<m:0,0,2.2>')
check('antiparallel spin',           'Fe<m:0,0,-2.2>')
check('AFM pair',                    'Fe<m:0,0,2.2>~>Fe<m:0,0,-2.2>')
check('zero moment (paramagnetic)',  'Fe<m:0,0,0>')
check('magnetic space group',        '[[mag:Pnma]] Fe<m:0,0,2.2>')
check('spin multiplicity still works', 'O=O<s:3>')

# ============================================================
print('\n--- T1.3: Property Tags ---')
check('bandgap direct',              '[Prop:Bandgap(1.2,eV)] C')
check('bandgap indirect',            '[Prop:Bandgap(1.2,eV,indirect)] C')
check('formation energy',            '[Prop:Eform(-3.5,eV_atom)] C')
check('superconducting Tc',          '[Prop:Tc(9.2,K,superconducting)] C')
check('magnetization',               '[Prop:Magnetization(2.2,muB)] Fe')
check('hardness',                    '[Prop:Hardness(8.5,GPa)] C')
check('multiple properties',         '[Prop:Bandgap(1.2,eV)] [Prop:Eform(-3.5,eV_atom)] C')

# ============================================================
print('\n--- T1.4: K-point Paths ---')
check('Gamma-X-M-Gamma path',        '[KPath:Gamma-X-M-Gamma] C')
check('hexagonal BZ path',           '[KPath:Gamma-M-K-Gamma] C')

# ============================================================
print('\n--- T2.1: Phase Transitions ---')
check('alpha-gamma Fe transition',   '[[xtal:BCC]] Fe >>(T=1185) [[xtal:FCC]] Fe')
check('BaTiO3 ferroelectric',        '[[xtal:Tetragonal]] BaTiO3 >>(T=393) [[xtal:Cubic]] BaTiO3')
check('pressure-induced transition', '[[xtal:Graphite]] C >>(P=15) [[xtal:Diamond]] C')
check('multi-condition transition',  '[[xtal:BCC]] Fe >>(T=1185,P=1) [[xtal:FCC]] Fe')
check('simple reaction still works', 'CC >> CCO')

# ============================================================
print('\n--- T2.3: Modulated Structures ---')
check('incommensurate q-vector',     '[[xtal:P1;q:0.31,0,0]] C')
check('commensurate q',              '[[xtal:P1;q:0.5,0,0]] C')
check('q with lattice params',       '[[xtal:P1;5,5,5,90,90,90;q:0.31,0,0]] C')

# ============================================================
print('\n--- T3: Specialized ---')
check('topological insulator Z2',    '[Topo:Z2(1)] Bi2Se3')
check('Chern insulator',             '[Topo:Chern(2)] C')
check('Weyl semimetal',              '[Topo:Weyl(1)] C')
check('graphene monolayer',          '[Layer:Graphene(1)] C')
check('Bernal bilayer',              '[Stack:AB(2,Graphene)] C')
check('twisted bilayer',             '[Stack:Twist(1.1)] C')
check('grain boundary Sigma5',       '[Grain:Sigma5(310)] Fe')
check('twin boundary Sigma3',        '[Grain:Sigma3(111)] Cu')
check('phonon mode',                 '[Phonon:Gamma_5(8.5,THz)] C')
check('elastic tensor',              '[Tensor:Elastic(C11=243,C12=146,GPa)] Fe')
check('liquid crystal nematic',      '[Mesophase:Nematic] C')
check('battery intercalation',       '[Site:Intercalation(0.5)] Li')
check('OER active site',             '[Site:Active(OER)] C')
check('AFM ordering type',           '[Magnet:Atype] Fe')
check('non-collinear magnetic',      '[Magnet:NonCollinear] Fe')
check('magnetic space group',        '[[msg:Pnma]] Fe')

# ============================================================
print('\n--- Combined: Real Materials ---')
check('doped Si (semiconductor)',    '[[xtal:Fd-3m]] Si-@0.25,0.25,0.25[_:Si<occ:0.001>][Df:P<occ:0.001>]')
check('NaCl with properties',        '[[xtal:Fm-3m]] [Prop:Bandgap(7,eV)] Na-@0.5,0.5,0.5Cl')
check('magnetic BCC Fe',             '[[xtal:Im-3m]] [Prop:Magnetization(2.2,muB)] Fe<m:0,0,2.2>')
check('YBCO superconductor',         '[[xtal:Pmmm]] [Prop:Tc(92,K,superconducting)] YBa2Cu3O7')
check('graphene with twist',         '[[mat:Graphene]] [Stack:Twist(1.1)] [Prop:Bandgap(0.01,eV)] C')

# ============================================================
print('\n' + '='*72)
n_pass = sum(1 for _, ok in results if ok)
n_fail = sum(1 for _, ok in results if not ok)
print(f'  V4.3 Materials Science: {n_pass} PASS, {n_fail} FAIL')
print('='*72)
