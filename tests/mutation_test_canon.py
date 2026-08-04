"""Mutation test on canonical.py -- are the pre-existing tests stronger?

SAFETY: Uses try/finally to ALWAYS restore the original file, even if the
test script crashes mid-run.  This prevents the script itself from
corrupting the codebase (which is what happened in the first run when
the script left M2_drop_stereo_marker applied permanently).
"""
import subprocess
import re
import sys
from pathlib import Path

# Portable paths -- resolve relative to this file's location
TESTS_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TESTS_DIR.parent
CANON_FILE = SCRIPT_DIR / 'script' / 'canonical.py'
PYTHON = sys.executable

original = CANON_FILE.read_text(encoding='utf-8')

# Mutations to apply to canonical.py
MUTATIONS = [
    ('M1_flip_E_Z_direction',
     "if bond.bond_dir in (1, 3):  # UP from begin->end\n                base = '\\\\' if is_reverse else '/'",
     "if bond.bond_dir in (1, 3):\n                base = '/' if is_reverse else '\\\\'  # MUTATED: flipped",
     'E/Z bond direction inverted (cis becomes trans and vice versa)'),

    ('M2_drop_stereo_marker',
     'chiral_sym = "@R" if bit == 0 else "@S"',
     'chiral_sym = ""  # MUTATED: drop stereo',
     'Stereo marker dropped -- all chiral centers become achiral'),

    ('M3_aromatic_ring_closure_wrong',
     'anubandha = ":" if is_arom else "-"',
     'anubandha = "-" if is_arom else ":"  # MUTATED: swapped',
     'Aromatic and aliphatic ring closures swapped'),

    ('M4_translation_vector_negated',
     'tx, ty, tz = -t[0], -t[1], -t[2]',
     'tx, ty, tz = t[0], t[1], t[2]  # MUTATED: not negated',
     'Periodic translation vector not negated when reversed'),

    ('M5_bond_order_2_to_3',
     "if bt in (BondType.DOUBLE, 2):     valence += 2",
     "if bt in (BondType.DOUBLE, 2):     valence += 3  # MUTATED",
     'Double bond counts as 3 in valence calculation'),
]


def run_tests():
    """Run the canonicalizer + stereo tests."""
    result = subprocess.run(
        [PYTHON, '-m', 'pytest',
         'tests/test_canonicalizer.py', 'tests/test_stereo_correctness.py',
         'tests/test_periodic.py', 'tests/test_v46_gaps.py::TestBridgeBondEdgeCases',
         'tests/test_v46_gaps.py::TestMacrocyclicStereo',
         'tests/test_advanced_features.py',
         '-v', '--tb=line', '-q', '--no-header'],
        cwd=str(SCRIPT_DIR),
        capture_output=True, text=True, timeout=120,
    )
    out = result.stdout + result.stderr
    m = re.search(r'(\d+) passed', out)
    passed = int(m.group(1)) if m else 0
    m = re.search(r'(\d+) failed', out)
    failed = int(m.group(1)) if m else 0
    m = re.search(r'(\d+) error', out)
    errors = int(m.group(1)) if m else 0
    return passed, failed, errors, out


print('=' * 72)
print('  MUTATION TESTING: canonical.py + stereo + periodic tests')
print('=' * 72)

# Baseline
print('\n[baseline] Running tests on unmodified code...')
try:
    p, f, e, out = run_tests()
    print(f'  Result: {p} passed, {f} failed, {e} errors')
    if f > 0 or e > 0:
        print('  WARNING: tests fail on baseline! Aborting.')
        print(out[-1500:])
        sys.exit(1)
except Exception as exn:
    print(f'  baseline crashed: {exn}')
    sys.exit(1)

results = []
try:
    for name, find, replace, desc in MUTATIONS:
        print(f'\n[{name}] {desc}')
        if find not in original:
            print(f'  SKIP: mutation target not found')
            results.append((name, 'SKIP', 0, 0, 0, desc))
            continue
        CANON_FILE.write_text(original.replace(find, replace, 1), encoding='utf-8')
        try:
            p, f, e, out = run_tests()
        except subprocess.TimeoutExpired:
            p, f, e, out = 0, 0, 1, 'TIMEOUT'
        status = 'CAUGHT' if (f + e) > 0 else 'MISSED'
        print(f'  Result: {p} passed, {f} failed, {e} errors -> {status}')
        if f + e == 0:
            print(f'  !!! BUG INTRODUCED BUT TESTS STILL PASS = GAMING !!!')
        else:
            failed_tests = re.findall(r'FAILED (tests/\S+::\S+)', out)
            for ft in failed_tests[:3]:
                print(f'    caught by: {ft}')
        results.append((name, status, p, f, e, desc))
finally:
    # ALWAYS restore original, even on crash
    CANON_FILE.write_text(original, encoding='utf-8')
    print(f'\n[restore] {CANON_FILE.name} restored to original')

# Summary
print('\n' + '=' * 72)
print('  MUTATION TESTING SUMMARY (canonical.py)')
print('=' * 72)
caught = sum(1 for r in results if r[1] == 'CAUGHT')
missed = sum(1 for r in results if r[1] == 'MISSED')
print(f'  Caught: {caught}   Missed: {missed}   Total: {len(results)}')
print()
for name, status, p, f, e, desc in results:
    print(f'  {status:7s} {name:35s} ({f+e} failures)  {desc}')
