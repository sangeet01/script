"""
Mutation Testing: Are the V4.6 graft copolymer tests real or gaming?

We deliberately introduce 5 bugs into graft_expander.py (one at a time),
re-run the test suite, and check whether tests catch each bug.  If a bug
goes undetected, the tests are gaming.  If all bugs are caught, the tests
are real.

Mutations applied:
  M1: Off-by-one in backbone branch point assignment (use i+1 instead of i)
  M2: Wrong bond order for graft attachment (DOUBLE instead of SINGLE)
  M3: Skip the last graft copy (range(1, M) -> range(1, M-1))
  M4: Don't remove the seam bond (skip _remove_bond_between call)
  M5: Attach surplus grafts to FIRST backbone atom instead of LAST

For each mutation:
  - Apply the change to a copy of graft_expander.py
  - Run pytest tests/test_v46_gaps.py::TestGraftCopolymer
  - Record: # tests passed, # tests failed, which tests caught it
  - Revert the change

A mutation is "caught" if at least one test fails.
A mutation is "missed" if all tests still pass -- that is gaming.

SAFETY: Uses try/finally to ALWAYS restore the original file, even if the
script crashes mid-run.
"""
import subprocess
import re
import sys
from pathlib import Path

# Portable paths -- resolve relative to this file's location
TESTS_DIR = Path(__file__).resolve().parent
SCRIPT_DIR = TESTS_DIR.parent
GRAFT_FILE = SCRIPT_DIR / 'script' / 'graft_expander.py'
PYTHON = sys.executable

# Save original
original = GRAFT_FILE.read_text(encoding='utf-8')

# Define mutations: (name, find, replace, description)
MUTATIONS = [
    ('M1_off_by_one_branch_point',
     'bb_branch_points.append(new_head)',
     'bb_branch_points.append(new_head + 0)  # MUTATED: skip first branch point\n        if len(bb_branch_points) > 1: bb_branch_points[-1] = bb_branch_points[0]',
     'Off-by-one: surplus grafts attach to first backbone atom instead of last'),

    ('M2_wrong_attachment_bond_order',
     'mol.add_bond(target_bp, gr_head, BondType.SINGLE, bond_class="")',
     'mol.add_bond(target_bp, gr_head, BondType.DOUBLE, bond_class="")  # MUTATED',
     'Graft attachment uses DOUBLE bond instead of SINGLE'),

    ('M3_skip_last_graft_copy',
     'for i in range(1, M):',
     'for i in range(1, max(1, M-1)):  # MUTATED: skip last graft',
     'Skip the last graft copy (M-1 grafts instead of M)'),

    ('M4_skip_seam_removal',
     '_remove_bond_between(mol, bb_end, gr_start)',
     'pass  # MUTATED: skip seam bond removal',
     "Don't remove the seam bond between backbone tail and graft head"),

    ('M5_wrong_surplus_attachment',
     'target_bp = bb_branch_points[-1] if bb_branch_points else bb_start',
     'target_bp = bb_branch_points[0] if bb_branch_points else bb_start  # MUTATED',
     'Surplus grafts attach to FIRST backbone atom instead of LAST'),
]


def run_tests():
    """Run the graft copolymer tests, return (passed, failed, output)."""
    result = subprocess.run(
        [PYTHON, '-m', 'pytest',
         'tests/test_v46_gaps.py::TestGraftCopolymer',
         '-v', '--tb=line', '-q'],
        cwd=str(SCRIPT_DIR),
        capture_output=True, text=True, timeout=60,
    )
    out = result.stdout + result.stderr
    m = re.search(r'(\d+) passed', out)
    passed = int(m.group(1)) if m else 0
    m = re.search(r'(\d+) failed', out)
    failed = int(m.group(1)) if m else 0
    return passed, failed, out


def apply_mutation(find, replace):
    """Apply a mutation if the find string exists. Return True if applied."""
    if find not in original:
        return False
    mutated = original.replace(find, replace, 1)
    GRAFT_FILE.write_text(mutated, encoding='utf-8')
    return True


print('=' * 72)
print('  MUTATION TESTING: Are the V4.6 graft tests real or gaming?')
print('=' * 72)

# Baseline: run tests on unmodified code
print('\n[baseline] Running tests on unmodified code...')
passed, failed, out = run_tests()
print(f'  Result: {passed} passed, {failed} failed')
if failed > 0:
    print('  WARNING: tests fail on baseline! Aborting.')
    sys.exit(1)

# Apply each mutation
results = []
try:
    for name, find, replace, desc in MUTATIONS:
        print(f'\n[{name}] {desc}')
        if not apply_mutation(find, replace):
            print(f'  SKIP: mutation target not found in source')
            results.append((name, 'SKIP', 0, 0, 'mutation target not found'))
            continue
        passed, failed, out = run_tests()
        status = 'CAUGHT' if failed > 0 else 'MISSED'
        print(f'  Result: {passed} passed, {failed} failed -> {status}')
        if failed == 0:
            print(f'  !!! BUG INTRODUCED BUT TESTS STILL PASS = GAMING !!!')
        else:
            failed_tests = re.findall(r'FAILED (tests/test_v46_gaps.py::\S+)', out)
            for ft in failed_tests[:3]:
                print(f'    caught by: {ft}')
        results.append((name, status, passed, failed, desc))
finally:
    # ALWAYS restore, even on crash
    GRAFT_FILE.write_text(original, encoding='utf-8')
    print(f'\n[restore] {GRAFT_FILE.name} restored to original')

# Summary
print('\n' + '=' * 72)
print('  MUTATION TESTING SUMMARY')
print('=' * 72)
caught = sum(1 for r in results if r[1] == 'CAUGHT')
missed = sum(1 for r in results if r[1] == 'MISSED')
skipped = sum(1 for r in results if r[1] == 'SKIP')
print(f'  Total mutations:    {len(results)}')
print(f'  Caught (real test): {caught}')
print(f'  Missed (gaming):    {missed}')
print(f'  Skipped:            {skipped}')
print()
if missed == 0 and caught > 0:
    print(f'  VERDICT: Tests are REAL -- every injected bug was caught.')
elif missed > 0:
    print(f'  VERDICT: Tests have GAPS -- {missed} bug(s) went undetected.')
    print(f'  The tests are partially gaming: they confirm the happy path')
    print(f'  but miss real failure modes.')
else:
    print(f'  VERDICT: Inconclusive (no mutations applied).')
print()
for name, status, p, f, desc in results:
    print(f'  {status:7s} {name:35s} ({f} test failures)  {desc}')
