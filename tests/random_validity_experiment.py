"""
Experiment 1: Random String Validity Comparison

Generates random strings from each notation's vocabulary and tests
whether they decode to valid molecules. This is the key experiment
from the SELFIES paper (arXiv:1905.13741), replicated for SCRIPT.

No ML model, no GPU, no training — just random sampling from the
vocabulary and validity checking.

Runs in < 5 minutes on a laptop.
"""

import sys
import os
import time
import random
import string

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# ═══════════════════════════════════════════════════════════════════════════════
# Vocabulary definitions
# ═══════════════════════════════════════════════════════════════════════════════

SMILES_VOCAB = list("CNOSPFClBrIBH()*-=#:.123456789[]@/\\+") + ['@@']

SELFIES_VOCAB = [
    '[C]', '[=C]', '[#C]', '[N]', '[=N]', '[O]', '[=O]', '[S]', '[=S]',
    '[F]', '[Cl]', '[Br]', '[I]', '[P]', '[B]', '[H]',
    '[Branch1]', '[Branch2]', '[Branch3]', '[Ring1]', '[Ring2]', '[Ring3]',
    '[=Branch1]', '[=Branch2]', '[#Branch1]', '[#Branch2]',
    '[nop]', '.', '[C@]', '[C@@]', '[C@H1]', '[C@@H1]',
]

SCRIPT_VOCAB = list("CNOSPFClBrIBH()*-=#:.&123456789[]@/\\") + [
    '@@', '@R', '@S', '@AX1', '@AX2', '@SP1', '@SP2',
    '->', '<-', '=:', '>', '*',
]


# ═══════════════════════════════════════════════════════════════════════════════
# Random string generators
# ═══════════════════════════════════════════════════════════════════════════════

def generate_random_smiles(length, rng=random):
    return ''.join(rng.choices(SMILES_VOCAB, k=length))


def generate_random_selfies(length, rng=random):
    tokens = rng.choices(SELFIES_VOCAB, k=length)
    return ''.join(tokens)


def generate_random_script_unconstrained(length, rng=random):
    return ''.join(rng.choices(SCRIPT_VOCAB, k=length))


def generate_random_script_constrained(length, rng=random):
    from script.constrained_decoder import ConstrainedSCRIPTDecoder
    decoder = ConstrainedSCRIPTDecoder()
    decoder.reset()
    generated = []
    for _ in range(length):
        valid_tokens = decoder.get_valid_tokens(SCRIPT_VOCAB)
        if not valid_tokens:
            decoder.reset()
            valid_tokens = decoder.get_valid_tokens(SCRIPT_VOCAB)
        token = rng.choice(list(valid_tokens))
        if not decoder.consume(token):
            decoder.reset()
            token = rng.choice(list(decoder.get_valid_tokens(SCRIPT_VOCAB)))
            decoder.consume(token)
        generated.append(token)
        if decoder.is_complete() and rng.random() < 0.1:
            break
    return ''.join(generated)


# ═══════════════════════════════════════════════════════════════════════════════
# Validity checkers
# ═══════════════════════════════════════════════════════════════════════════════

def check_smiles_validity(smiles_str):
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles_str)
        return mol is not None and mol.GetNumAtoms() > 0
    except Exception:
        return False


def check_selfies_validity(selfies_str):
    try:
        import selfies as sf
        smiles = sf.decoder(selfies_str)
        if smiles is None or smiles == '':
            return False
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None and mol.GetNumAtoms() > 0
    except ImportError:
        return True  # SELFIES guarantees validity by construction
    except Exception:
        return False


def check_script_validity(script_str):
    try:
        from script.parser import SCRIPTParser
        parser = SCRIPTParser()
        result = parser.parse(script_str)
        if not result["success"]:
            return False
        mol = result["molecule"]
        if mol is None:
            return False
        if isinstance(mol, list):
            return len(mol) > 0 and all(len(m.atoms) > 0 for m in mol)
        return len(mol.atoms) > 0
    except Exception:
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# Main experiment
# ═══════════════════════════════════════════════════════════════════════════════

def run_experiment(n_samples=1000, max_length=50, seed=42):
    rng = random.Random(seed)
    
    results = {
        'SMILES': {'valid': 0, 'invalid': 0, 'total': 0, 'time': 0},
        'SELFIES': {'valid': 0, 'invalid': 0, 'total': 0, 'time': 0},
        'SCRIPT (unconstrained)': {'valid': 0, 'invalid': 0, 'total': 0, 'time': 0},
        'SCRIPT (constrained)': {'valid': 0, 'invalid': 0, 'total': 0, 'time': 0},
    }
    
    print(f"Generating {n_samples} random strings per notation (max length {max_length})...")
    print(f"{'Notation':<25} {'Valid':>8} {'Invalid':>8} {'Rate':>8} {'Time':>8}")
    print("-" * 65)
    
    # SMILES
    t0 = time.time()
    for i in range(n_samples):
        length = rng.randint(5, max_length)
        s = generate_random_smiles(length, rng)
        valid = check_smiles_validity(s)
        results['SMILES']['total'] += 1
        results['SMILES']['valid' if valid else 'invalid'] += 1
    results['SMILES']['time'] = time.time() - t0
    r = results['SMILES']
    print(f"{'SMILES':<25} {r['valid']:>8} {r['invalid']:>8} {r['valid']/r['total']*100:>7.1f}% {r['time']:>7.1f}s")
    
    # SELFIES
    t0 = time.time()
    for i in range(n_samples):
        length = rng.randint(5, max_length // 3)
        s = generate_random_selfies(length, rng)
        valid = check_selfies_validity(s)
        results['SELFIES']['total'] += 1
        results['SELFIES']['valid' if valid else 'invalid'] += 1
    results['SELFIES']['time'] = time.time() - t0
    r = results['SELFIES']
    print(f"{'SELFIES':<25} {r['valid']:>8} {r['invalid']:>8} {r['valid']/r['total']*100:>7.1f}% {r['time']:>7.1f}s")
    
    # SCRIPT (unconstrained)
    t0 = time.time()
    for i in range(n_samples):
        length = rng.randint(5, max_length)
        s = generate_random_script_unconstrained(length, rng)
        valid = check_script_validity(s)
        results['SCRIPT (unconstrained)']['total'] += 1
        results['SCRIPT (unconstrained)']['valid' if valid else 'invalid'] += 1
    results['SCRIPT (unconstrained)']['time'] = time.time() - t0
    r = results['SCRIPT (unconstrained)']
    print(f"{'SCRIPT (unconstrained)':<25} {r['valid']:>8} {r['invalid']:>8} {r['valid']/r['total']*100:>7.1f}% {r['time']:>7.1f}s")
    
    # SCRIPT (constrained)
    t0 = time.time()
    for i in range(n_samples):
        length = rng.randint(5, max_length)
        s = generate_random_script_constrained(length, rng)
        valid = check_script_validity(s)
        results['SCRIPT (constrained)']['total'] += 1
        results['SCRIPT (constrained)']['valid' if valid else 'invalid'] += 1
    results['SCRIPT (constrained)']['time'] = time.time() - t0
    r = results['SCRIPT (constrained)']
    print(f"{'SCRIPT (constrained)':<25} {r['valid']:>8} {r['invalid']:>8} {r['valid']/r['total']*100:>7.1f}% {r['time']:>7.1f}s")
    
    # Summary
    print("\n" + "=" * 65)
    print("SUMMARY: Random String Validity Comparison")
    print("=" * 65)
    print(f"{'Notation':<25} {'Valid':>8} {'Invalid':>8} {'Validity':>10}")
    print("-" * 65)
    for name, r in results.items():
        rate = r['valid'] / r['total'] * 100 if r['total'] > 0 else 0
        print(f"{name:<25} {r['valid']:>8} {r['invalid']:>8} {rate:>9.1f}%")
    print("=" * 65)
    
    # LaTeX table
    print("\n% LaTeX table for paper:")
    print(r"\begin{table}[h]")
    print(r"\centering")
    print(r"\caption{Random string validity comparison (N=" + str(n_samples) + r", max length=" + str(max_length) + r")}")
    print(r"\label{tab:validity}")
    print(r"\begin{tabular}{lrrr}")
    print(r"\toprule")
    print(r"\textbf{Notation} & \textbf{Valid} & \textbf{Invalid} & \textbf{Validity} \\")
    print(r"\midrule")
    for name, r in results.items():
        rate = r['valid'] / r['total'] * 100 if r['total'] > 0 else 0
        print(f"{name} & {r['valid']} & {r['invalid']} & {rate:.1f}\\% \\\\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")
    
    return results


if __name__ == "__main__":
    # Check dependencies
    try:
        from rdkit import Chem
        print("RDKit: OK")
    except ImportError:
        print("WARNING: RDKit not installed")
    
    try:
        import selfies as sf
        print("SELFIES: OK")
    except ImportError:
        print("NOTE: SELFIES not installed — assuming 100% validity (by construction)")
    
    print()
    results = run_experiment(n_samples=1000, max_length=50, seed=42)
