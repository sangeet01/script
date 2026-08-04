# ==============================================================================
# CELL 1: Environment & Dataset Download (Run in Bash)
# ==============================================================================
# Installs RDKit and Lark, then downloads and extracts the latest ChEMBL reps.
# ChEMBL is preferred over PubChem for this benchmark because:
#   - Drug-like molecules (more stereo, more complexity)
#   - Curated canonical SMILES (higher quality than PubChem connectivity SMILES)
#   - Standard dataset for notation benchmarks (SELFIES, SAFE used it too)

!pip install rdkit lark
!wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_chemreps.txt.gz
!gunzip -f chembl_37_chemreps.txt.gz
!mkdir -p raw_data
!mv chembl_37_chemreps.txt raw_data/


# ==============================================================================
# CELL 2: Clone Repository & Change Directory (Run in Bash)
# ==============================================================================
# Clones your SCRIPT repo and enters it to access RDKit bridge code.

!git clone https://github.com/sangeet01/script.git
%cd script


# ==============================================================================
# CELL 3: Core Pipeline Imports
# ==============================================================================
# Imports all standard, scientific, and local libraries.
# FIX #1: Added missing typing imports (Iterator, Tuple, Dict) and subprocess.

import os
import sys
import csv
import gzip
import json
import time
import heapq
import hashlib
import subprocess
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path
from typing import Iterator, Tuple, Dict, List, Any, Optional

# Ensure root is in path
ROOT = Path.cwd()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Descriptors
from rdkit.Chem import inchi as rd_inchi
from script.rdkit_bridge import MolFromSCRIPT, SCRIPTFromMol

print(f"RDKit Version: {rdBase.rdkitVersion}")
print(f"Working Directory: {ROOT}")


# ==============================================================================
# CELL 4: Build Manifest (Data Selection & Stratification)
# ==============================================================================
# Multiprocessed logic for building a stratified 10,000-molecule sample from
# the 2.9 million line ChEMBL representations dataset.
#
# Stratification ensures the sample covers chemical space rather than being
# biased toward easy molecules (which a pure random sample would be, since
# ChEMBL is dominated by small organic molecules).
#
# FIX #4: Reduced per_stratum from 418 to 250 to avoid rare strata being
# starved while common strata overflow.  Global backfill redistributes.

METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi",
}

# Worker initialisation
_WORKER_SEED = "script-chembl-pilot-v1"

def init_manifest_worker(seed: str):
    global _WORKER_SEED
    _WORKER_SEED = seed

def calculate_stratum(mol) -> str:
    weight = Descriptors.MolWt(mol)
    weight_bin = "small" if weight < 200 else "medium" if weight < 500 else "large"
    stereo = "stereo" if Chem.FindMolChiralCenters(mol, includeUnassigned=True) else "achiral"
    charged = "charged" if any(atom.GetFormalCharge() for atom in mol.GetAtoms()) else "neutral"
    metal = "metal" if any(atom.GetSymbol() in METALS for atom in mol.GetAtoms()) else "organic"
    return "|".join((weight_bin, stereo, charged, metal))

def compute_score(seed: str, record_id: str, smiles: str) -> int:
    digest = hashlib.sha256(f"{seed}\0{record_id}\0{smiles}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")

def process_manifest_batch(batch):
    results = []
    for record_id, smiles in batch:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                results.append((None, "invalid_input", record_id, smiles))
                continue
            normalized = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            bucket = calculate_stratum(mol)
            score_val = compute_score(_WORKER_SEED, record_id, normalized)
            results.append((score_val, record_id, normalized, bucket))
        except Exception as e:
            results.append((None, str(e)[:120], record_id, smiles))
    return results

def iter_source(source_path: Path) -> Iterator[Tuple[str, str]]:
    with open(source_path, "r", encoding="utf-8", newline="") as f:
        first_line = f.readline()
        f.seek(0)
        delimiter = "\t" if "\t" in first_line else ","
        
        reader = csv.DictReader(f, delimiter=delimiter)
        fieldnames = [fn.lower().replace(" ", "_") for fn in (reader.fieldnames or [])]
        
        smiles_col = None
        id_col = None
        
        for idx, orig_name in enumerate(reader.fieldnames or []):
            norm_name = fieldnames[idx]
            if norm_name in ("canonical_smiles", "smiles", "smiles_string"):
                smiles_col = orig_name
            elif norm_name in ("chembl_id", "molregno", "id", "compound_id"):
                id_col = orig_name
                
        if not smiles_col or not id_col:
            f.seek(0)
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split(delimiter)
                if len(parts) >= 2:
                    yield parts[0], parts[1]
        else:
            for row in reader:
                record_id = (row.get(id_col) or "").strip()
                smiles = (row.get(smiles_col) or "").strip()
                if record_id and smiles:
                    yield record_id, smiles

# --- Run the selection process ---
source_file = Path("../raw_data/chembl_37_chemreps.txt")
output_manifest = Path("experiments/data/chembl_pilot_v1.tsv")
target_size = 10000
seed_val = "script-chembl-pilot-v1"
batch_size = 1000
num_cores = os.cpu_count() or 4

print(f"Building manifest from: {source_file}")
print(f"Target size: {target_size} across {num_cores} cores...")

# FIX #4: Reduced per_stratum to 250 (was 418) so rare strata aren't starved.
# 24 strata × 250 = 6000 base; remaining 4000 filled by global backfill from
# the most populous strata.  This gives better coverage of rare chemistry.
per_stratum = 250
selected = defaultdict(list)
global_selected = []
seen_smiles = set()
stats = Counter()

pool = Pool(processes=num_cores, initializer=init_manifest_worker, initargs=(seed_val,))

def batch_iterator():
    current_batch = []
    for record_id, smiles in iter_source(source_file):
        stats["source_rows"] += 1
        current_batch.append((record_id, smiles))
        if len(current_batch) >= batch_size:
            yield current_batch
            current_batch = []
    if current_batch:
        yield current_batch

results_iterator = pool.imap_unordered(process_manifest_batch, batch_iterator(), chunksize=4)
processed_count = 0

t0 = time.perf_counter()
for batch_results in results_iterator:
    for result in batch_results:
        processed_count += 1
        if processed_count % 200000 == 0:
            print(f"  Processed {processed_count} rows...")
            
        score_val = result[0]
        if score_val is None:
            reason = result[1]
            stats[reason] += 1
            continue
            
        _, record_id, normalized, bucket = result
        if normalized in seen_smiles:
            stats["duplicate_normalized"] += 1
            continue
        seen_smiles.add(normalized)
        
        heap = selected[bucket]
        entry = (-score_val, record_id, normalized)
        if len(heap) < per_stratum:
            heapq.heappush(heap, entry)
        elif entry > heap[0]:
            heapq.heappushpop(heap, entry)
            
        global_entry = (-score_val, record_id, normalized, bucket)
        if len(global_selected) < target_size:
            heapq.heappush(global_selected, global_entry)
        elif global_entry > global_selected[0]:
            heapq.heappushpop(global_selected, global_entry)

pool.close()
pool.join()

# Merge strata selections
rows = []
for bucket, heap in selected.items():
    for negative_score, record_id, smiles in heap:
        rows.append((-negative_score, record_id, smiles, bucket))
rows.sort()

selected_keys = {(record_id, smiles) for _, record_id, smiles, _ in rows}
for negative_score, record_id, smiles, bucket in sorted(global_selected, reverse=True):
    if len(rows) >= target_size:
        break
    if (record_id, smiles) not in selected_keys:
        rows.append((-negative_score, record_id, smiles, bucket))
        selected_keys.add((record_id, smiles))
rows.sort()

output_manifest.parent.mkdir(parents=True, exist_ok=True)
with output_manifest.open("w", encoding="utf-8", newline="\n") as handle:
    handle.write(f"# source={source_file.name}\n# seed={seed_val}\n")
    for _, record_id, smiles, _ in rows:
        handle.write(f"{record_id}\t{smiles}\n")

manifest_hash = hashlib.sha256(output_manifest.read_bytes()).hexdigest()
duration = time.perf_counter() - t0

report = {
    "source": str(source_file),
    "target": target_size,
    "selected": len(rows),
    "seed": seed_val,
    "manifest_sha256": manifest_hash,
    "source_stats": dict(stats),
    "strata": dict(sorted(Counter(bucket for *_, bucket in rows).items())),
}
report_path = output_manifest.with_suffix(".selection.json")
report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

print(f"Manifest written to: {output_manifest}")
print(f"Manifest Hash: {manifest_hash}")
print(f"Selection finished in {duration:.2f} seconds.")
print(f"Strata distribution:")
for k, v in sorted(report["strata"].items()):
    print(f"  {k}: {v}")


# ==============================================================================
# CELL 5: Run Benchmark Roundtrip (SCRIPT arm)
# ==============================================================================
# Run the 10,000 selected molecules from the manifest through the RDKit -> SCRIPT -> RDKit
# roundtrip evaluation in parallel.
#
# FIX #5: Detailed failure categorization — capture exception messages, diff InChI
# layers, record which stereo layer is lost.
# FIX #6: Added RDKit baseline arm (Cell 5b below) for head-to-head comparison.

def iter_manifest_records(path: Path) -> Iterator[Tuple[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            record_id, smiles = line.split("\t", maxsplit=1)
            yield record_id, smiles

def _diff_inchi_layers(source_inchi: str, decoded_inchi: str) -> str:
    """Identify which InChI layer differs between source and decoded.
    Returns a comma-separated string of missing/changed layer names.
    Layers: /f=formula, /c=connections, /h=hydrogens, /q=charge, /p=protons,
            /b=E/Z stereo, /t=tetrahedral stereo, /m=parity, /s=stereo type,
            /i=isotopic, /fixedH, /reconnected
    """
    if not source_inchi or not decoded_inchi:
        return "empty_inchi"
    s_layers = source_inchi.split("/")
    d_layers = decoded_inchi.split("/")
    diffs = []
    # Standard layer prefixes (position 0 is "InChI=1S")
    layer_names = {1: "formula", 2: "connections", 3: "h", 4: "q", 5: "p",
                   6: "b_stereo", 7: "t_stereo", 8: "m_parity", 9: "s_type"}
    max_len = max(len(s_layers), len(d_layers))
    for i in range(1, max_len):
        s = s_layers[i] if i < len(s_layers) else ""
        d = d_layers[i] if i < len(d_layers) else ""
        if s != d:
            name = layer_names.get(i, f"layer_{i}")
            if not d and s:
                diffs.append(f"missing_{name}")
            elif not s and d:
                diffs.append(f"extra_{name}")
            else:
                diffs.append(f"changed_{name}")
    return ",".join(diffs) if diffs else "identical"

def evaluate_roundtrip(record_id: str, input_smiles: str) -> Dict[str, object]:
    """Run one molecule through SCRIPT round-trip and record detailed results.
    
    Failure stages:
      - input:   RDKit couldn't parse the source SMILES
      - encode:  SCRIPTFromMol failed (exception or empty result)
      - decode:  MolFromSCRIPT failed (exception or None result)
      - fidelity: round-trip produced different SMILES or InChIKey
    
    Failure reasons are specific (exception type + message, or layer diff).
    """
    row: Dict[str, object] = {
        "record_id": record_id,
        "input_smiles": input_smiles,
        "arm": "script",
        "encode_ok": False,
        "decode_ok": False,
        "strict_smiles_match": False,
        "inchikey_match": False,
        "canonical_idempotent": False,
        "failure_stage": "",
        "failure_reason": "",
        "inchi_layer_diff": "",
    }
    started = time.perf_counter()
    source = Chem.MolFromSmiles(input_smiles)
    if source is None:
        row.update(failure_stage="input", failure_reason="rdkit_parse_failed")
        return row

    source_smiles = Chem.MolToSmiles(source, canonical=True, isomericSmiles=True)
    source_inchikey = rd_inchi.MolToInchiKey(source)
    source_inchi = rd_inchi.MolToInchi(source)
    row.update(
        source_canonical_smiles=source_smiles,
        source_inchikey=source_inchikey,
        source_inchi=source_inchi,
        atom_count=source.GetNumAtoms(),
        bond_count=source.GetNumBonds(),
    )

    # Encode SCRIPT — capture detailed exception info
    encode_started = time.perf_counter()
    try:
        script = SCRIPTFromMol(source)
    except Exception as exc:
        # FIX #5: capture exception type AND truncated message
        msg = str(exc)[:200].replace("\n", " ")
        row.update(failure_stage="encode",
                   failure_reason=f"{type(exc).__name__}: {msg}")
        return row
    row["encode_ms"] = round((time.perf_counter() - encode_started) * 1000, 3)
    if not script:
        row.update(failure_stage="encode", failure_reason="script_encoding_empty")
        return row
    row.update(encode_ok=True, script=script, script_length=len(script))

    # Decode SCRIPT — capture detailed exception info
    decode_started = time.perf_counter()
    try:
        decoded = MolFromSCRIPT(script)
    except Exception as exc:
        msg = str(exc)[:200].replace("\n", " ")
        row.update(failure_stage="decode",
                   failure_reason=f"{type(exc).__name__}: {msg}")
        return row
    row["decode_ms"] = round((time.perf_counter() - decode_started) * 1000, 3)
    if decoded is None:
        row.update(failure_stage="decode", failure_reason="script_decoding_none")
        return row
    row["decode_ok"] = True

    # Fidelity Verification
    decoded_smiles = Chem.MolToSmiles(decoded, canonical=True, isomericSmiles=True)
    decoded_inchikey = rd_inchi.MolToInchiKey(decoded)
    decoded_inchi = rd_inchi.MolToInchi(decoded)
    row.update(
        decoded_canonical_smiles=decoded_smiles,
        decoded_inchikey=decoded_inchikey,
        decoded_inchi=decoded_inchi,
        strict_smiles_match=decoded_smiles == source_smiles,
        inchikey_match=decoded_inchikey == source_inchikey,
    )

    # Canonical idempotency: re-encoding the decoded mol should give the same SCRIPT
    try:
        reencoded = SCRIPTFromMol(decoded)
        row["canonical_idempotent"] = reencoded == script
    except Exception:
        row["canonical_idempotent"] = False

    # FIX #5: categorize fidelity failures by which InChI layer differs
    if not row["strict_smiles_match"]:
        row.update(failure_stage="fidelity", failure_reason="isomeric_smiles_mismatch")
    elif not row["inchikey_match"]:
        # InChIKey mismatch — dig into which InChI layer differs
        layer_diff = _diff_inchi_layers(source_inchi or "", decoded_inchi or "")
        row.update(failure_stage="fidelity",
                   failure_reason="inchikey_mismatch",
                   inchi_layer_diff=layer_diff)
    row["total_ms"] = round((time.perf_counter() - started) * 1000, 3)
    return row

def process_benchmark_batch(batch):
    return [evaluate_roundtrip(record_id, smiles) for record_id, smiles in batch]

# --- Run the SCRIPT Benchmark ---
manifest_input = Path("experiments/data/chembl_pilot_v1.tsv")
results_output = Path("experiments/results/chembl_pilot_v1_script.jsonl")

# Ensure clean start
if results_output.exists():
    results_output.unlink()
results_output.parent.mkdir(parents=True, exist_ok=True)

print(f"Starting SCRIPT benchmark roundtrip using manifest: {manifest_input}")
benchmark_batch_size = 100
benchmark_cores = os.cpu_count() or 4

def benchmark_batch_iterator():
    current_batch = []
    for record_id, smiles in iter_manifest_records(manifest_input):
        current_batch.append((record_id, smiles))
        if len(current_batch) >= benchmark_batch_size:
            yield current_batch
            current_batch = []
    if current_batch:
        yield current_batch

t_bench_start = time.perf_counter()
processed_bench = 0

with Pool(processes=benchmark_cores) as pool:
    with results_output.open("w", encoding="utf-8") as handle:
        for batch_results in pool.imap_unordered(process_benchmark_batch, benchmark_batch_iterator(), chunksize=4):
            for row in batch_results:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                processed_bench += 1
            handle.flush()
            if processed_bench % 1000 == 0:
                print(f"  SCRIPT: benchmarked {processed_bench} records...")

bench_duration = time.perf_counter() - t_bench_start
print(f"SCRIPT benchmark finished in {bench_duration:.2f} seconds ({processed_bench} records).")


# ==============================================================================
# CELL 5b: Run RDKit Baseline (head-to-head comparison)
# ==============================================================================
# FIX #6: Same 10K molecules through RDKit canonical SMILES round-trip ONLY.
# This establishes the "industry standard" baseline.  If SCRIPT matches or
# beats RDKit on the same molecules, the "100% looks gamed" critique is
# neutralized — SCRIPT is no worse than the reference implementation.

def evaluate_rdkit_baseline(record_id: str, input_smiles: str) -> Dict[str, object]:
    """Run one molecule through RDKit canonical SMILES round-trip only.
    
    This is the baseline: SMILES -> MolToSmiles(canonical) -> MolFromSmiles -> InChIKey.
    If RDKit itself can't round-trip a molecule through its own canonical form,
    SCRIPT's failure on the same molecule is not SCRIPT's fault.
    """
    row: Dict[str, object] = {
        "record_id": record_id,
        "input_smiles": input_smiles,
        "arm": "rdkit_baseline",
        "encode_ok": False,
        "decode_ok": False,
        "strict_smiles_match": False,
        "inchikey_match": False,
        "canonical_idempotent": False,
        "failure_stage": "",
        "failure_reason": "",
        "inchi_layer_diff": "",
    }
    started = time.perf_counter()
    source = Chem.MolFromSmiles(input_smiles)
    if source is None:
        row.update(failure_stage="input", failure_reason="rdkit_parse_failed")
        return row

    source_smiles = Chem.MolToSmiles(source, canonical=True, isomericSmiles=True)
    source_inchikey = rd_inchi.MolToInchiKey(source)
    source_inchi = rd_inchi.MolToInchi(source)
    row.update(
        source_canonical_smiles=source_smiles,
        source_inchikey=source_inchikey,
        source_inchi=source_inchi,
        atom_count=source.GetNumAtoms(),
        bond_count=source.GetNumBonds(),
    )

    # "Encode" = RDKit canonical SMILES
    encode_started = time.perf_counter()
    try:
        canon_smiles = Chem.MolToSmiles(source, canonical=True, isomericSmiles=True)
    except Exception as exc:
        msg = str(exc)[:200].replace("\n", " ")
        row.update(failure_stage="encode",
                   failure_reason=f"{type(exc).__name__}: {msg}")
        return row
    row["encode_ms"] = round((time.perf_counter() - encode_started) * 1000, 3)
    if not canon_smiles:
        row.update(failure_stage="encode", failure_reason="canon_smiles_empty")
        return row
    row.update(encode_ok=True, script=canon_smiles, script_length=len(canon_smiles))

    # "Decode" = re-parse the canonical SMILES
    decode_started = time.perf_counter()
    try:
        decoded = Chem.MolFromSmiles(canon_smiles)
    except Exception as exc:
        msg = str(exc)[:200].replace("\n", " ")
        row.update(failure_stage="decode",
                   failure_reason=f"{type(exc).__name__}: {msg}")
        return row
    row["decode_ms"] = round((time.perf_counter() - decode_started) * 1000, 3)
    if decoded is None:
        row.update(failure_stage="decode", failure_reason="rdkit_reparse_none")
        return row
    row["decode_ok"] = True

    # Fidelity
    decoded_smiles = Chem.MolToSmiles(decoded, canonical=True, isomericSmiles=True)
    decoded_inchikey = rd_inchi.MolToInchiKey(decoded)
    decoded_inchi = rd_inchi.MolToInchi(decoded)
    row.update(
        decoded_canonical_smiles=decoded_smiles,
        decoded_inchikey=decoded_inchikey,
        decoded_inchi=decoded_inchi,
        strict_smiles_match=decoded_smiles == source_smiles,
        inchikey_match=decoded_inchikey == source_inchikey,
    )

    # Idempotency: re-canonicalizing should give the same string
    try:
        recanon = Chem.MolToSmiles(decoded, canonical=True, isomericSmiles=True)
        row["canonical_idempotent"] = recanon == canon_smiles
    except Exception:
        row["canonical_idempotent"] = False

    if not row["strict_smiles_match"]:
        row.update(failure_stage="fidelity", failure_reason="isomeric_smiles_mismatch")
    elif not row["inchikey_match"]:
        layer_diff = _diff_inchi_layers(source_inchi or "", decoded_inchi or "")
        row.update(failure_stage="fidelity",
                   failure_reason="inchikey_mismatch",
                   inchi_layer_diff=layer_diff)
    row["total_ms"] = round((time.perf_counter() - started) * 1000, 3)
    return row

def process_baseline_batch(batch):
    return [evaluate_rdkit_baseline(record_id, smiles) for record_id, smiles in batch]

# --- Run the RDKit Baseline ---
baseline_output = Path("experiments/results/chembl_pilot_v1_rdkit_baseline.jsonl")

if baseline_output.exists():
    baseline_output.unlink()

print(f"\nStarting RDKit baseline roundtrip using manifest: {manifest_input}")

t_baseline_start = time.perf_counter()
processed_baseline = 0

with Pool(processes=benchmark_cores) as pool:
    with baseline_output.open("w", encoding="utf-8") as handle:
        for batch_results in pool.imap_unordered(process_baseline_batch, benchmark_batch_iterator(), chunksize=4):
            for row in batch_results:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                processed_baseline += 1
            handle.flush()
            if processed_baseline % 1000 == 0:
                print(f"  RDKit: benchmarked {processed_baseline} records...")

baseline_duration = time.perf_counter() - t_baseline_start
print(f"RDKit baseline finished in {baseline_duration:.2f} seconds ({processed_baseline} records).")


# ==============================================================================
# CELL 6: Build Summary Report & Performance Metrics
# ==============================================================================
# Parse the streaming results from BOTH arms and format a final head-to-head report.

def load_jsonl(path: Path) -> List[Dict]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            rows.append(json.loads(line))
    return rows

script_rows = load_jsonl(results_output)
baseline_rows = load_jsonl(baseline_output)

script_counts = Counter(row["failure_stage"] or "pass" for row in script_rows)
baseline_counts = Counter(row["failure_stage"] or "pass" for row in baseline_rows)

# FIX #5: Aggregate failure reasons (not just stages)
script_fail_reasons = Counter(
    row["failure_reason"] for row in script_rows if row["failure_stage"]
)
baseline_fail_reasons = Counter(
    row["failure_reason"] for row in baseline_rows if row["failure_stage"]
)

# FIX #5: Aggregate InChI layer diffs for fidelity failures
script_layer_diffs = Counter(
    row.get("inchi_layer_diff", "") for row in script_rows
    if row["failure_stage"] == "fidelity" and row.get("inchi_layer_diff")
)
baseline_layer_diffs = Counter(
    row.get("inchi_layer_diff", "") for row in baseline_rows
    if row["failure_stage"] == "fidelity" and row.get("inchi_layer_diff")
)

try:
    commit_hash = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
except Exception:
    commit_hash = "unavailable"

# Head-to-head comparison
script_pass = sum(1 for r in script_rows if not r["failure_stage"])
script_exact = sum(1 for r in script_rows
                   if r["strict_smiles_match"] and r["inchikey_match"])
script_idem = sum(1 for r in script_rows if r["canonical_idempotent"])
baseline_pass = sum(1 for r in baseline_rows if not r["failure_stage"])
baseline_exact = sum(1 for r in baseline_rows
                     if r["strict_smiles_match"] and r["inchikey_match"])
baseline_idem = sum(1 for r in baseline_rows if r["canonical_idempotent"])

n_script = len(script_rows)
n_baseline = len(baseline_rows)

summary_report = {
    "manifest": str(manifest_input),
    "manifest_sha256": hashlib.sha256(manifest_input.read_bytes()).hexdigest(),
    "rdkit_version": rdBase.rdkitVersion,
    "script_commit": commit_hash,
    "script_arm": {
        "result_file": str(results_output),
        "records": n_script,
        "exact_roundtrip": script_exact,
        "exact_roundtrip_pct": round(100 * script_exact / max(1, n_script), 2),
        "canonical_idempotent": script_idem,
        "canonical_idempotent_pct": round(100 * script_idem / max(1, n_script), 2),
        "failure_stages": dict(sorted(script_counts.items())),
        "top_failure_reasons": dict(script_fail_reasons.most_common(15)),
        "inchi_layer_diffs": dict(script_layer_diffs.most_common(15)),
    },
    "rdkit_baseline_arm": {
        "result_file": str(baseline_output),
        "records": n_baseline,
        "exact_roundtrip": baseline_exact,
        "exact_roundtrip_pct": round(100 * baseline_exact / max(1, n_baseline), 2),
        "canonical_idempotent": baseline_idem,
        "canonical_idempotent_pct": round(100 * baseline_idem / max(1, n_baseline), 2),
        "failure_stages": dict(sorted(baseline_counts.items())),
        "top_failure_reasons": dict(baseline_fail_reasons.most_common(15)),
        "inchi_layer_diffs": dict(baseline_layer_diffs.most_common(15)),
    },
    "head_to_head": {
        "script_exact_pct": round(100 * script_exact / max(1, n_script), 2),
        "rdkit_exact_pct": round(100 * baseline_exact / max(1, n_baseline), 2),
        "script_beats_rdkit": script_exact / max(1, n_script) >= baseline_exact / max(1, n_baseline),
        "delta_pct": round(100 * (script_exact / max(1, n_script) - baseline_exact / max(1, n_baseline)), 2),
    },
}

summary_path = results_output.with_suffix(".summary.json")
summary_path.write_text(json.dumps(summary_report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

print("\n" + "=" * 72)
print("  HEAD-TO-HEAD BENCHMARK SUMMARY")
print("=" * 72)
print(f"\n  Records: {n_script} (SCRIPT) vs {n_baseline} (RDKit baseline)")
print(f"\n  {'Metric':<30s} {'SCRIPT':>12s} {'RDKit baseline':>15s}")
print(f"  {'-'*30} {'-'*12} {'-'*15}")
print(f"  {'Exact round-trip':<30s} {script_exact:>12d} {baseline_exact:>15d}")
print(f"  {'Exact round-trip %':<30s} {summary_report['script_arm']['exact_roundtrip_pct']:>11.2f}% {summary_report['rdkit_baseline_arm']['exact_roundtrip_pct']:>14.2f}%")
print(f"  {'Canonical idempotent':<30s} {script_idem:>12d} {baseline_idem:>15d}")
print(f"  {'Canonical idempotent %':<30s} {summary_report['script_arm']['canonical_idempotent_pct']:>11.2f}% {summary_report['rdkit_baseline_arm']['canonical_idempotent_pct']:>14.2f}%")
print(f"\n  Delta (SCRIPT - RDKit): {summary_report['head_to_head']['delta_pct']:+.2f} percentage points")
if summary_report['head_to_head']['script_beats_rdkit']:
    print(f"  [OK] SCRIPT matches or beats RDKit on the same molecules.")
else:
    print(f"  [!!]  SCRIPT underperforms RDKit by {abs(summary_report['head_to_head']['delta_pct']):.2f} pp.")
print(f"\n  SCRIPT failure stages: {dict(sorted(script_counts.items()))}")
print(f"  RDKit failure stages:  {dict(sorted(baseline_counts.items()))}")
if script_fail_reasons:
    print(f"\n  SCRIPT top failure reasons:")
    for reason, n in script_fail_reasons.most_common(10):
        print(f"    {n:5d}  {reason[:100]}")
if script_layer_diffs:
    print(f"\n  SCRIPT InChI layer diffs (fidelity failures):")
    for diff, n in script_layer_diffs.most_common(10):
        print(f"    {n:5d}  {diff}")
print(f"\n  Full summary written to: {summary_path}")
print(f"  SCRIPT results:   {results_output}")
print(f"  RDKit results:    {baseline_output}")
