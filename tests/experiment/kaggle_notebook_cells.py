# ==============================================================================
# CELL 1: Environment & Dataset Download (Run in Bash)
# ==============================================================================
# Installs RDKit and Lark, then downloads and extracts the latest ChEMBL reps.

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

import os
import sys
import csv
import gzip
import json
import time
import heapq
import hashlib
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

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
# CELL 4: Build Manifest (Data Selection & Stratification - Run in Python)
# ==============================================================================
# Multiprocessed logic for building a stratified 10,000-molecule sample from
# the 2.9 million line ChEMBL representations dataset.

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
            results.append((None, str(e), record_id, smiles))
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

per_stratum = max(1, (target_size + 23) // 24)
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
            heapq.heapreplace(heap, entry)
            
        global_entry = (-score_val, record_id, normalized, bucket)
        if len(global_selected) < target_size:
            heapq.heappush(global_selected, global_entry)
        elif global_entry > global_selected[0]:
            heapq.heapreplace(global_selected, global_entry)

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


# ==============================================================================
# CELL 5: Run Benchmark Roundtrip (Run in Python)
# ==============================================================================
# Run the 10,000 selected molecules from the manifest through the RDKit -> SCRIPT -> RDKit
# roundtrip evaluation in parallel.

def iter_manifest_records(path: Path) -> Iterator[Tuple[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            record_id, smiles = line.split("\t", maxsplit=1)
            yield record_id, smiles

def evaluate_roundtrip(record_id: str, input_smiles: str) -> Dict[str, object]:
    row: Dict[str, object] = {
        "record_id": record_id,
        "input_smiles": input_smiles,
        "encode_ok": False,
        "decode_ok": False,
        "strict_smiles_match": False,
        "inchikey_match": False,
        "canonical_idempotent": False,
        "failure_stage": "",
        "failure_reason": "",
    }
    started = time.perf_counter()
    source = Chem.MolFromSmiles(input_smiles)
    if source is None:
        row.update(failure_stage="input", failure_reason="rdkit_parse_failed")
        return row

    source_smiles = Chem.MolToSmiles(source, canonical=True, isomericSmiles=True)
    source_inchikey = rd_inchi.MolToInchiKey(source)
    row.update(
        source_canonical_smiles=source_smiles,
        source_inchikey=source_inchikey,
        atom_count=source.GetNumAtoms(),
        bond_count=source.GetNumBonds(),
    )

    # Encode SCRIPT
    encode_started = time.perf_counter()
    try:
        script = SCRIPTFromMol(source)
    except Exception as exc:
        row.update(failure_stage="encode", failure_reason=type(exc).__name__)
        return row
    row["encode_ms"] = round((time.perf_counter() - encode_started) * 1000, 3)
    if not script:
        row.update(failure_stage="encode", failure_reason="script_encoding_failed")
        return row
    row.update(encode_ok=True, script=script, script_length=len(script))

    # Decode SCRIPT
    decode_started = time.perf_counter()
    try:
        decoded = MolFromSCRIPT(script)
    except Exception as exc:
        row.update(failure_stage="decode", failure_reason=type(exc).__name__)
        return row
    row["decode_ms"] = round((time.perf_counter() - decode_started) * 1000, 3)
    if decoded is None:
        row.update(failure_stage="decode", failure_reason="script_decoding_failed")
        return row
    row["decode_ok"] = True

    # Fidelity Verification
    decoded_smiles = Chem.MolToSmiles(decoded, canonical=True, isomericSmiles=True)
    decoded_inchikey = rd_inchi.MolToInchiKey(decoded)
    row.update(
        decoded_canonical_smiles=decoded_smiles,
        decoded_inchikey=decoded_inchikey,
        strict_smiles_match=decoded_smiles == source_smiles,
        inchikey_match=decoded_inchikey == source_inchikey,
    )

    try:
        reencoded = SCRIPTFromMol(decoded)
        row["canonical_idempotent"] = reencoded == script
    except Exception:
        row["canonical_idempotent"] = False

    if not row["strict_smiles_match"]:
        row.update(failure_stage="fidelity", failure_reason="isomeric_smiles_mismatch")
    elif not row["inchikey_match"]:
        row.update(failure_stage="fidelity", failure_reason="inchikey_mismatch")
    row["total_ms"] = round((time.perf_counter() - started) * 1000, 3)
    return row

def process_benchmark_batch(batch):
    return [evaluate_roundtrip(record_id, smiles) for record_id, smiles in batch]

# --- Run the Benchmark ---
manifest_input = Path("experiments/data/chembl_pilot_v1.tsv")
results_output = Path("experiments/results/chembl_pilot_v1.jsonl")

# Ensure clean start
if results_output.exists():
    results_output.unlink()
results_output.parent.mkdir(parents=True, exist_ok=True)

print(f"Starting benchmark roundtrip using manifest: {manifest_input}")
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
                print(f"  Benchmarked {processed_bench} records...")

bench_duration = time.perf_counter() - t_bench_start
print(f"Benchmark finished in {bench_duration:.2f} seconds.")


# ==============================================================================
# CELL 6: Build Summary Report & Performance Metrics
# ==============================================================================
# Parse the streaming results and format a final report.

rows = []
with results_output.open("r", encoding="utf-8") as handle:
    for line in handle:
        rows.append(json.loads(line))

counts = Counter(row["failure_stage"] or "pass" for row in rows)
try:
    commit_hash = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
except Exception:
    commit_hash = "unavailable"

summary_report = {
    "manifest": str(manifest_input),
    "manifest_sha256": hashlib.sha256(manifest_input.read_bytes()).hexdigest(),
    "result_file": str(results_output),
    "records": len(rows),
    "rdkit_version": rdBase.rdkitVersion,
    "script_commit": commit_hash,
    "exact_roundtrip": sum(bool(row["strict_smiles_match"]) and bool(row["inchikey_match"]) for row in rows),
    "canonical_idempotent": sum(bool(row["canonical_idempotent"]) for row in rows),
    "failure_stages": dict(sorted(counts.items())),
}

summary_path = results_output.with_suffix(".summary.json")
summary_path.write_text(json.dumps(summary_report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

print(json.dumps(summary_report, indent=2, sort_keys=True))
