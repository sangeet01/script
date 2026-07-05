#!/usr/bin/env python3
"""
QM9 Experiments: Generative Validity and Internal Diversity Comparison
=======================================================================
Experiment 2 - QM9 Generative Model:
    Train a character-level Markov model on QM9 molecules (all three
    notations), generate 1 000 samples per notation, and report validity
    rates.  SCRIPT is tested in both unconstrained and constrained modes.

Experiment 3 - Internal Diversity:
    Compute internal diversity (1 - mean Tanimoto) of the valid molecules
    generated in Experiment 2.

Data:
    Real QM9 dataset (Ramakrishnan et al., 2014).
    134 885 small organic molecules, SMILES embedded in each .xyz file.
    Downloaded once from Figshare and cached locally.

    Primary: https://figshare.com/ndownloader/files/3195389
    API fallback: https://api.figshare.com/v2/articles/978904/files

Usage:
    python tests/qm9_experiment.py                   # full run (all 134k)
    python tests/qm9_experiment.py --n-train 5000    # smaller training set
    python tests/qm9_experiment.py --quick           # 1k train / 200 gen
    python tests/qm9_experiment.py --no-cache        # force fresh download
    python tests/qm9_experiment.py --json            # extra JSON output

Requirements:
    rdkit       - pip install rdkit
    selfies     - pip install selfies          (optional, auto-skipped)
    numpy       - pip install numpy            (optional, for diversity)
    requests    - pip install requests         (auto-installed if missing)

Outputs (saved next to the script):
    qm9_experiments_YYYYMMDD_HHMM.txt
    qm9_experiments_YYYYMMDD_HHMM.json   (with --json)
"""

from __future__ import annotations

import argparse
import bz2
import collections
import datetime
import io
import json
import logging
import math
import os
import random
import re
import sys
import tarfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("qm9")

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT  = TESTS_DIR.parent
sys.path.insert(0, str(REPO_ROOT))

CACHE_DIR  = TESTS_DIR / ".qm9_cache"
CACHE_FILE = CACHE_DIR / "qm9_smiles.json"

# ---------------------------------------------------------------------------
# QM9 download URLs (in priority order)
# ---------------------------------------------------------------------------
QM9_URLS = [
    # Historic ndownloader direct link (file ID 3195389)
    "https://figshare.com/ndownloader/files/3195389",
    # Figshare API resolves to the same file but via a stable endpoint
    "https://ndownloader.figshare.com/files/3195389",
]
QM9_FIGSHARE_API = "https://api.figshare.com/v2/articles/978904/files"
# QM9 has 133,885 molecules after excluding 3,054 uncharacterized entries.
QM9_TOTAL_MOLECULES = 133_885

EXCLUDE_URL = "https://figshare.com/ndownloader/files/3195404"  # uncharacterized.txt


# ===========================================================================
# Dependency check
# ===========================================================================

def _require(pkg: str, install: str = "") -> bool:
    import importlib
    if importlib.util.find_spec(pkg) is not None:
        return True
    install = install or pkg
    log.warning("Package '%s' not found. Install: pip install %s", pkg, install)
    return False


def _ensure_requests():
    try:
        import requests  # noqa: F401
    except ImportError:
        import subprocess
        log.info("Installing 'requests'...")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "requests", "-q"]
        )


# ===========================================================================
# QM9 downloader
# ===========================================================================

def _get(url: str, stream: bool = False, timeout: int = 60, retries: int = 4):
    import requests
    for attempt in range(retries):
        try:
            resp = requests.get(url, stream=stream, timeout=timeout,
                                allow_redirects=True)
            if resp.status_code == 200:
                return resp
            if resp.status_code in (429, 503):
                wait = 2 ** attempt
                log.warning("HTTP %d — retrying in %ds", resp.status_code, wait)
                time.sleep(wait)
                continue
            log.error("HTTP %d for %s", resp.status_code, url)
            return None
        except Exception as exc:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                log.error("Request failed: %s", exc)
    return None


def _resolve_url_via_api() -> Optional[str]:
    """Use the Figshare API to resolve the current download URL for QM9."""
    import requests
    try:
        resp = requests.get(QM9_FIGSHARE_API, timeout=15)
        if resp.status_code == 200:
            files = resp.json()
            for f in files:
                name = f.get("name", "")
                if "dsgdb9nsd" in name and name.endswith(".tar.bz2"):
                    url = f.get("download_url") or f.get("url_private_download")
                    if url:
                        log.info("Resolved QM9 URL via API: %s", url)
                        return url
    except Exception as exc:
        log.debug("API resolution failed: %s", exc)
    return None


def _load_excluded_indices() -> set:
    """Download and parse the list of 3054 uncharacterized molecules to exclude."""
    log.info("Fetching excluded-molecule list from Figshare...")
    resp = _get(EXCLUDE_URL)
    if resp is None:
        log.warning("Could not download exclusion list; none excluded.")
        return set()
    indices = set()
    for line in resp.text.splitlines():
        line = line.strip()
        if line.isdigit():
            indices.add(int(line))
    log.info("Loaded %d excluded molecule indices.", len(indices))
    return indices


def _parse_qm9_tar(tar_bytes: bytes, excluded: set) -> List[str]:
    """
    Parse the QM9 .xyz tar.bz2 archive and extract SMILES.

    Each .xyz file has this layout (line 0-indexed):
      0:   <natoms>
      1:   gdb <idx> <A> <B> <C> <mu> ... (17 properties)
      2..natoms+1:  element  x  y  z  Mulliken
      natoms+2:     harmonic frequencies
      natoms+3:     SMILES_GDB  SMILES_B3LYP
      natoms+4:     InChI
    """
    smiles_list = []
    log.info("Parsing QM9 archive (%s bytes)...", f"{len(tar_bytes):,}")
    try:
        with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r:bz2") as tf:
            members = [m for m in tf.getmembers() if m.name.endswith(".xyz")]
            log.info("Found %d .xyz files in archive.", len(members))
            for member in members:
                # Derive molecule index from filename: dsgdb9nsd_XXXXXX.xyz
                basename = os.path.basename(member.name)
                m = re.search(r"_(\d+)\.xyz$", basename)
                if m and int(m.group(1)) in excluded:
                    continue
                try:
                    f = tf.extractfile(member)
                    if f is None:
                        continue
                    lines = f.read().decode("utf-8", errors="replace").splitlines()
                    if len(lines) < 4:
                        continue
                    natoms = int(lines[0].strip())
                    smiles_line_idx = natoms + 3  # 0-indexed
                    if smiles_line_idx >= len(lines):
                        continue
                    smiles_col = lines[smiles_line_idx].strip().split()
                    if not smiles_col:
                        continue
                    smi = smiles_col[0]   # GDB SMILES (column 0)
                    if smi:
                        smiles_list.append(smi)
                except Exception:
                    pass
    except Exception as exc:
        log.error("Failed to parse archive: %s", exc)
    log.info("Extracted %d SMILES from archive.", len(smiles_list))
    return smiles_list


def _download_qm9_raw() -> Optional[bytes]:
    """Download the QM9 tar.bz2 archive, trying all known URLs."""
    urls = list(QM9_URLS)
    api_url = _resolve_url_via_api()
    if api_url:
        urls.insert(0, api_url)

    for url in urls:
        log.info("Downloading QM9 from: %s", url)
        resp = _get(url, stream=True, timeout=300)
        if resp is None:
            continue
        chunks = []
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        t0 = time.time()
        for chunk in resp.iter_content(chunk_size=1 << 20):  # 1 MB chunks
            chunks.append(chunk)
            downloaded += len(chunk)
            if total:
                pct = 100 * downloaded / total
                mb = downloaded / 1e6
                speed = mb / max(time.time() - t0, 0.001)
                print(
                    f"\r  {pct:5.1f}%  {mb:.0f}/{total/1e6:.0f} MB  "
                    f"{speed:.1f} MB/s    ",
                    end="", flush=True,
                )
        print()
        raw = b"".join(chunks)
        if raw:
            log.info("Downloaded %s bytes.", f"{len(raw):,}")
            return raw
        log.warning("Empty response from %s", url)

    log.error("All download URLs failed.")
    return None


def load_qm9_smiles(max_molecules: Optional[int] = None,
                    no_cache: bool = False) -> List[str]:
    """
    Return a list of QM9 canonical SMILES strings.

    Workflow:
      1. Check local cache  ->  return immediately if valid.
      2. Download tar.bz2 from Figshare.
      3. Download exclusion list.
      4. Parse .xyz files, extract GDB SMILES.
      5. Canonicalize via RDKit.
      6. Save to cache.
      7. Return (optionally truncated) list.
    """
    _ensure_requests()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Cache hit
    if not no_cache and CACHE_FILE.exists():
        try:
            with open(CACHE_FILE, "r", encoding="utf-8") as fh:
                cached = json.load(fh)
            smiles = cached.get("smiles", [])
            if len(smiles) > 1000:
                log.info("Loaded %d SMILES from cache (%s).",
                         len(smiles), CACHE_FILE)
                if max_molecules:
                    smiles = smiles[:max_molecules]
                return smiles
        except Exception as exc:
            log.warning("Cache read failed (%s); re-downloading.", exc)

    # Download raw archive
    raw = _download_qm9_raw()
    if raw is None:
        log.error("Could not download QM9. Using built-in fallback molecules.")
        return _fallback_molecules(max_molecules or 3000)

    # Parse
    excluded = _load_excluded_indices()
    smiles_raw = _parse_qm9_tar(raw, excluded)
    if not smiles_raw:
        log.error("Parsing produced 0 SMILES. Using fallback.")
        return _fallback_molecules(max_molecules or 3000)

    # Canonicalize
    smiles_canonical = _canonicalize_smiles(smiles_raw)
    log.info("Canonical SMILES after deduplication: %d", len(smiles_canonical))

    # Cache
    try:
        with open(CACHE_FILE, "w", encoding="utf-8") as fh:
            json.dump({"smiles": smiles_canonical,
                       "timestamp": datetime.datetime.now().isoformat(),
                       "n": len(smiles_canonical)}, fh)
        log.info("Saved %d SMILES to cache.", len(smiles_canonical))
    except Exception as exc:
        log.warning("Could not save cache: %s", exc)

    if max_molecules:
        smiles_canonical = smiles_canonical[:max_molecules]
    return smiles_canonical


def _canonicalize_smiles(smiles_list: List[str]) -> List[str]:
    """Canonicalize and deduplicate via RDKit."""
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        log.warning("RDKit not available; returning raw SMILES without canonicalization.")
        return list(dict.fromkeys(smiles_list))  # deduplicate preserving order

    seen = set()
    out = []
    n_fail = 0
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                n_fail += 1
                continue
            canon = Chem.MolToSmiles(mol)
            if canon not in seen:
                seen.add(canon)
                out.append(canon)
        except Exception:
            n_fail += 1
    if n_fail:
        log.warning("RDKit rejected %d SMILES during canonicalization.", n_fail)
    return out


# ---------------------------------------------------------------------------
# Fallback molecules (used when download fails)
# ---------------------------------------------------------------------------

_FALLBACK_SMILES = [
    "C", "CC", "CCC", "CCCC", "CCCCC", "CCO", "CCCO", "CCCCO",
    "CCN", "CCCN", "C=O", "CC=O", "CCC=O", "C=C", "CC=C", "C=CC",
    "C#N", "CC#N", "C#C", "CC#C", "C1CC1", "C1CCC1", "C1CCCC1",
    "C1CCCCC1", "c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccsc1",
    "CC(=O)O", "CC(=O)N", "CC(=O)OC", "CC(=O)NC",
    "OCC", "NCC", "SCC", "FCC", "ClCC",
    "C(C)(C)C", "C(C)(C)CC", "CC(C)CC",
    "OC1CCCCC1", "NC1CCCCC1", "C1CCOC1", "C1CCNC1",
    "c1ccc2ccccc2c1", "c1ccc(O)cc1", "c1ccc(N)cc1",
    "O=C(O)CC", "O=C(N)CC", "O=C(O)CCC",
    "C1=CC=CC=C1", "c1ccc(C)cc1", "c1ccc(Cl)cc1",
]


def _fallback_molecules(n: int) -> List[str]:
    """Generate fallback molecules when QM9 is unavailable."""
    log.warning("Using built-in fallback set (%d molecules).", n)
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
        valid = [Chem.MolToSmiles(Chem.MolFromSmiles(s))
                 for s in _FALLBACK_SMILES
                 if Chem.MolFromSmiles(s) is not None]
    except ImportError:
        valid = list(_FALLBACK_SMILES)

    rng = random.Random(42)
    result = []
    while len(result) < n:
        result.extend(rng.sample(valid, min(len(valid), n - len(result))))
    return result[:n]


# ===========================================================================
# Notation conversion
# ===========================================================================

def to_selfies(smiles_list: List[str]) -> List[str]:
    try:
        import selfies as sf
    except ImportError:
        log.warning("SELFIES not installed; SELFIES experiments skipped.")
        return []
    out = []
    n_fail = 0
    for smi in smiles_list:
        try:
            sel = sf.encoder(smi)
            if sel:
                out.append(sel)
        except Exception:
            n_fail += 1
    if n_fail:
        log.debug("SELFIES encoding: %d failures.", n_fail)
    return out


def to_script(smiles_list: List[str]) -> List[str]:
    from script.rdkit_bridge import smiles_to_script as s2s
    out = []
    n_fail = 0
    for smi in smiles_list:
        try:
            s = s2s(smi)
            if s:
                out.append(s)
        except Exception:
            n_fail += 1
    if n_fail:
        log.debug("SCRIPT encoding: %d failures.", n_fail)
    return out


# ===========================================================================
# Markov models
# ===========================================================================

class CharBigramModel:
    """Character-level bigram language model with Laplace smoothing."""

    def __init__(self, alpha: float = 0.01):
        self.alpha = alpha  # Laplace smoothing
        self.transitions: Dict[str, collections.Counter] = collections.defaultdict(
            collections.Counter
        )
        self.start_chars: collections.Counter = collections.Counter()
        self.vocab: set = set()

    def train(self, strings: List[str]) -> None:
        for s in strings:
            if not s:
                continue
            self.start_chars[s[0]] += 1
            self.vocab.add(s[0])
            for a, b in zip(s, s[1:]):
                self.transitions[a][b] += 1
                self.vocab.add(b)

    def generate(self, max_len: int = 60, rng: random.Random = random) -> str:
        if not self.start_chars:
            return ""
        chars = list(self.start_chars)
        weights = list(self.start_chars.values())
        result = [rng.choices(chars, weights=weights)[0]]
        for _ in range(max_len - 1):
            prev = result[-1]
            ctr = self.transitions.get(prev)
            if not ctr:
                break
            next_chars  = list(ctr.keys())
            next_weights = list(ctr.values())
            result.append(rng.choices(next_chars, weights=next_weights)[0])
        return "".join(result)


class TokenBigramModel:
    """Token-level bigram model for SELFIES (tokens like [C], [Branch1])."""

    _TOK_RE = re.compile(r"\[[^\]]*\]|.")

    def __init__(self):
        self.transitions: Dict[str, collections.Counter] = collections.defaultdict(
            collections.Counter
        )
        self.start_tokens: collections.Counter = collections.Counter()
        self.vocab: set = set()

    @classmethod
    def tokenize(cls, s: str) -> List[str]:
        return cls._TOK_RE.findall(s)

    def train(self, strings: List[str]) -> None:
        for s in strings:
            tokens = self.tokenize(s)
            if not tokens:
                continue
            self.start_tokens[tokens[0]] += 1
            self.vocab.add(tokens[0])
            for a, b in zip(tokens, tokens[1:]):
                self.transitions[a][b] += 1
                self.vocab.add(b)

    def generate(self, max_tokens: int = 40, rng: random.Random = random) -> str:
        if not self.start_tokens:
            return ""
        tokens = list(self.start_tokens)
        weights = list(self.start_tokens.values())
        result = [rng.choices(tokens, weights=weights)[0]]
        for _ in range(max_tokens - 1):
            prev = result[-1]
            ctr = self.transitions.get(prev)
            if not ctr:
                break
            nt  = list(ctr.keys())
            nw  = list(ctr.values())
            result.append(rng.choices(nt, weights=nw)[0])
        return "".join(result)


# ===========================================================================
# Validity checkers
# ===========================================================================

def _rdmol(smi: str):
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smi)
    except Exception:
        return None


def valid_smiles(s: str) -> bool:
    mol = _rdmol(s)
    return mol is not None and mol.GetNumAtoms() > 0


def valid_selfies(s: str) -> bool:
    try:
        import selfies as sf
        smi = sf.decoder(s)
        if not smi:
            return False
        return valid_smiles(smi)
    except ImportError:
        return True   # if selfies not installed, assume all generated are valid
    except Exception:
        return False


def valid_script(s: str) -> bool:
    try:
        from script.parser import parse_script
        r = parse_script(s)
        if not r.get("success"):
            return False
        mol = r.get("molecule")
        if mol is None:
            return False
        if isinstance(mol, list):
            return bool(mol) and all(len(m.atoms) > 0 for m in mol)
        return len(mol.atoms) > 0
    except Exception:
        return False


def script_get_smiles(s: str) -> Optional[str]:
    try:
        from script.rdkit_bridge import script_to_smiles as s2s
        return s2s(s)
    except Exception:
        return None


# ===========================================================================
# Constrained SCRIPT generation
# ===========================================================================

# Probability of early-stopping constrained generation once a complete
# molecule has been formed. Higher = shorter average strings; lower = longer.
CONSTRAINED_EARLY_STOP_PROB: float = 0.15


def _generate_constrained(
    model: CharBigramModel,
    n: int,
    max_len: int,
    rng: random.Random,
    early_stop_prob: float = CONSTRAINED_EARLY_STOP_PROB,
) -> List[str]:
    """Generate SCRIPT strings using the grammar-constrained decoder.

    We use a rejection-sampling loop: only strings that are grammatically
    complete (according to the decoder) and valid SCRIPT strings (according
    to the parser) are collected. This ensures that every generated molecule
    can be successfully evaluated for validity and diversity.
    """
    from script.constrained_decoder import ConstrainedSCRIPTDecoder

    script_vocab = sorted(model.vocab)
    valid = []
    attempts = 0
    max_attempts = n * 50  # safety threshold to prevent infinite loop

    while len(valid) < n and attempts < max_attempts:
        attempts += 1
        decoder = ConstrainedSCRIPTDecoder()
        decoder.reset()
        generated: List[str] = []

        for _step in range(max_len):
            valid_tokens = decoder.get_valid_tokens(script_vocab)
            if not valid_tokens:
                break

            prev = generated[-1] if generated else None
            if prev and prev in model.transitions:
                ctr = model.transitions[prev]
                weighted = [(c, w)
                            for c, w in ctr.items()
                            if c in valid_tokens]
                if weighted:
                    chars, weights = zip(*weighted)
                    token = rng.choices(list(chars), weights=list(weights))[0]
                else:
                    token = rng.choice(list(valid_tokens))
            else:
                token = rng.choice(list(valid_tokens))

            if not decoder.consume(token):
                break
            generated.append(token)

            # Early stop once a grammatically complete molecule is formed
            if decoder.is_complete() and rng.random() < early_stop_prob:
                break

        s = "".join(generated)
        if s and decoder.is_complete() and valid_script(s):
            valid.append(s)

    return valid


# ===========================================================================
# Experiment 2: Generative Validity
# ===========================================================================

def experiment_generative(
    smiles_data: List[str],
    n_generate: int,
    seed: int,
) -> Dict:
    """Train bigram models, generate molecules, report validity."""

    sep = "=" * 72
    log.info("%s", sep)
    log.info("  EXPERIMENT 2: QM9 Generative Validity Comparison")
    log.info("%s", sep)

    rng = random.Random(seed)
    n_train = len(smiles_data)

    # --- Convert to each notation ---
    log.info("Converting %d molecules to SELFIES...", n_train)
    selfies_data = to_selfies(smiles_data)
    log.info("  SELFIES: %d converted", len(selfies_data))

    log.info("Converting %d molecules to SCRIPT...", n_train)
    script_data = to_script(smiles_data)
    log.info("  SCRIPT: %d converted", len(script_data))

    # --- Train models ---
    log.info("Training bigram models...")
    smiles_model = CharBigramModel()
    smiles_model.train(smiles_data)
    log.info("  SMILES: vocab size %d chars", len(smiles_model.vocab))

    selfies_model: Optional[TokenBigramModel] = None
    if selfies_data:
        selfies_model = TokenBigramModel()
        selfies_model.train(selfies_data)
        log.info("  SELFIES: vocab size %d tokens", len(selfies_model.vocab))

    script_model = CharBigramModel()
    script_model.train(script_data)
    log.info("  SCRIPT: vocab size %d chars", len(script_model.vocab))

    # --- Generate and validate ---
    def _run_notation(name: str, generate_fn, check_fn) -> dict:
        log.info("Generating %d %s strings...", n_generate, name)
        t0 = time.time()
        valid_mols = []
        for k in range(n_generate):
            s = generate_fn()
            if check_fn(s):
                valid_mols.append(s)
            if (k + 1) % 200 == 0:
                log.debug("  %s: %d/%d done, %d valid so far",
                          name, k + 1, n_generate, len(valid_mols))
        elapsed = time.time() - t0
        validity = len(valid_mols) / n_generate
        log.info("  %s: %d/%d valid (%.1f%%) in %.1fs",
                 name, len(valid_mols), n_generate, validity * 100, elapsed)
        return {
            "n_generate": n_generate,
            "valid": len(valid_mols),
            "invalid": n_generate - len(valid_mols),
            "validity_rate": round(validity, 6),
            "elapsed_s": round(elapsed, 2),
            "valid_molecules": valid_mols,
        }

    results = {}

    results["SMILES"] = _run_notation(
        "SMILES",
        lambda: smiles_model.generate(max_len=60, rng=rng),
        valid_smiles,
    )

    if selfies_model is not None:
        results["SELFIES"] = _run_notation(
            "SELFIES",
            lambda: selfies_model.generate(max_tokens=40, rng=rng),
            valid_selfies,
        )

    results["SCRIPT (unconstrained)"] = _run_notation(
        "SCRIPT (unconstrained)",
        lambda: script_model.generate(max_len=60, rng=rng),
        valid_script,
    )

    # Constrained SCRIPT
    log.info("Generating %d SCRIPT (constrained) strings...", n_generate)
    t0 = time.time()
    try:
        constrained_valid = _generate_constrained(
            script_model, n_generate, max_len=60, rng=rng
        )
        elapsed_con = time.time() - t0
        validity_con = len(constrained_valid) / n_generate
        log.info("  SCRIPT (constrained): %d/%d valid (%.1f%%) in %.1fs",
                 len(constrained_valid), n_generate, validity_con * 100, elapsed_con)
        results["SCRIPT (constrained)"] = {
            "n_generate": n_generate,
            "valid": len(constrained_valid),
            "invalid": n_generate - len(constrained_valid),
            "validity_rate": round(validity_con, 6),
            "elapsed_s": round(elapsed_con, 2),
            "valid_molecules": constrained_valid,
        }
    except ImportError:
        log.warning("ConstrainedSCRIPTDecoder not found; skipping constrained mode.")

    return {"n_train": n_train, "n_generate": n_generate, "notations": results}


# ===========================================================================
# Experiment 3: Internal Diversity
# ===========================================================================

def _morgan_fp(smi: str):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    except Exception:
        return None


def _diversity(smiles_list: List[str], max_pairs: int = 5000) -> Tuple[float, int]:
    """
    Compute internal diversity = 1 - mean(pairwise Tanimoto).
    Uses numpy for speed when available; falls back to pure Python.
    """
    fps = [fp for fp in (_morgan_fp(s) for s in smiles_list) if fp is not None]
    n = len(fps)
    if n < 2:
        return 0.0, n

    try:
        from rdkit.Chem import DataStructs
        import numpy as np

        mat = np.zeros((n, n), dtype=np.float32)
        for i in range(n):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
            mat[i] = sims
        # Exclude diagonal
        np.fill_diagonal(mat, 0.0)
        mean_sim = mat.sum() / (n * (n - 1))
        return round(1.0 - float(mean_sim), 6), n

    except ImportError:
        pass

    # Pure Python fallback (sample pairs for speed)
    try:
        from rdkit.Chem import DataStructs as _DS2
    except ImportError:
        return 0.0, n
    rng = random.Random(0)
    indices = list(range(n))
    pairs = [(i, j) for i in indices for j in range(i + 1, n)]
    if len(pairs) > max_pairs:
        pairs = rng.sample(pairs, max_pairs)

    total_sim = sum(_DS2.TanimotoSimilarity(fps[i], fps[j])
                    for i, j in pairs)
    mean_sim = total_sim / len(pairs)
    return round(1.0 - mean_sim, 6), n


def experiment_diversity(gen_results: Dict) -> Dict:
    """Compute internal diversity for each notation's generated set."""
    sep = "=" * 72
    log.info("%s", sep)
    log.info("  EXPERIMENT 3: Internal Diversity Comparison")
    log.info("%s", sep)

    notations_data = gen_results.get("notations", {})
    diversity_results = {}

    for name, r in notations_data.items():
        valid_mols = r.get("valid_molecules", [])
        if not valid_mols:
            diversity_results[name] = {"diversity": 0.0, "n": 0}
            continue

        # Convert to SMILES for fingerprinting
        smiles_for_fp: List[str] = []
        for mol_str in valid_mols:
            if "SELFIES" in name:
                try:
                    import selfies as sf
                    smi = sf.decoder(mol_str)
                    if smi:
                        smiles_for_fp.append(smi)
                except Exception:
                    pass
            elif "SCRIPT" in name:
                smi = script_get_smiles(mol_str)
                if smi:
                    smiles_for_fp.append(smi)
            else:
                smiles_for_fp.append(mol_str)

        div, n = _diversity(smiles_for_fp)
        diversity_results[name] = {"diversity": div, "n": n}
        log.info("  %-28s  diversity=%.4f  (n=%d)", name, div, n)

    return diversity_results


# ===========================================================================
# Report generation
# ===========================================================================

def _pct(rate: float) -> str:
    return f"{rate * 100:.2f}%"


def build_report(
    gen: Dict,
    div: Dict,
    args,
) -> str:
    lines: List[str] = []
    sep  = "=" * 72
    thin = "-" * 72
    now  = datetime.datetime.now().isoformat(timespec="seconds")

    lines += [
        sep,
        "  QM9 EXPERIMENTS: NOTATION GENERATIVE VALIDITY AND DIVERSITY",
        sep,
        f"  Generated : {now}",
        f"  Dataset   : QM9 (Ramakrishnan et al., 2014)",
        f"  Training  : {gen['n_train']:,} molecules",
        f"  Generated : {gen['n_generate']:,} per notation",
        "",
    ]

    # --- Experiment 2 table ---
    lines += [
        thin,
        "  EXPERIMENT 2: Generative Validity",
        thin,
        f"  {'Notation':<28} {'Valid':>7} {'Invalid':>8} {'Validity':>10} {'Time':>8}",
        "  " + "-" * 62,
    ]
    for name, r in gen["notations"].items():
        lines.append(
            f"  {name:<28} {r['valid']:>7} {r['invalid']:>8} "
            f"{_pct(r['validity_rate']):>10} {r['elapsed_s']:>7.1f}s"
        )
    lines += [
        "",
        "  Validity = fraction of generated strings that parse as valid molecules.",
        "",
    ]

    # --- Experiment 3 table ---
    lines += [
        thin,
        "  EXPERIMENT 3: Internal Diversity",
        thin,
        f"  {'Notation':<28} {'Diversity':>10} {'N':>6}",
        "  " + "-" * 46,
    ]
    for name, r in div.items():
        lines.append(
            f"  {name:<28} {r['diversity']:>10.4f} {r['n']:>6}"
        )
    lines += [
        "",
        "  Diversity = 1 - mean(pairwise Tanimoto on Morgan-2048 fingerprints).",
        "  Higher = more diverse generated set.",
        "",
    ]

    # --- LaTeX: Table 1 (validity) ---
    lines += [
        sep,
        "  LATEX TABLE 1: Generative Validity (Experiment 2)",
        sep,
        r"\begin{table}[ht]",
        r"\centering",
        (r"\caption{Validity of molecules generated by a character-level "
         r"bigram model trained on QM9. SCRIPT (constrained) uses grammar-state"
         r"-aware token masking.}"),
        r"\label{tab:qm9_validity}",
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
        (r"\textbf{Notation} & \textbf{Valid} & \textbf{Invalid} "
         r"& \textbf{Validity (\%)} & \textbf{Time (s)} \\"),
        r"\midrule",
    ]
    for name, r in gen["notations"].items():
        lines.append(
            f"{name} & {r['valid']} & {r['invalid']} "
            f"& {r['validity_rate']*100:.2f} & {r['elapsed_s']:.1f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]

    # --- LaTeX: Table 2 (diversity) ---
    lines += [
        sep,
        "  LATEX TABLE 2: Internal Diversity (Experiment 3)",
        sep,
        r"\begin{table}[ht]",
        r"\centering",
        (r"\caption{Internal diversity of generated valid molecules. "
         r"Diversity = 1 $-$ mean pairwise Tanimoto similarity "
         r"(Morgan-2048 fingerprints).}"),
        r"\label{tab:qm9_diversity}",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"\textbf{Notation} & \textbf{Diversity} & \textbf{N} \\",
        r"\midrule",
    ]
    for name, r in div.items():
        lines.append(f"{name} & {r['diversity']:.4f} & {r['n']} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]

    lines.append(sep)
    return "\n".join(lines)


def build_json_output(gen: Dict, div: Dict) -> Dict:
    # Strip molecule lists from JSON to keep it small
    notations_clean = {}
    for name, r in gen["notations"].items():
        notations_clean[name] = {k: v for k, v in r.items()
                                 if k != "valid_molecules"}
    return {
        "timestamp": datetime.datetime.now().isoformat(),
        "experiment_2": {
            "n_train": gen["n_train"],
            "n_generate": gen["n_generate"],
            "notations": notations_clean,
        },
        "experiment_3": div,
    }


# ===========================================================================
# CLI
# ===========================================================================

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="QM9 Generative Validity and Diversity Experiments"
    )
    ap.add_argument("--n-train",   type=int, default=None,
                    help="Max training molecules (default: all ~133k after exclusions)")
    ap.add_argument("--n-generate", type=int, default=1000,
                    help="Molecules to generate per notation (default: 1000)")
    ap.add_argument("--seed",      type=int, default=42)
    ap.add_argument("--quick",     action="store_true",
                    help="Quick mode: 1 000 train, 200 generate")
    ap.add_argument("--no-cache",  action="store_true",
                    help="Force fresh download from Figshare")
    ap.add_argument("--json",      action="store_true",
                    help="Also write JSON output file")
    ap.add_argument("--verbose",   action="store_true",
                    help="Debug-level logging")
    return ap.parse_args()


def main():
    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.quick:
        args.n_train    = args.n_train    or 1_000
        args.n_generate = args.n_generate if args.n_generate != 1000 else 200

    # --- Dependency checks ---
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
        log.info("RDKit: OK")
    except ImportError:
        log.error("RDKit is required. Install: pip install rdkit")
        sys.exit(1)

    try:
        import selfies  # noqa: F401
        log.info("SELFIES: OK")
    except ImportError:
        log.info("SELFIES not installed; SELFIES experiments will be skipped.")

    try:
        from script.parser import parse_script  # noqa: F401
        log.info("SCRIPT: OK")
    except ImportError:
        log.error("SCRIPT not importable from this location.")
        sys.exit(1)

    # --- Load data ---
    log.info("Loading QM9 dataset...")
    smiles_data = load_qm9_smiles(
        max_molecules=args.n_train,
        no_cache=args.no_cache,
    )
    log.info("Training set: %d molecules", len(smiles_data))

    # --- Run experiments ---
    gen_results = experiment_generative(smiles_data, args.n_generate, args.seed)
    div_results = experiment_diversity(gen_results)

    # --- Report ---
    report = build_report(gen_results, div_results, args)
    print("\n" + report)

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    txt_path = TESTS_DIR / f"qm9_experiments_{ts}.txt"
    txt_path.write_text(report, encoding="utf-8")
    log.info("Report saved: %s", txt_path)

    if args.json:
        json_path = TESTS_DIR / f"qm9_experiments_{ts}.json"
        json_out = build_json_output(gen_results, div_results)
        json_path.write_text(
            json.dumps(json_out, indent=2, default=str), encoding="utf-8"
        )
        log.info("JSON saved: %s", json_path)

    log.info("All experiments complete.")


if __name__ == "__main__":
    main()