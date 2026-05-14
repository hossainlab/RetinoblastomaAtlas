"""
data_exploration.py
===================
Data extraction, format detection, sample loading, and integrity checks.

Functions
---------
Extraction
    inspect_tar              — peek inside tar without extracting
    extract_tar              — universal GEO tar extractor (all 5 structures)
    extract_tar_if_needed    — skip if already extracted
    extract_multiple_tars    — merge multiple tars into one folder

Format & discovery
    detect_format            — detect data format in a folder
    discover_samples         — find all sample folders under a root dir

Loading
    load_sample              — load one sample (any format)
    load_all_samples         — load + concatenate all samples

Data integrity
    check_raw_counts         — verify data is raw integer counts
    detect_species           — infer human vs mouse from gene names
    convert_ensembl_to_symbols — ENSG IDs → gene symbols (requires mygene)
    handle_duplicate_genes   — resolve duplicate var_names
    inject_metadata          — join sample metadata CSV into adata.obs
    load_if_not_exists       — resume from existing h5ad
    diagnose_sample          — debug a failed sample folder

Example
-------
    from utils.data_exploration import (
        inspect_tar, extract_tar_if_needed,
        load_all_samples, check_raw_counts, detect_species,
    )
"""

import os
import re
import shutil
import tarfile
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0


# ─────────────────────────────────────────────────────────────────────────────
# PRIVATE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _strip_tar_ext(name: str) -> str:
    """Strip .tar.gz / .tar.bz2 / .tar.xz / .tar from a filename."""
    for ext in (".tar.gz", ".tar.bz2", ".tar.xz", ".tar"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _clean_gsm_name(archive_name: str) -> str:
    """Strip tar extension + trailing GEO count-file suffixes from archive name.

    Example: GSM5139852_RB01_rep1_counts.mtx.tsv.tar.gz → GSM5139852_RB01_rep1
    """
    name = _strip_tar_ext(archive_name)
    name = re.sub(
        r'[._](counts?|matrix|barcodes?|features?|genes?|raw|filtered)'
        r'([._].*)?$',
        '', name, flags=re.IGNORECASE,
    )
    return name or _strip_tar_ext(archive_name)


def _flatten_sample_dir(sample_dir: Path) -> None:
    """Move data files from deeply nested subdirs up to sample_dir root.

    Handles inner tars that stored absolute/user paths
    (e.g. home/user/Aspera/RB_raw/sample_name/).
    """
    _DATA_EXTS = {".gz", ".h5", ".h5ad", ".csv", ".tsv", ".txt", ".loom", ".mtx"}

    # Already flat — data files present at root
    if any(f.is_file() and f.suffix in _DATA_EXTS for f in sample_dir.iterdir()):
        return

    # Find deepest directory that contains data files
    data_dir = None
    for d in sorted(sample_dir.rglob("*"), key=lambda p: len(p.parts), reverse=True):
        if not d.is_dir():
            continue
        if any(f.is_file() and f.suffix in _DATA_EXTS for f in d.iterdir()):
            data_dir = d
            break

    if data_dir is None:
        return

    # Move files to sample_dir root
    for f in list(data_dir.iterdir()):
        if f.is_file():
            target = sample_dir / f.name
            if not target.exists():
                shutil.move(str(f), str(target))

    # Remove intermediate subdirectory tree
    for child in list(sample_dir.iterdir()):
        if child.is_dir():
            shutil.rmtree(str(child))


def _load_mtx(folder: str, sample_id: str) -> ad.AnnData:
    folder   = Path(folder)
    genes    = folder / "genes.tsv"
    features = folder / "features.tsv"
    if genes.exists() and not features.exists():
        shutil.copy(str(genes), str(features))
        print(f"  [{sample_id}] Renamed genes.tsv → features.tsv (CellRanger v2)")
    return sc.read_10x_mtx(str(folder), var_names="gene_symbols", cache=False, gex_only=True)


def _load_h5(folder: str, sample_id: str) -> ad.AnnData:
    files    = list(Path(folder).glob("*.h5"))
    filtered = [f for f in files if "filtered" in f.name.lower()]
    target   = filtered[0] if filtered else files[0]
    print(f"  [{sample_id}] h5 file: {target.name}")
    return sc.read_10x_h5(str(target))


def _load_h5ad(folder: str, sample_id: str) -> ad.AnnData:
    files = list(Path(folder).glob("*.h5ad"))
    adata = sc.read_h5ad(str(files[0]))
    if adata.X is not None:
        sample_vals = (adata.X[:5, :5].toarray()
                       if sparse.issparse(adata.X) else adata.X[:5, :5])
        if np.any((sample_vals > 0) & (sample_vals < 1)):
            print(f"  [{sample_id}] WARNING: non-integer values — may be pre-normalized")
    return adata


def _load_csv(folder: str, sample_id: str) -> ad.AnnData:
    files = (list(Path(folder).glob("*.csv.gz")) + list(Path(folder).glob("*.csv")) +
             list(Path(folder).glob("*.tsv.gz")) + list(Path(folder).glob("*.tsv")))
    sep = "\t" if "tsv" in files[0].name else ","
    df  = pd.read_csv(str(files[0]), index_col=0, sep=sep)
    if df.shape[0] > df.shape[1]:
        print(f"  [{sample_id}] Transposing: genes×cells → cells×genes")
        df = df.T
    return ad.AnnData(X=sparse.csr_matrix(df.values.astype(np.float32)),
                      obs=pd.DataFrame(index=df.index),
                      var=pd.DataFrame(index=df.columns))


def _load_txt(folder: str, sample_id: str) -> ad.AnnData:
    files = (list(Path(folder).glob("*.tsv.gz")) + list(Path(folder).glob("*.tsv")) +
             list(Path(folder).glob("*.txt.gz")) + list(Path(folder).glob("*.txt")))
    if not files:
        raise FileNotFoundError(f"No tsv/txt file found in {folder}")
    print(f"  [{sample_id}] Loading TSV: {files[0].name}")
    df = pd.read_csv(str(files[0]), sep="\t", index_col=0)
    if df.shape[0] > df.shape[1]:
        print(f"  [{sample_id}] Transposing: genes×cells → cells×genes")
        df = df.T
    return ad.AnnData(X=sparse.csr_matrix(df.values.astype(np.float32)),
                      obs=pd.DataFrame(index=df.index),
                      var=pd.DataFrame(index=df.columns))


def _load_loom(folder: str, sample_id: str) -> ad.AnnData:
    files = list(Path(folder).glob("*.loom"))
    return sc.read_loom(str(files[0]))


_LOADERS = {
    "mtx":  _load_mtx,
    "h5":   _load_h5,
    "h5ad": _load_h5ad,
    "csv":  _load_csv,
    "txt":  _load_txt,
    "loom": _load_loom,
}


# ─────────────────────────────────────────────────────────────────────────────
# EXTRACTION
# ─────────────────────────────────────────────────────────────────────────────

def inspect_tar(tar_path: str, n: int = 20) -> str:
    """
    Peek inside a tar archive without extracting.

    Detects GEO structure type so you know what extract_tar will do.
    Always run this before extract_tar.

    Parameters
    ----------
    tar_path : str   Path to .tar file.
    n        : int   Entries to preview (default 20).

    Returns
    -------
    str  'nested_archives' | 'subfolders' | 'flat_gsm_prefix' |
         'single_sample'  | 'unknown'

    Example
    -------
    >>> inspect_tar("../data/GSE168434_RAW.tar")
    # [inspect] Structure A — NESTED ARCHIVES: 20 per-sample .tar.gz
    """
    _DATA_EXTS = {".mtx", ".h5", ".h5ad", ".csv", ".tsv", ".txt", ".loom"}

    with tarfile.open(tar_path) as tar:
        members = tar.getmembers()

    print(f"[inspect] {Path(tar_path).name}  —  {len(members)} entries total")
    print(f"[inspect] First {min(n, len(members))} entries:")
    for m in members[:n]:
        tag = "DIR " if m.isdir() else f"{m.size/1024:>8.1f} KB"
        print(f"  {tag}  {m.name}")
    if len(members) > n:
        print(f"  ... ({len(members) - n} more)")

    names         = [m.name for m in members]
    is_archive    = lambda x: any(x.endswith(e) for e in (".tar.gz", ".tar.bz2", ".tar.xz"))
    gsm_match     = lambda x: re.match(r"^(GSM\d+)", Path(x).name)
    archive_files = [x for x in names if is_archive(x)]
    subdir_files  = [x for x in names if "/" in x]
    data_files    = [x for x in names if Path(x).suffix in _DATA_EXTS or
                     any(x.endswith(e) for e in (".mtx.gz", ".tsv.gz", ".txt.gz"))]
    gsm_files     = [x for x in names if gsm_match(x)]

    if archive_files:
        structure = "nested_archives"
        print(f"\n[inspect] Structure A — NESTED ARCHIVES: {len(archive_files)} per-sample .tar.gz")
        print(f"  → extract each inner archive to its own folder")
    elif subdir_files:
        top_dirs  = sorted({x.split("/")[0] for x in subdir_files})
        structure = "subfolders"
        print(f"\n[inspect] Structure B — SUBFOLDERS: {len(top_dirs)} top-level dirs")
        for d in top_dirs[:5]: print(f"  {d}/")
    elif gsm_files:
        gsm_ids   = sorted({gsm_match(x).group(1) for x in gsm_files if gsm_match(x)})
        structure = "flat_gsm_prefix"
        print(f"\n[inspect] Structure C — FLAT GSM PREFIX: {len(gsm_ids)} unique GSM IDs")
        for g in gsm_ids[:5]: print(f"  {g}")
    elif data_files:
        structure = "single_sample"
        print(f"\n[inspect] Structure D — SINGLE SAMPLE: {len(data_files)} data files at root")
    else:
        structure = "unknown"
        print(f"\n[inspect] Structure E — UNKNOWN: manual inspection recommended")

    return structure


def extract_tar(tar_path: str, out_dir: str) -> str:
    """
    Universal GEO tar extractor. Handles all 5 common GEO structures.

    Structures
    ----------
    A) Nested archives  tar → GSM*.tar.gz → each to its own folder
    B) Subfolders       tar → sample_dir/files → unpack directly
    C) Flat GSM prefix  tar → GSM001_f1, GSM001_f2... → group by GSM ID
    D) Single sample    tar → data files at root → one folder
    E) Unknown          flat unpack with warning

    Run inspect_tar() first.

    Parameters
    ----------
    tar_path : str   Path to outer .tar file.
    out_dir  : str   Parent directory (created if absent).

    Returns
    -------
    str  Path to folder containing all sample subfolders.

    Example
    -------
    >>> inspect_tar("../data/GSE168434_RAW.tar")
    >>> data_dir = extract_tar("../data/GSE168434_RAW.tar", "../data/raw/")
    # → ../data/raw/GSE168434_RAW/
    #     ├── GSM5139852_RB01_rep1/
    #     └── GSM5139853_RB01_rep2/
    """
    tar_path = Path(tar_path)
    dest     = Path(out_dir) / _strip_tar_ext(tar_path.name)
    dest.mkdir(parents=True, exist_ok=True)
    print(f"[extract] {tar_path.name} → {dest}")

    with tarfile.open(str(tar_path)) as tar:
        members = tar.getmembers()
    names      = [m.name for m in members]
    is_archive = lambda x: any(x.endswith(e) for e in (".tar.gz", ".tar.bz2", ".tar.xz"))
    gsm_match  = lambda x: re.match(r"^(GSM\d+)", Path(x).name)

    archive_files = [x for x in names if is_archive(x)]
    subdir_files  = [x for x in names if "/" in x]
    gsm_files     = [x for x in names if gsm_match(x)]

    if archive_files:
        print(f"[extract] Structure A: {len(archive_files)} nested archives")
        with tarfile.open(str(tar_path)) as tar:
            tar.extractall(str(dest))
        for archive in (sorted(dest.glob("*.tar.gz")) +
                        sorted(dest.glob("*.tar.bz2")) +
                        sorted(dest.glob("*.tar.xz"))):
            sample_dir = dest / _clean_gsm_name(archive.name)
            sample_dir.mkdir(exist_ok=True)
            with tarfile.open(str(archive)) as inner:
                inner.extractall(str(sample_dir))
            archive.unlink()
            _flatten_sample_dir(sample_dir)
            print(f"  ✓ {sample_dir.name}")

    elif subdir_files:
        top_dirs = sorted({x.split("/")[0] for x in subdir_files})
        print(f"[extract] Structure B: {len(top_dirs)} subfolders")
        with tarfile.open(str(tar_path)) as tar:
            tar.extractall(str(dest))
        for d in top_dirs: print(f"  ✓ {d}/")

    elif gsm_files:
        gsm_ids = sorted({gsm_match(x).group(1) for x in gsm_files if gsm_match(x)})
        print(f"[extract] Structure C: flat GSM files → {len(gsm_ids)} samples")
        with tarfile.open(str(tar_path)) as tar:
            tar.extractall(str(dest))
        gsm_groups: dict = {}
        for f in dest.iterdir():
            if not f.is_file(): continue
            m = gsm_match(f.name)
            if m:
                gsm_groups.setdefault(m.group(1), []).append(f)
        for gsm_id, files in sorted(gsm_groups.items()):
            sample_dir = dest / gsm_id
            sample_dir.mkdir(exist_ok=True)
            for f in files:
                f.rename(sample_dir / f.name)
            print(f"  ✓ {gsm_id}  ({len(files)} files)")

    else:
        print(f"[extract] Structure D: flat/single-sample")
        with tarfile.open(str(tar_path)) as tar:
            tar.extractall(str(dest))
        print(f"  WARNING: no GSM prefix — treated as single sample")

    n_files = sum(1 for f in dest.rglob("*") if f.is_file())
    print(f"\n[extract] Done — {n_files} files in {dest}")
    return str(dest)


def extract_tar_if_needed(tar_path: str, out_dir: str) -> str:
    """
    Extract only if destination folder doesn't exist or is empty.

    Use when re-running a notebook — avoids re-extracting large archives.

    Parameters
    ----------
    tar_path : str   Path to .tar file.
    out_dir  : str   Parent output directory.

    Returns
    -------
    str  Path to extraction folder.

    Example
    -------
    >>> data_dir = extract_tar_if_needed("../data/GSE168434_RAW.tar", "../data/raw/")
    # Second run: [extract] Already extracted → ... (skipping)
    """
    dest = Path(out_dir) / _strip_tar_ext(Path(tar_path).name)
    if dest.exists() and any(dest.iterdir()):
        n_files = sum(1 for f in dest.rglob("*") if f.is_file())
        print(f"[extract] Already extracted → {dest}  ({n_files} files, skipping)")
        return str(dest)
    return extract_tar(tar_path, out_dir)


def extract_multiple_tars(tar_paths: list, out_dir: str) -> str:
    """
    Extract multiple tar files from one study into a single combined folder.

    All samples land together for unified loading with load_all_samples().

    Parameters
    ----------
    tar_paths : list[str]   List of .tar file paths.
    out_dir   : str         Parent output directory.

    Returns
    -------
    str  Path to combined folder.

    Example
    -------
    >>> data_dir = extract_multiple_tars(
    ...     ["../data/GSE001_batch1.tar", "../data/GSE001_batch2.tar"],
    ...     "../data/raw/",
    ... )
    """
    if not tar_paths:
        raise ValueError("tar_paths is empty")

    combined = Path(out_dir) / _strip_tar_ext(Path(tar_paths[0]).name)
    combined.mkdir(parents=True, exist_ok=True)

    for tar_path in tar_paths:
        print(f"\n[extract] Processing {Path(tar_path).name}...")
        tmp_parent = Path(out_dir)
        tmp_name   = _strip_tar_ext(Path(tar_path).name)
        tmp        = tmp_parent / tmp_name

        extract_tar(tar_path, str(tmp_parent))

        src = tmp if tmp.exists() else combined
        if src != combined:
            for item in src.iterdir():
                target = combined / item.name
                if not target.exists():
                    item.rename(target)
                else:
                    print(f"  WARNING: {item.name} already exists — skipping")
            if not any(src.iterdir()):
                src.rmdir()

    n_files = sum(1 for f in combined.rglob("*") if f.is_file())
    print(f"\n[extract] Combined: {n_files} files in {combined}")
    return str(combined)


# ─────────────────────────────────────────────────────────────────────────────
# FORMAT DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def detect_format(folder: str) -> str:
    """
    Detect the single-cell data format inside a folder.

    Priority order: h5ad > h5 > mtx > loom > csv > tsv/txt > unknown

    Parameters
    ----------
    folder : str   Path to the sample folder.

    Returns
    -------
    str  'h5ad' | 'h5' | 'mtx' | 'loom' | 'csv' | 'txt' | 'unknown'

    Example
    -------
    >>> detect_format("./raw/GSM1234567")
    'mtx'
    >>> detect_format("./raw/GSM9999999_unknown_sample")
    'unknown'
    """
    files = [f.lower() for f in os.listdir(folder)]
    if any(f.endswith(".h5ad") for f in files):                          return "h5ad"
    if any(f.endswith(".h5")   for f in files):                          return "h5"
    if any(f in ("matrix.mtx.gz", "matrix.mtx") or
           f.endswith(".mtx.gz") or f.endswith(".mtx")
           for f in files):                                               return "mtx"
    if any(f.endswith(".loom")  for f in files):                         return "loom"
    if any(f.endswith(".csv.gz") or f.endswith(".csv") for f in files):  return "csv"
    if any(f.endswith(".tsv.gz") or f.endswith(".tsv") for f in files):  return "txt"
    if any(f.endswith(".txt.gz") or f.endswith(".txt") for f in files):  return "txt"
    return "unknown"


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE DISCOVERY + LOADING
# ─────────────────────────────────────────────────────────────────────────────

def discover_samples(root_dir: str) -> dict:
    """
    Walk directory tree, find all leaf folders containing data files.

    Parameters
    ----------
    root_dir : str   Root of the extracted archive.

    Returns
    -------
    dict  {sample_id: folder_path}

    Example
    -------
    >>> samples = discover_samples("../data/raw/GSE168434_RAW")
    >>> list(samples.keys())[:3]
    ['GSM5139852_RB01_rep1', 'GSM5139853_RB01_rep2', ...]
    """
    sample_dirs: dict = {}
    root = Path(root_dir)

    for folder in sorted(root.rglob("*")):
        if not folder.is_dir() or folder == root:
            continue
        has_data = any(
            f.suffix in {".gz", ".h5", ".h5ad", ".csv", ".tsv", ".txt", ".loom", ".mtx"}
            for f in folder.iterdir() if f.is_file()
        )
        if has_data:
            sample_dirs[folder.name] = str(folder)

    print(f"[discover] Found {len(sample_dirs)} sample(s)")
    for sid, path in sample_dirs.items():
        print(f"  {sid} → {path}")
    return sample_dirs


def load_sample(folder: str, sample_id: str) -> ad.AnnData | None:
    """
    Load a single sample from disk.

    Handles: format detection, CSR conversion, barcode uniquification,
    raw count preservation in layers['counts'].

    Parameters
    ----------
    folder    : str   Path to the sample folder.
    sample_id : str   Injected into adata.obs['sample'] and barcode prefix.

    Returns
    -------
    AnnData | None  None if format unknown or loading failed.

    Example
    -------
    >>> adata = load_sample("./raw/GSM5139852_RB01_rep1", "GSM5139852")
    # [GSM5139852] Format: txt
    # [GSM5139852] 3842 cells × 33538 genes | sparsity 94.2%
    """
    fmt = detect_format(folder)
    if fmt == "unknown":
        print(f"  [{sample_id}] SKIP — unknown format. Files: {os.listdir(folder)}")
        return None

    print(f"  [{sample_id}] Format: {fmt}")
    try:
        adata = _LOADERS[fmt](folder, sample_id)
    except Exception as exc:
        print(f"  [{sample_id}] FAILED: {exc}")
        return None

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs["sample"] = sample_id
    adata.obs_names     = [f"{sample_id}_{bc}" for bc in adata.obs_names]

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    elif adata.X.format != "csr":
        adata.X = adata.X.tocsr()

    adata.layers["counts"] = adata.X.copy()

    sparsity = 1 - adata.X.nnz / (adata.n_obs * adata.n_vars)
    print(f"  [{sample_id}] {adata.n_obs} cells × {adata.n_vars} genes | sparsity {sparsity:.1%}")
    return adata


def load_all_samples(root_dir: str) -> ad.AnnData:
    """
    Discover, load, and concatenate all samples under root_dir.

    Parameters
    ----------
    root_dir : str   Root directory of extracted data.

    Returns
    -------
    AnnData  Combined object. adata.obs['sample'] identifies each cell's origin.
             adata.layers['counts'] preserves raw integer counts.

    Example
    -------
    >>> adata = load_all_samples("../data/raw/GSE168434_RAW")
    >>> adata.obs['sample'].value_counts()
    """
    sample_dirs  = discover_samples(root_dir)
    adatas: dict = {}
    failed: list = []

    print("\n[load] Loading samples...")
    for sid, path in sample_dirs.items():
        obj = load_sample(path, sid)
        if obj is not None:
            adatas[sid] = obj
        else:
            failed.append(sid)

    if not adatas:
        raise RuntimeError("No samples loaded. Check directory structure.")

    print(f"\n[load] Loaded: {len(adatas)} | Failed: {len(failed)}")
    if failed:
        print(f"  Failed: {failed}")

    combined = ad.concat(adatas.values(), join="outer", fill_value=0)
    combined.obs_names_make_unique()

    if not sparse.issparse(combined.X):
        combined.X = sparse.csr_matrix(combined.X)
    combined.X = combined.X.tocsr()

    print(f"\n[load] Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")
    print(combined.obs["sample"].value_counts().to_string())
    return combined


# ─────────────────────────────────────────────────────────────────────────────
# DATA INTEGRITY
# ─────────────────────────────────────────────────────────────────────────────

def check_raw_counts(adata: ad.AnnData) -> bool:
    """
    Verify adata.X contains raw integer counts, not normalized values.

    Pre-normalized data produces incorrect QC metrics and MAD thresholds.
    Run immediately after loading, before any QC.

    Tests
    -----
    1. All values non-negative
    2. All values integer (within 0.001 tolerance)
    3. Max value < 1,000,000 (plausible for raw counts)
    4. >30% zeros (expected sparse structure)

    Returns
    -------
    bool  True = likely raw counts, False = likely normalized.

    Example
    -------
    >>> assert check_raw_counts(adata), "Data is not raw counts — stop pipeline"
    """
    X = adata.X[:200, :200]
    if sparse.issparse(X):
        X = X.toarray()

    all_nonneg    = bool(np.all(X >= 0))
    all_int       = bool(np.allclose(X, np.round(X), atol=1e-3))
    max_plausible = bool(X.max() < 1_000_000)
    has_zeros     = bool((X == 0).mean() > 0.3)
    is_raw        = all_nonneg and all_int and max_plausible and has_zeros

    print("[check_raw] Counts check:")
    print(f"  All non-negative: {all_nonneg}")
    print(f"  All integers:     {all_int}")
    print(f"  Max plausible:    {max_plausible}  (max={X.max():.1f})")
    print(f"  Has zeros:        {has_zeros}  ({(X==0).mean():.1%} zeros in sample)")
    print(f"  → Likely raw:     {is_raw}")
    if not is_raw:
        print("  WARNING: may be pre-normalized — QC metrics will be unreliable")
    return is_raw


def detect_species(adata: ad.AnnData) -> str:
    """
    Infer species from gene naming conventions in adata.var_names.

    Heuristic scores
    ----------------
    Human: uppercase MT- genes, ENSG IDs, RPS/RPL prefix
    Mouse: lowercase mt- genes, ENSMUSG IDs, Rps/Rpl prefix

    Returns
    -------
    str  'human' | 'mouse' | 'unknown'

    Example
    -------
    >>> species = detect_species(adata)
    >>> if species == "human":
    ...     compute_qc_metrics_human(adata)
    """
    genes = adata.var_names

    human_score = (genes.str.startswith("MT-").sum() +
                   genes.str.startswith("ENSG").sum() +
                   genes.str.match(r"^RPS|^RPL").sum())
    mouse_score = (genes.str.startswith("mt-").sum() +
                   genes.str.startswith("ENSMUSG").sum() +
                   genes.str.match(r"^Rps|^Rpl").sum())

    print("[detect_species] Evidence scores:")
    print(f"  Human: {human_score}  (MT- genes, ENSG IDs, RPS/RPL)")
    print(f"  Mouse: {mouse_score}  (mt- genes, ENSMUSG IDs, Rps/Rpl)")

    if   human_score > mouse_score and human_score > 5: species = "human"
    elif mouse_score > human_score and mouse_score > 5: species = "mouse"
    else:                                               species = "unknown"

    print(f"  → Detected: {species}")
    if species == "unknown":
        print("  Set species manually: compute_qc_metrics_human/mouse(adata)")
    return species


def convert_ensembl_to_symbols(adata: ad.AnnData, species: str = "human") -> ad.AnnData:
    """
    Convert ENSEMBL gene IDs to gene symbols via mygene.

    Required when var_names contain ENSG/ENSMUSG IDs instead of
    gene symbols. QC gene-set patterns (MT-, RPS, etc.) need symbols.

    Parameters
    ----------
    adata   : AnnData with ENSEMBL IDs as var_names.
    species : 'human' or 'mouse'

    Returns
    -------
    AnnData  var_names replaced by symbols where available.
             Unmapped IDs kept as-is.

    Example
    -------
    >>> adata.var_names[:2].tolist()
    ['ENSG00000000003', 'ENSG00000000005']
    >>> adata = convert_ensembl_to_symbols(adata, species="human")
    >>> adata.var_names[:2].tolist()
    ['TSPAN6', 'TNMD']
    """
    try:
        import mygene
    except ImportError:
        raise ImportError("pip install mygene")

    mg     = mygene.MyGeneInfo()
    ids    = adata.var_names.tolist()
    taxon  = 9606 if species == "human" else 10090

    print(f"[ensembl] Querying {len(ids)} IDs via mygene...")
    result  = mg.querymany(ids, scopes="ensembl.gene", fields="symbol",
                           species=taxon, as_dataframe=True, verbose=False)
    mapping = result["symbol"].dropna().to_dict() if "symbol" in result.columns else {}

    n_mapped   = sum(1 for i in ids if i in mapping)
    new_names  = [mapping.get(i, i) for i in ids]
    print(f"[ensembl] Mapped: {n_mapped}/{len(ids)} genes")

    adata.var_names = new_names
    adata.var_names_make_unique()
    return adata


def handle_duplicate_genes(adata: ad.AnnData, strategy: str = "sum") -> ad.AnnData:
    """
    Resolve duplicate gene names in var_names.

    Strategies
    ----------
    'sum'    : sum counts of duplicated genes (default — no signal lost)
    'keep'   : keep first occurrence only (fast, slight loss)
    'unique' : append _1, _2 suffix (no loss, no merging)

    Parameters
    ----------
    adata    : AnnData
    strategy : 'sum' | 'keep' | 'unique'

    Returns
    -------
    AnnData

    Example
    -------
    >>> n_dups = adata.var_names.duplicated().sum()
    >>> print(f"Duplicates: {n_dups}")
    >>> adata = handle_duplicate_genes(adata, strategy="sum")
    # [dedup] 47 duplicates | strategy: sum → 33491 genes remaining
    """
    n_dups = adata.var_names.duplicated().sum()
    if n_dups == 0:
        print("[dedup] No duplicate gene names found")
        return adata

    print(f"[dedup] {n_dups} duplicates | strategy: {strategy}")

    if strategy == "unique":
        adata.var_names_make_unique()
        print(f"[dedup] Suffixes appended → {adata.n_vars} unique names")
        return adata

    if strategy == "keep":
        mask  = ~adata.var_names.duplicated(keep="first")
        adata = adata[:, mask].copy()
        print(f"[dedup] First kept → {adata.n_vars} genes remaining")
        return adata

    if strategy == "sum":
        X  = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
        df = df.T.groupby(level=0).sum().T
        return ad.AnnData(X=sparse.csr_matrix(df.values.astype(np.float32)),
                          obs=adata.obs.copy(),
                          var=pd.DataFrame(index=df.columns))

    raise ValueError(f"strategy must be 'sum', 'keep', or 'unique'. Got '{strategy}'")


def inject_metadata(adata: ad.AnnData,
                    metadata_path: str,
                    on: str = "sample") -> ad.AnnData:
    """
    Join experimental/clinical metadata CSV into adata.obs by sample ID.

    Parameters
    ----------
    adata         : AnnData with adata.obs[on] containing sample IDs.
    metadata_path : CSV or TSV with a column matching `on`.
    on            : obs column to join on (default 'sample').

    Returns
    -------
    AnnData  With metadata columns added to adata.obs.

    Example
    -------
    CSV format (sample_metadata.csv):
        sample,condition,patient,batch
        GSM5139852,tumor,RB01,1
        GSM5139853,normal,RB01,1

    >>> adata = inject_metadata(adata, "../data/sample_metadata.csv")
    >>> adata.obs[['sample', 'condition', 'patient']].head()
    """
    sep  = "\t" if metadata_path.endswith(".tsv") else ","
    meta = pd.read_csv(metadata_path, sep=sep)

    if on not in meta.columns:
        raise ValueError(f"Column '{on}' not found. Available: {meta.columns.tolist()}")

    meta     = meta.set_index(on)
    new_cols = [c for c in meta.columns if c not in adata.obs.columns]
    print(f"[metadata] Injecting {len(new_cols)} columns: {new_cols}")

    adata.obs = adata.obs.join(meta[new_cols], on=on, how="left")
    missing   = adata.obs[new_cols[0]].isna().sum() if new_cols else 0
    print(f"  {'All cells matched ✓' if not missing else f'WARNING: {missing} cells unmatched'}")
    return adata


def load_if_not_exists(h5ad_path: str, compute_fn, *args, **kwargs) -> ad.AnnData:
    """
    Load existing h5ad if present; otherwise run compute_fn to create it.

    Avoids recomputing expensive steps when re-running notebooks.

    Parameters
    ----------
    h5ad_path  : str       Expected .h5ad output path.
    compute_fn : callable  Function that returns AnnData (e.g. run_pipeline).
    *args, **kwargs        Forwarded to compute_fn if recomputing.

    Returns
    -------
    AnnData

    Example
    -------
    >>> adata = load_if_not_exists(
    ...     "../results/qc_filtered.h5ad",
    ...     run_pipeline,
    ...     tar_path="../data/GSE168434_RAW.tar",
    ...     species="human",
    ...     out_dir="../results/",
    ... )
    # File exists:   [resume] Loaded: 70287 × 28717
    # File missing:  [resume] Running pipeline...
    """
    if Path(h5ad_path).exists():
        print(f"[resume] Found: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        print(f"[resume] Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        return adata
    print(f"[resume] {h5ad_path} not found — running pipeline...")
    return compute_fn(*args, **kwargs)


def diagnose_sample(folder: str, sample_id: str = "sample") -> None:
    """
    Debug a sample that failed to load.

    Prints folder existence, all files with sizes, detected format,
    and a preview of the first text file found.

    Parameters
    ----------
    folder    : str   Path to the failed sample folder.
    sample_id : str   Label for print output.

    Example
    -------
    >>> diagnose_sample("../data/raw/GSE168434_RAW/GSM5139860_failed")
    # [diagnose] Files (1):
    #   GSM5139860_counts.mtx.tsv   1024.3 KB
    # Detected format: txt
    # Preview — first 3 rows shown
    """
    folder = Path(folder)
    print(f"[diagnose] {sample_id}  →  {folder}")
    print(f"  Exists: {folder.exists()}")

    if not folder.exists():
        print("  ERROR: folder does not exist")
        return

    files = sorted(folder.iterdir())
    print(f"  Files ({len(files)}):")
    for f in files:
        print(f"    {f.name:<55} {f.stat().st_size/1024:>8.1f} KB")

    fmt = detect_format(str(folder))
    print(f"\n  Detected format: {fmt}")
    if fmt == "unknown":
        print("  HINT: no recognised data file. Check extensions above.")
        return

    for f in files:
        if any(f.name.endswith(e) for e in (".tsv", ".csv", ".txt", ".tsv.gz", ".csv.gz")):
            try:
                df = pd.read_csv(str(f), sep="\t", nrows=3, header=None)
                print(f"\n  Preview — {f.name}:")
                print(df.to_string(index=False))
            except Exception as e:
                print(f"\n  Preview failed for {f.name}: {e}")
            break
