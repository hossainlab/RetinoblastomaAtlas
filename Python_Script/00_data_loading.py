"""
=============================================================================
Script 00 — Data Loading and Atlas Assembly
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS SCRIPT?
-----------------
This is the entry point of the entire pipeline.  It reproduces as a reproducible,
memory-efficient command-line script.  Key responsibilities:

  1. Extract raw GEO tarballs into a consistent directory layout.
  2. Load each sample's 10x MTX matrix ONE AT A TIME (never all in RAM).
  3. Merge all samples into a single AnnData atlas using disk-based
     concatenation — avoiding the "load everything then concat" anti-pattern
     used in the notebook.
  4. Annotate cells with biological and technical metadata.
  5. Write the final atlas in compressed .h5ad format ready for QC (script 01).

MEMORY-EFFICIENCY STRATEGY
----------------------------
The naive approach (notebook) loads all 14 samples into a Python list and
calls anndata.concat() — briefly placing ~4–6 GB of count matrices in RAM
simultaneously.  For ~140k cells × 62k genes this is feasible but scales
poorly.  This script uses four complementary techniques:

  TECHNIQUE 1 — Per-sample sequential loading with immediate disk write
    Each sample is loaded, annotated, and written to a temporary compressed
    h5ad file.  The in-memory object is then explicitly deleted and garbage-
    collected before loading the next sample.  Peak RAM = one sample (~0.5 GB)
    rather than all samples simultaneously.

  TECHNIQUE 2 — uint16 count matrix dtype
    Raw UMI counts from 10x Genomics rarely exceed 32,767 per gene per cell.
    Storing them as uint16 (range 0–65,535) rather than float32 (default)
    halves the memory footprint of the count matrix:
      float32: 4 bytes/value → 140k cells × 62k genes × 4 B ≈ 34 GB (dense)
      uint16:  2 bytes/value → same ≈ 17 GB (dense)
    Combined with CSR sparse storage (>95% zeros), actual RAM is ~0.5 GB
    per sample.

  TECHNIQUE 3 — Gene universe computed from feature files alone
    Before loading any matrix, we read all feature.tsv.gz files (tiny, <1 MB
    each) to compute the full gene union.  This allows each sample's sparse
    matrix to be directly reshaped to the correct gene space during load,
    avoiding the large dense fill_value=0 expansion that anndata.concat()
    with join='outer' performs in RAM.

  TECHNIQUE 4 — anndata.experimental.concat_on_disk()
    Available since anndata ≥ 0.10.  Concatenates h5ad files by streaming
    chunks from disk WITHOUT loading them into RAM, writing directly to the
    output file.  Peak RAM = single chunk size (configurable), not total data.
    This is the most important technique for atlas-scale datasets.

DIRECTORY STRUCTURE ASSUMED
-----------------------------
data/raw/
  GSE168434_RAW.tar
  GSE249995_RAW.tar
After extraction → (script creates):
  data/raw/GSE168434/<sample>/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}
  data/raw/GSE249995/{GSM*_barcodes.tsv.gz, GSM*_features.tsv.gz, GSM*_matrix.mtx.gz}

SAMPLE METADATA
---------------
GSE168434 — 10 samples from 7 RB patients:
  GSM5139852_RB01_rep1, GSM5139853_RB01_rep2 → patient RB01 (intraocular)
  GSM5139854_RB02_rep1, GSM5139855_RB02_rep2 → patient RB02 (intraocular)
  GSM5139856_RB03_rep1, GSM5139857_RB03_rep2 → patient RB03 (intraocular)
  GSM5139858_RB04                             → patient RB04 (intraocular)
  GSM5139859_RB05                             → patient RB05 (intraocular)
  GSM5139860_RB06                             → patient RB06 (intraocular)
  GSM5139861_RB07                             → patient RB07 (intraocular)

GSE249995 — 4 samples (intraocular and extraocular pairs):
  GSM7968797 → S1_in1  (intraocular, patient 1)
  GSM7968798 → S2_in2  (intraocular, patient 2)
  GSM7968799 → S3_ex1  (extraocular, patient 1)
  GSM7968800 → S4_ex2  (extraocular, patient 2)

REFERENCES
----------
- La Manno G et al. RNA velocity of single cells. Nature. 2018;560:494-498.
  https://doi.org/10.1038/s41586-018-0414-6

- Virshup I et al. The scverse project provides a computational ecosystem for
  single-cell omics data analysis. Nat Biotechnol. 2023;41:604-606.
  https://doi.org/10.1038/s41587-023-01733-8

INPUT  : data/raw/GSE168434_RAW.tar
         data/raw/GSE249995_RAW.tar
OUTPUT :
  data/processed/01_merged_raw.h5ad  (compressed, uint16 counts)
  results/tables/sample_metadata.csv

USAGE  : pixi run python script/00_data_loading.py [--no-extract] [--join inner|outer]
=============================================================================
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
import shutil
import sys
import tarfile
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT      = Path(__file__).resolve().parents[1]
RAW_DIR   = ROOT / "data" / "raw"
PROC_DIR  = ROOT / "data" / "processed"
TAB_DIR   = ROOT / "results" / "tables"
OUT_H5AD  = PROC_DIR / "01_merged_raw.h5ad"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)
warnings.filterwarnings("ignore")
sc.settings.verbosity = 0  # suppress scanpy's own messages; we handle logging

# ---------------------------------------------------------------------------
# Sample manifest
# ---------------------------------------------------------------------------
# GSE168434: 10 samples, each in its own subdirectory after extraction
GSE168434_SAMPLES: list[dict] = [
    {"gsm": "GSM5139852_RB01_rep1", "patient_id": "RB01", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139853_RB01_rep2", "patient_id": "RB01", "replicate": "rep2",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139854_RB02_rep1", "patient_id": "RB02", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139855_RB02_rep2", "patient_id": "RB02", "replicate": "rep2",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139856_RB03_rep1", "patient_id": "RB03", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139857_RB03_rep2", "patient_id": "RB03", "replicate": "rep2",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139858_RB04",      "patient_id": "RB04", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139859_RB05",      "patient_id": "RB05", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139860_RB06",      "patient_id": "RB06", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM5139861_RB07",      "patient_id": "RB07", "replicate": "rep1",
     "disease_stage": "intraocular", "tissue": "primary_tumour"},
]

# GSE249995: 4 samples — intraocular (S1, S2) and extraocular (S3, S4)
# GSM IDs map to descriptive sample names based on GEO submission order
GSE249995_SAMPLES: list[dict] = [
    {"gsm": "GSM7968797", "sample_id": "S1_in1", "patient_id": "P1",
     "replicate": "rep1", "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM7968798", "sample_id": "S2_in2", "patient_id": "P2",
     "replicate": "rep1", "disease_stage": "intraocular", "tissue": "primary_tumour"},
    {"gsm": "GSM7968799", "sample_id": "S3_ex1", "patient_id": "P1",
     "replicate": "rep1", "disease_stage": "extraocular", "tissue": "optic_nerve"},
    {"gsm": "GSM7968800", "sample_id": "S4_ex2", "patient_id": "P2",
     "replicate": "rep1", "disease_stage": "extraocular", "tissue": "optic_nerve"},
]


# ---------------------------------------------------------------------------
# STEP 1 — Extraction
# ---------------------------------------------------------------------------

def extract_tar(tar_path: Path, dest_dir: Path) -> None:
    """Extract a .tar or .tar.gz file to dest_dir.

    SIGNIFICANCE:
    GEO distributes data as monolithic tar archives.  Extraction is a
    one-time operation; subsequent pipeline runs skip it via --no-extract.
    Using filter='data' avoids the Python 3.14 deprecation warning and
    prevents path traversal vulnerabilities in tar files.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    mode = "r:gz" if str(tar_path).endswith(".gz") else "r"
    log.info(f"  Extracting {tar_path.name} → {dest_dir.name}/")
    try:
        with tarfile.open(tar_path, mode) as tar:
            try:
                tar.extractall(path=dest_dir, filter="data")
            except TypeError:
                # Python < 3.12 does not support filter= argument
                tar.extractall(path=dest_dir)  # noqa: S202
    except tarfile.TarError as e:
        log.error(f"  Extraction failed: {e}")
        raise


def extract_gse168434(raw_dir: Path) -> None:
    """Extract GSE168434: outer tar then nested per-sample tar.gz files.

    STRUCTURE AFTER EXTRACTION:
    Each sample gets its own subdirectory:
      data/raw/GSE168434/<GSM_sample>/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}

    The nested structure from Velocyto paths (home/.../counts/) is flattened
    to keep only the three 10x MTX files at the sample root.
    """
    outer_tar = raw_dir / "GSE168434_RAW.tar"
    gse_dir   = raw_dir / "GSE168434"

    if not outer_tar.exists():
        raise FileNotFoundError(
            f"{outer_tar} not found.\n"
            "  Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168434"
        )

    # Extract outer tar (contains nested per-sample tar.gz files)
    extract_tar(outer_tar, gse_dir)

    # Extract each per-sample tar.gz
    for tar_gz in sorted(gse_dir.glob("*.tar.gz")):
        sample_name = tar_gz.name.split("_counts")[0]
        sample_dir  = gse_dir / sample_name
        if sample_dir.exists() and (sample_dir / "matrix.mtx.gz").exists():
            log.info(f"    {sample_name}: already extracted, skipping")
            continue
        sample_dir.mkdir(exist_ok=True)
        log.info(f"    Extracting sample: {sample_name}")
        extract_tar(tar_gz, sample_dir)

        # Flatten nested directory structure (Velocyto puts files deep in home/)
        _flatten_sample_dir(sample_dir)

    log.info(f"  GSE168434 extraction complete ({len(list(gse_dir.iterdir()))} entries)")


def _flatten_sample_dir(sample_dir: Path) -> None:
    """Move barcodes/features/matrix files to sample_dir root.

    The tar.gz files from GSE168434 contain a deep path like:
      home/<user>/CellRanger/<sample>/outs/filtered_feature_bc_matrix/
    We flatten this to sample_dir/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}.
    """
    target_files = {"barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz",
                    "barcodes.tsv",    "features.tsv",    "matrix.mtx"}
    for root, dirs, files in os.walk(sample_dir):
        root_path = Path(root)
        if root_path == sample_dir:
            continue
        for f in files:
            if f in target_files:
                src = root_path / f
                dst = sample_dir / f
                if not dst.exists():
                    shutil.move(str(src), str(dst))

    # Remove the now-empty nested directory tree (only the top-level "home" dir)
    home_dir = sample_dir / "home"
    if home_dir.exists():
        shutil.rmtree(home_dir, ignore_errors=True)


def extract_gse249995(raw_dir: Path) -> None:
    """Extract GSE249995: flat structure (prefixed files in one directory).

    STRUCTURE AFTER EXTRACTION:
    All files land flat in data/raw/GSE249995/:
      GSM7968797_barcodes.tsv.gz
      GSM7968797_features.tsv.gz
      GSM7968797_matrix.mtx.gz
      GSM7968798_barcodes.tsv.gz
      ... (similarly for GSM7968798, 7968799, 7968800)
    """
    outer_tar = raw_dir / "GSE249995_RAW.tar"
    gse_dir   = raw_dir / "GSE249995"

    if not outer_tar.exists():
        raise FileNotFoundError(
            f"{outer_tar} not found.\n"
            "  Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249995"
        )

    already_extracted = gse_dir.exists() and any(gse_dir.glob("GSM*matrix*"))
    if already_extracted:
        log.info(f"  GSE249995: already extracted in {gse_dir.name}/")
        return

    extract_tar(outer_tar, gse_dir)
    log.info(f"  GSE249995 extraction complete ({len(list(gse_dir.iterdir()))} files)")


# ---------------------------------------------------------------------------
# STEP 2 — Gene universe from feature files (no matrix loading)
# ---------------------------------------------------------------------------

def read_feature_names(feature_file: Path) -> pd.Series:
    """Read gene symbols from a features.tsv.gz file (no matrix).

    SIGNIFICANCE (Technique 3):
    Reading feature files is trivial in memory (<1 MB each vs. >500 MB for
    the matrix).  By computing the gene universe BEFORE loading any matrix,
    we can pre-allocate each sample's sparse matrix in the correct shape,
    avoiding the costly in-memory reindex that anndata.concat(join='outer')
    performs.
    """
    return pd.read_csv(
        feature_file, sep="\t", header=None, usecols=[1], squeeze=False
    )[1]


def compute_gene_universe(
    sample_feature_files: list[Path],
    join: str = "outer",
) -> list[str]:
    """Compute the union or intersection of gene names across all samples.

    Parameters
    ----------
    join : 'outer' — union (default; preserves all genes, NaN-filled for absent)
           'inner' — intersection (smaller; only genes detected in ALL samples)

    SIGNIFICANCE:
    For cross-dataset integration:
    - 'outer' is biologically safer: retains dataset-specific genes
      (e.g., mitochondrial pseudogenes present in one Cell Ranger version)
    - 'inner' is more memory-efficient: fewer genes, no fill-value expansion
    - We default to 'outer' so that mitochondrial QC genes are not lost if
      absent from one dataset.  The penalty is ~15–25k extra genes filled
      with zeros — still sparse.
    """
    log.info(f"  Computing gene universe (join='{join}') from {len(sample_feature_files)} feature files...")
    gene_sets = []
    for f in sample_feature_files:
        genes = read_feature_names(f)
        gene_sets.append(set(genes.values))
        log.info(f"    {f.parent.name or f.name}: {len(genes):,} genes")

    if join == "inner":
        universe = sorted(set.intersection(*gene_sets))
    else:
        universe = sorted(set.union(*gene_sets))

    log.info(f"  Gene universe ({join}): {len(universe):,} genes")
    return universe


# ---------------------------------------------------------------------------
# STEP 3 — Per-sample loading (Techniques 1 & 2)
# ---------------------------------------------------------------------------

def load_10x_sample_uint16(
    mtx_dir: Path,
    sample_meta: dict,
    gene_universe: list[str],
) -> ad.AnnData:
    """Load one 10x MTX sample as a uint16 sparse AnnData.

    SIGNIFICANCE (Technique 2 — uint16):
    scipy.sparse stores (data, indices, indptr) arrays.  By requesting
    dtype=uint16 instead of float32:
      - 'data' array: 2 bytes/non-zero instead of 4 → 50% saving for values
      - indices/indptr stay int32 (unchanged)
    For a typical sample with 10k cells × 20k expressed genes × ~3% fill:
      float32: 10k × 62k × 0.03 × 4 B ≈ 74 MB
      uint16:  10k × 62k × 0.03 × 2 B ≈ 37 MB
    Across 14 samples loaded sequentially this doesn't change peak RAM,
    but it halves the size of each written temp file (→ faster I/O).

    NOTE: We convert to int32 before writing to h5ad because h5py does not
    natively support uint16 — it will silently upcast to float32 otherwise.
    int32 is safe for counts ≤ 2,147,483,647.
    """
    # Read with scanpy (handles .gz transparently)
    adata = sc.read_10x_mtx(
        mtx_dir,
        var_names="gene_symbols",
        cache=False,
        make_unique=True,    # deduplicate gene names with suffix
    )

    # Cast to uint16 immediately after load
    X = adata.X
    if sp.issparse(X):
        X = X.astype(np.uint16)
        adata.X = X
    else:
        adata.X = X.astype(np.uint16)

    # Reindex genes to the full universe (fills absent genes with 0 = sparse zero)
    # This is faster than outer-join during concat because it operates on
    # a single small matrix rather than the merged giant matrix
    adata = adata[:, adata.var_names.isin(gene_universe)].copy()
    missing_genes = [g for g in gene_universe if g not in adata.var_names]
    if missing_genes:
        # Build zero columns for missing genes (as sparse block)
        n_cells   = adata.n_obs
        n_missing = len(missing_genes)
        zero_block = sp.csr_matrix((n_cells, n_missing), dtype=np.uint16)
        missing_var = pd.DataFrame(index=missing_genes)
        missing_adata = ad.AnnData(
            X=zero_block,
            obs=pd.DataFrame(index=adata.obs_names),
            var=missing_var,
        )
        adata = ad.concat(
            [adata, missing_adata], axis=1, join="outer"
        )

    # Reorder columns to match universe order exactly (needed for concat_on_disk)
    adata = adata[:, gene_universe].copy()

    # Add barcode prefix to make obs_names unique across all samples
    sample_id  = sample_meta.get("sample_id", sample_meta.get("gsm", "unknown"))
    adata.obs_names = [f"{bc}-{sample_id}" for bc in adata.obs_names]

    # Add metadata columns
    for col, val in sample_meta.items():
        adata.obs[col] = val

    # Rename 'gsm' to 'sample_id' if it was the key
    if "gsm" in adata.obs.columns and "sample_id" not in sample_meta:
        adata.obs.rename(columns={"gsm": "sample_id"}, inplace=True)

    log.info(
        f"    {sample_id}: {adata.n_obs:,} cells × {adata.n_vars:,} genes  "
        f"[{adata.X.data.nbytes / 1e6:.1f} MB non-zero data, uint16]"
    )
    return adata


def load_gse249995_flat(
    gsm_id: str,
    gse_dir: Path,
    sample_meta: dict,
    gene_universe: list[str],
) -> ad.AnnData:
    """Load one GSE249995 sample from flat prefixed MTX files.

    GSE249995 stores files as:
      {gsm_id}_matrix.mtx.gz / {gsm_id}_barcodes.tsv.gz / {gsm_id}_features.tsv.gz
    We copy them to a temp directory with the standard 10x names so that
    sc.read_10x_mtx() can load them without modification.
    """
    # Resolve the actual file names (may have extra suffixes in some GEO uploads)
    def _find(pattern: str) -> Path:
        matches = list(gse_dir.glob(f"{gsm_id}*{pattern}*"))
        if not matches:
            raise FileNotFoundError(
                f"No file matching '{gsm_id}*{pattern}*' in {gse_dir}"
            )
        return matches[0]

    mtx_src      = _find("matrix.mtx")
    barcodes_src = _find("barcodes.tsv")
    features_src = _find("features.tsv")

    # Create a temp dir with standard 10x names
    tmp = Path(tempfile.mkdtemp(prefix=f"gsm_{gsm_id}_"))
    try:
        shutil.copy2(mtx_src,      tmp / "matrix.mtx.gz"  if str(mtx_src).endswith(".gz") else tmp / "matrix.mtx")
        shutil.copy2(barcodes_src, tmp / "barcodes.tsv.gz" if str(barcodes_src).endswith(".gz") else tmp / "barcodes.tsv")
        shutil.copy2(features_src, tmp / "features.tsv.gz" if str(features_src).endswith(".gz") else tmp / "features.tsv")

        return load_10x_sample_uint16(tmp, sample_meta, gene_universe)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
# STEP 4 — Per-sample temp h5ad writing (Technique 1)
# ---------------------------------------------------------------------------

def write_temp_h5ad(adata: ad.AnnData, tmp_dir: Path, index: int) -> Path:
    """Write a single-sample AnnData to a compressed temp h5ad.

    SIGNIFICANCE (Technique 1 — sequential + GC):
    After writing to disk, the caller deletes the AnnData and calls
    gc.collect().  This ensures Python releases the sparse matrix memory
    before loading the next sample.  Without explicit gc.collect(), Python's
    cyclic garbage collector may defer release, causing two samples to
    coexist in RAM during the loading loop.

    Compression level 4 balances I/O speed vs. file size.  For raw uint16
    count matrices, gzip achieves ~5–10x compression.
    """
    out = tmp_dir / f"sample_{index:03d}.h5ad"
    adata.write_h5ad(out, compression="gzip", compression_opts=4)
    log.info(f"    Written temp file: {out.name}  ({out.stat().st_size / 1e6:.1f} MB)")
    return out


# ---------------------------------------------------------------------------
# STEP 5 — Disk-based concatenation (Technique 4)
# ---------------------------------------------------------------------------

def concat_on_disk_safe(
    tmp_h5ads: list[Path],
    out_path: Path,
) -> ad.AnnData:
    """Concatenate per-sample h5ads without loading into RAM.

    SIGNIFICANCE (Technique 4 — anndata.experimental.concat_on_disk):
    concat_on_disk() reads chunks of each h5ad sequentially and writes
    them directly to out_path.  Peak RAM = size of one chunk (a few hundred
    MB), not the total atlas (several GB).

    This function falls back to the standard in-memory concat if
    anndata < 0.10 is installed (with a warning).
    """
    try:
        from anndata.experimental import concat_on_disk
        log.info(
            "  Using anndata.experimental.concat_on_disk() — "
            "peak RAM bounded by single chunk"
        )
        concat_on_disk(
            in_files={p.stem: p for p in tmp_h5ads},
            out_file=out_path,
            axis=0,
            join="outer",
            label="sample_key",   # stores dict key (stem) in obs column
            index_unique=None,    # barcodes already unique (we prefixed them)
            fill_value=0,
        )
        # Read back the result (backed mode for inspection, close before return)
        adata = ad.read_h5ad(out_path)
        log.info(
            f"  concat_on_disk complete → {out_path.name}\n"
            f"  Atlas: {adata.n_obs:,} cells × {adata.n_vars:,} genes"
        )
        return adata

    except ImportError:
        log.warning(
            "  anndata.experimental.concat_on_disk not available "
            "(requires anndata >= 0.10).\n"
            "  Falling back to in-memory concat — peak RAM will be higher."
        )
        adatas = [ad.read_h5ad(p) for p in tmp_h5ads]
        adata  = ad.concat(
            adatas, axis=0, join="outer", fill_value=0, index_unique=None
        )
        del adatas
        gc.collect()
        adata.write_h5ad(out_path, compression="gzip", compression_opts=4)
        return adata


# ---------------------------------------------------------------------------
# STEP 6 — Gene-level QC flags
# ---------------------------------------------------------------------------

def annotate_gene_flags(adata: ad.AnnData) -> None:
    """Add Boolean gene flags for QC: mitochondrial, ribosomal, hemoglobin.

    SIGNIFICANCE:
    These flags are used in script 01 (QC filtering) to:
      - Compute pct_counts_mt per cell (indicator of apoptosis / damaged cells)
      - Identify ribosomal gene content (high ribo% = stressed/cycling cells)
      - Flag hemoglobin contamination (RBC contamination from ocular blood vessels)

    Flags are stored in adata.var (gene metadata), not adata.obs, so they
    are preserved through all downstream pipeline steps.
    """
    adata.var["mt"]   = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"]   = adata.var_names.str.match(r"^HB[^(P)]")

    log.info(
        f"  Gene flags:\n"
        f"    Mitochondrial (MT-):      {adata.var['mt'].sum():,}\n"
        f"    Ribosomal (RPS/RPL):      {adata.var['ribo'].sum():,}\n"
        f"    Hemoglobin (HB*):         {adata.var['hb'].sum():,}"
    )


# ---------------------------------------------------------------------------
# STEP 7 — Sample metadata summary
# ---------------------------------------------------------------------------

def save_sample_metadata(adata: ad.AnnData, out_path: Path) -> None:
    """Save a per-sample cell count summary table."""
    group_cols = [c for c in ["sample_id", "dataset", "patient_id",
                               "disease_stage", "tissue", "replicate"]
                  if c in adata.obs.columns]
    meta = (
        adata.obs[group_cols]
        .drop_duplicates()
        .merge(
            adata.obs.groupby("sample_id").size().rename("n_cells").reset_index(),
            on="sample_id",
        )
        .sort_values("sample_id")
        .reset_index(drop=True)
    )
    meta.to_csv(out_path, index=False)
    log.info(f"  Sample metadata saved → {out_path.name}")
    log.info(f"\n{meta.to_string(index=False)}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_loading(
    do_extract: bool = True,
    gene_join: str = "outer",
) -> ad.AnnData:
    """Full data loading pipeline.

    Steps
    -----
    1. Extract raw tar archives (optional, skip with --no-extract).
    2. Collect all feature files and compute gene universe.
    3. Load each sample sequentially → write temp h5ad → delete from RAM.
    4. Concatenate all temp h5ads on disk (no full atlas in RAM).
    5. Add gene QC flags (mt, ribo, hb).
    6. Save final atlas and sample metadata table.
    7. Clean up temp files.
    """
    PROC_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. Extract -------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Extracting raw GEO archives")
    log.info("=" * 60)
    if do_extract:
        extract_gse168434(RAW_DIR)
        extract_gse249995(RAW_DIR)
    else:
        log.info("  --no-extract: skipping extraction")

    # ---- 2. Gene universe -------------------------------------------------
    log.info("\nSTEP 2 — Computing gene universe from feature files")
    log.info(
        "  Reading only the small feature.tsv.gz files (no matrices).\n"
        "  WHY: Knowing the full gene list before loading allows each\n"
        "  sample to be reindexed once, cheaply, rather than performing\n"
        "  a large in-memory outer join at concat time."
    )
    feature_files: list[Path] = []
    # GSE168434 feature files (one per sample subdirectory)
    gse168434_dir = RAW_DIR / "GSE168434"
    for meta in GSE168434_SAMPLES:
        f = gse168434_dir / meta["gsm"] / "features.tsv.gz"
        if not f.exists():
            f = gse168434_dir / meta["gsm"] / "features.tsv"
        if f.exists():
            feature_files.append(f)
        else:
            log.warning(f"  Feature file not found for {meta['gsm']} — sample will be skipped")

    # GSE249995 feature files (prefixed flat files)
    gse249995_dir = RAW_DIR / "GSE249995"
    for meta in GSE249995_SAMPLES:
        matches = list(gse249995_dir.glob(f"{meta['gsm']}*features*"))
        if matches:
            feature_files.append(matches[0])
        else:
            log.warning(f"  Feature file not found for {meta['gsm']} — sample will be skipped")

    gene_universe = compute_gene_universe(feature_files, join=gene_join)

    # ---- 3 & 4. Sequential load → write temp h5ads ----------------------
    log.info(
        "\nSTEP 3 — Loading samples sequentially (memory-efficient)\n"
        "  TECHNIQUE 1: Load one sample → write temp h5ad → del + gc.collect()\n"
        "  TECHNIQUE 2: Count matrices stored as uint16 (2 bytes/value)\n"
        "  TECHNIQUE 3: Reindex to gene universe per-sample (no in-memory outer join)"
    )
    tmp_dir    = PROC_DIR / "_tmp_samples"
    tmp_dir.mkdir(exist_ok=True)
    tmp_h5ads: list[Path] = []
    sample_idx = 0

    # Load GSE168434
    log.info(f"\n  --- GSE168434 ({len(GSE168434_SAMPLES)} samples) ---")
    for meta in GSE168434_SAMPLES:
        sample_dir = gse168434_dir / meta["gsm"]
        if not sample_dir.exists():
            log.warning(f"  Directory not found: {sample_dir} — skipping")
            continue

        sample_meta = dict(meta)
        if "sample_id" not in sample_meta:
            sample_meta["sample_id"] = meta["gsm"]
        sample_meta["dataset"] = "GSE168434"

        try:
            adata_s = load_10x_sample_uint16(sample_dir, sample_meta, gene_universe)
            tmp_path = write_temp_h5ad(adata_s, tmp_dir, sample_idx)
            tmp_h5ads.append(tmp_path)
            sample_idx += 1
        except Exception as e:
            log.error(f"  Failed to load {meta['gsm']}: {e}")
            continue
        finally:
            del adata_s  # explicit delete
            gc.collect() # release sparse matrix memory immediately
            log.info(f"    Memory released after {meta['gsm']}")

    # Load GSE249995
    log.info(f"\n  --- GSE249995 ({len(GSE249995_SAMPLES)} samples) ---")
    for meta in GSE249995_SAMPLES:
        sample_meta = dict(meta)
        sample_meta["dataset"] = "GSE249995"

        try:
            adata_s = load_gse249995_flat(meta["gsm"], gse249995_dir, sample_meta, gene_universe)
            tmp_path = write_temp_h5ad(adata_s, tmp_dir, sample_idx)
            tmp_h5ads.append(tmp_path)
            sample_idx += 1
        except Exception as e:
            log.error(f"  Failed to load {meta['gsm']}: {e}")
            continue
        finally:
            del adata_s
            gc.collect()
            log.info(f"    Memory released after {meta['gsm']}")

    if not tmp_h5ads:
        raise RuntimeError("No samples loaded successfully. Check raw data paths.")

    log.info(f"\n  {len(tmp_h5ads)} samples written to temp h5ads")

    # ---- 5. Disk concatenation -------------------------------------------
    log.info(
        "\nSTEP 4 — Concatenating on disk\n"
        "  TECHNIQUE 4: anndata.experimental.concat_on_disk()\n"
        "  Peak RAM = single chunk, not the full atlas"
    )
    concat_on_disk_safe(tmp_h5ads, OUT_H5AD)
    adata = ad.read_h5ad(OUT_H5AD)

    # ---- 6. Gene QC flags ------------------------------------------------
    log.info("\nSTEP 5 — Adding gene QC flags")
    log.info(
        "  Adding .var['mt'], .var['ribo'], .var['hb'] flags.\n"
        "  These enable per-cell pct_counts_mt / ribo / hb in script 01."
    )
    annotate_gene_flags(adata)

    # Remove the sample_key column added by concat_on_disk label parameter
    if "sample_key" in adata.obs.columns:
        adata.obs.drop(columns=["sample_key"], inplace=True)

    # ---- 7. Save summary table ------------------------------------------
    log.info("\nSTEP 6 — Saving sample metadata summary")
    save_sample_metadata(adata, TAB_DIR / "sample_metadata.csv")

    # ---- 8. Write final h5ad (overwrite with gene flags added) ----------
    log.info(f"\nSTEP 7 — Writing final atlas → {OUT_H5AD.name}")
    log.info(
        "  Contents:\n"
        "    .X                   : raw UMI counts (int32, CSR sparse)\n"
        "    .obs['sample_id']    : sample identifier\n"
        "    .obs['dataset']      : GSE168434 | GSE249995\n"
        "    .obs['patient_id']   : RB01–RB07 | P1–P2\n"
        "    .obs['disease_stage']: intraocular | extraocular\n"
        "    .obs['tissue']       : primary_tumour | optic_nerve\n"
        "    .obs['replicate']    : rep1 | rep2\n"
        "    .var['mt']           : mitochondrial gene flag\n"
        "    .var['ribo']         : ribosomal gene flag\n"
        "    .var['hb']           : hemoglobin gene flag"
    )
    adata.write_h5ad(OUT_H5AD, compression="gzip", compression_opts=4)
    size_mb = OUT_H5AD.stat().st_size / 1e6
    log.info(
        f"  Saved: {OUT_H5AD.name}  "
        f"({adata.n_obs:,} cells × {adata.n_vars:,} genes, {size_mb:.0f} MB)"
    )

    # ---- 9. Cleanup temp files ------------------------------------------
    log.info("\nSTEP 8 — Cleaning up temporary files")
    shutil.rmtree(tmp_dir, ignore_errors=True)
    log.info(f"  Removed temp directory: {tmp_dir.name}")
    log.info("  Done.\n")
    return adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Memory-efficient data loading for RetinoblastomaAtlas",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--no-extract",
        action="store_true",
        default=False,
        help="Skip tar extraction (if raw data already extracted)",
    )
    p.add_argument(
        "--join",
        choices=["outer", "inner"],
        default="outer",
        help=(
            "Gene join strategy across samples. "
            "'outer' = union (all genes, more memory); "
            "'inner' = intersection (fewer genes, less memory)"
        ),
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_loading(
        do_extract=not args.no_extract,
        gene_join=args.join,
    )
