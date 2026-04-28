"""
=============================================================================
Script 05 — Copy Number Variation (CNV) Inference
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
Retinoblastoma is initiated by biallelic inactivation of the RB1 tumour
suppressor (chromosome 13q14), but progression to extraocular disease
involves additional somatic copy number alterations (SCNAs). The most
consistently reported SCNAs in high-risk RB are:

  - Chromosome 6p amplification (E2F3, DEK, ID4) — associated with
    aggressive tumour behaviour and poor differentiation
  - Chromosome 1q gain (MDM4) — linked to p53 pathway evasion
  - Chromosome 16q loss — recurrent in metastatic disease
  - Chromosome 2p amplification (MYCN) — found in subtype 2 / MYCN-amplified
    RB with heritable predisposition and high proliferation index

Inferring CNVs from scRNA-seq is valuable because:
  1. It DISTINGUISHES tumour cells from non-malignant stromal/immune cells
     without needing matched WGS data — essential when only scRNA-seq is
     available (as in GSE168434 and GSE249995).
  2. It identifies INTRA-TUMORAL HETEROGENEITY: subclones defined by distinct
     SCNA profiles may differ in invasive potential, matching the
     intraocular → extraocular progression axis in this atlas.
  3. It validates cell-type annotation: cells annotated as Cone_precursor or
     Tumour_proliferating should show recurrent CNV patterns, while immune
     and glial cells should appear diploid.

METHODS — inferCNV vs. CopyKAT vs. scCNV-Score
-------------------------------------------------
We use inferCNPy (Python port of Broad Institute's inferCNV approach) via
the `cnvspy` package.  The method works by:
  1. Computing the mean expression of genes in consecutive chromosomal
     windows (101 genes) across cells.
  2. Normalising by the expression of reference (non-malignant) cells
     (Muller glia, microglia, endothelial cells annotated in script 04).
  3. Smoothing across windows and thresholding to call amplifications and
     deletions.

LIMITATIONS
-----------
scRNA-seq CNV inference has lower resolution than WGS/WES; it detects
large arm-level events (>5 Mb) reliably but misses focal alterations.
The signal also depends on the quality of the reference cell set.

REFERENCES
----------
- Tirosh I et al. Dissecting the multicellular ecosystem of metastatic
  melanoma by single-cell RNA-seq. Science. 2016;352(6282):189-196.
  https://doi.org/10.1126/science.aad0501

- Macosko EZ et al. Highly parallel genome-wide expression profiling of
  individual cells using nanoliter droplets. Cell. 2015;161(5):1202-1214.
  https://doi.org/10.1016/j.cell.2015.05.002

- Gao R et al. Delineating copy number and clonal substructure in human
  tumors from single-cell transcriptomes. Nat Biotechnol. 2021;39:599-608.
  https://doi.org/10.1038/s41587-020-00795-2

- Deng M et al. Identification of key pathways and genes in retinoblastoma
  using bioinformatics analysis. Gene. 2019;680:97-106.
  https://doi.org/10.1016/j.gene.2018.08.079

INPUT  : data/processed/05_annotated.h5ad
OUTPUT :
  data/processed/06_cnv.h5ad
  results/figures/cnv_heatmap_chromosome.pdf
  results/figures/cnv_umap_score.pdf
  results/tables/cnv_scores_per_cell.csv

USAGE  : pixi run python script/05_copy_number_variation.py
=============================================================================
"""

from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data" / "processed" / "05_annotated.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "06_cnv.h5ad"
FIG_DIR  = ROOT / "results" / "figures"
TAB_DIR  = ROOT / "results" / "tables"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Reference (non-malignant) cell types used to normalise CNV signal.
# These should be diploid and not carry RB1 mutation.
REFERENCE_CELL_TYPES = [
    "Muller_glia",
    "Microglia_TAM",
    "Endothelial",
    "Astrocyte",
]

# Chromosomal windows of interest for RB-specific SCNAs.
# Used for targeted summary statistics.
RB_SCNA_REGIONS = {
    "chr2p_MYCN":  ("2", 0,   50_000_000),
    "chr6p_E2F3":  ("6", 0,   65_000_000),
    "chr13q_RB1":  ("13", 45_000_000, 70_000_000),
    "chr16q_loss": ("16", 65_000_000, 120_000_000),
    "chr1q_gain":  ("1",  120_000_000, 250_000_000),
}

WINDOW_SIZE = 101  # genes per chromosomal window (standard inferCNV setting)


# ---------------------------------------------------------------------------
# CNV inference (sliding-window approach)
# ---------------------------------------------------------------------------

def get_gene_positions(adata: sc.AnnData) -> pd.DataFrame:
    """Map genes to chromosomal positions using Ensembl IDs in adata.var.

    SIGNIFICANCE:
    Chromosomal gene order is fundamental to CNV inference: amplified regions
    appear as contiguous blocks of high expression when genes are ordered by
    genomic position.  Without positional information, windowed averaging
    cannot detect SCNAs.

    The approach reads 'chromosome', 'start', 'end' columns from adata.var
    if they exist (populated by earlier pipeline steps), or falls back to
    a pre-compiled gene → position table for GRCh38.
    """
    required_cols = {"chromosome", "start"}
    if required_cols.issubset(set(adata.var.columns)):
        log.info("  Using chromosomal positions from adata.var")
        positions = adata.var[["chromosome", "start"]].copy()
    else:
        log.warning(
            "  adata.var lacks 'chromosome'/'start' columns. "
            "Attempting to load pre-compiled gene positions from "
            "data/processed/gene_positions_grch38.csv"
        )
        pos_file = ROOT / "data" / "processed" / "gene_positions_grch38.csv"
        if not pos_file.exists():
            raise FileNotFoundError(
                f"{pos_file} not found. Please run:\n"
                "  pixi run python script/00_fetch_gene_positions.py\n"
                "to generate gene position annotations from Ensembl."
            )
        positions = pd.read_csv(pos_file, index_col=0)

    # Keep only autosomal chromosomes (1–22)
    valid_chroms = [str(i) for i in range(1, 23)]
    positions = positions[positions["chromosome"].isin(valid_chroms)]
    positions["chromosome"] = positions["chromosome"].astype(int)
    positions = positions.sort_values(["chromosome", "start"])
    return positions


def compute_cnv_matrix(
    adata: sc.AnnData,
    positions: pd.DataFrame,
    reference_mask: np.ndarray,
    window_size: int = WINDOW_SIZE,
) -> np.ndarray:
    """Compute per-cell CNV scores using sliding chromosomal windows.

    Algorithm
    ---------
    1. Restrict to genes present in both adata.var and positions.
    2. Sort genes by (chromosome, start_position).
    3. For each window of `window_size` consecutive genes:
       a. Compute per-cell mean expression (from .layers['scvi_normalized']).
       b. Subtract the per-window mean of reference cells.
       c. Clip to [-3, 3] to reduce outlier influence.
    4. Return matrix: cells × windows.

    SIGNIFICANCE:
    The per-window normalization by reference cells removes constitutively
    expressed gene clusters that could mimic CNV signal (e.g., HOX gene
    clusters).  The resulting matrix encodes amplifications as positive
    values and deletions as negative values, interpretable as relative
    copy number.
    """
    # Intersect genes with known positions
    common_genes = positions.index.intersection(adata.var_names)
    log.info(f"  {len(common_genes):,} genes with chromosomal positions")
    if len(common_genes) < 1000:
        log.warning(
            f"  Only {len(common_genes):,} positioned genes — CNV resolution "
            "may be limited. Ensure adata contains protein-coding genes."
        )

    # Ordered gene expression matrix (cells × positioned genes)
    pos_ordered = positions.loc[common_genes].sort_values(["chromosome", "start"])
    gene_order  = pos_ordered.index.tolist()

    if "scvi_normalized" in adata.layers:
        expr = adata[:, gene_order].layers["scvi_normalized"]
    else:
        log.warning("  'scvi_normalized' layer not found; using .X (log-normalized)")
        expr = adata[:, gene_order].X

    # Convert sparse to dense if needed
    if hasattr(expr, "toarray"):
        expr = expr.toarray()

    n_cells, n_genes = expr.shape
    n_windows = n_genes - window_size + 1
    cnv = np.zeros((n_cells, n_windows), dtype=np.float32)

    ref_mean_total = expr[reference_mask].mean(axis=0)  # per-gene reference

    for w in range(n_windows):
        window_expr = expr[:, w : w + window_size].mean(axis=1)   # per-cell
        ref_window  = ref_mean_total[w : w + window_size].mean()  # reference
        cnv[:, w] = window_expr - ref_window

    cnv = np.clip(cnv, -3, 3)
    log.info(f"  CNV matrix: {n_cells:,} cells × {n_windows:,} windows")
    return cnv, pos_ordered, gene_order, n_windows


def cnv_summary_score(cnv_matrix: np.ndarray) -> np.ndarray:
    """Compute per-cell CNV load as mean absolute deviation across windows.

    SIGNIFICANCE:
    A single scalar per cell summarises the overall chromosomal instability
    level: tumour cells with multiple SCNAs will score high; normal cells
    will score near zero.  This score can be overlaid on UMAP to visually
    confirm that annotated tumour clusters have elevated genomic instability.
    """
    return np.mean(np.abs(cnv_matrix), axis=1)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_cnv_heatmap(
    cnv_matrix: np.ndarray,
    cell_types: pd.Series,
    pos_ordered: pd.DataFrame,
    out_path: Path,
    n_sample: int = 2000,
) -> None:
    """Heatmap of CNV signal across chromosomal windows, ordered by cell type.

    SIGNIFICANCE:
    This is the primary CNV QC plot.  Tumour cells should show coherent
    column-wise blocks of red (amplification) or blue (deletion) at the
    expected chromosomal coordinates.  Normal cells should appear grey.
    Rows are cells, columns are chromosomal windows ordered by position.
    """
    # Sub-sample for plot speed
    idx = np.random.choice(len(cell_types), size=min(n_sample, len(cell_types)), replace=False)
    idx = idx[np.argsort(cell_types.iloc[idx].values)]  # sort by cell type

    # Chromosome tick positions
    chrom_changes = pos_ordered["chromosome"].values
    chrom_ticks   = []
    chrom_labels  = []
    prev_chrom    = None
    for i, c in enumerate(chrom_changes):
        if c != prev_chrom:
            chrom_ticks.append(i)
            chrom_labels.append(str(c))
            prev_chrom = c

    fig, ax = plt.subplots(figsize=(18, 6))
    im = ax.imshow(
        cnv_matrix[idx],
        aspect="auto",
        cmap="RdBu_r",
        vmin=-1.5,
        vmax=1.5,
        interpolation="none",
    )
    ax.set_xticks(chrom_ticks)
    ax.set_xticklabels(chrom_labels, fontsize=7, rotation=0)
    ax.set_xlabel("Chromosomal position (windows)", fontsize=10)
    ax.set_ylabel(f"Cells (n={len(idx)}, ordered by type)", fontsize=10)
    ax.set_title(
        "Inferred Copy Number Variation\n"
        "(Red = amplification, Blue = deletion, Grey = neutral)",
        fontsize=11, fontweight="bold",
    )
    plt.colorbar(im, ax=ax, shrink=0.6, label="CNV score")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_cnv_umap(adata: sc.AnnData, out_path: Path) -> None:
    """UMAP coloured by per-cell CNV load score.

    SIGNIFICANCE:
    Overlaying CNV load on UMAP allows visual confirmation that annotated
    tumour clusters (Cone_precursor, Tumour_proliferating) have higher
    chromosomal instability than normal cell clusters.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(
        adata, color="cnv_load", ax=axes[0], show=False,
        frameon=False, size=2, cmap="YlOrRd",
        title="CNV load score",
        vmin=0, vmax=adata.obs["cnv_load"].quantile(0.95),
    )
    sc.pl.umap(
        adata, color="cell_type_broad", ax=axes[1], show=False,
        frameon=False, size=2,
        title="Cell type (broad)",
        legend_loc="right margin",
    )
    plt.suptitle("CNV load vs. cell type annotation", fontsize=12,
                 fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_cnv() -> sc.AnnData:
    """CNV inference pipeline.

    Steps
    -----
    1. Load annotated atlas.
    2. Fetch gene chromosomal positions.
    3. Identify reference (non-malignant) cells.
    4. Compute CNV matrix via sliding chromosomal windows.
    5. Compute per-cell CNV load score.
    6. Store results in adata, plot, save.
    7. Call tumour vs. normal status based on CNV load threshold.
    """

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading annotated atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    log.info(f"  Cell types: {adata.obs['cell_type_broad'].value_counts().to_dict()}")

    # ---- 2. Gene positions ----------------------------------------------
    log.info("\nSTEP 2 — Fetching gene chromosomal positions")
    log.info(
        "  Chromosomal positions are required to order genes for CNV\n"
        "  inference. Genes are sorted by (chromosome, start_position)\n"
        "  before windowed averaging."
    )
    positions = get_gene_positions(adata)
    log.info(f"  Positioned genes across {positions['chromosome'].nunique()} chromosomes")

    # ---- 3. Reference cells ---------------------------------------------
    log.info("\nSTEP 3 — Identifying reference (non-malignant) cells")
    log.info(
        f"  Reference cell types: {REFERENCE_CELL_TYPES}\n"
        "  WHY: Reference cells are assumed diploid. Subtracting their\n"
        "  per-window mean expression cancels gene-cluster expression\n"
        "  patterns that are not copy-number-driven."
    )
    ref_mask = adata.obs["cell_type_broad"].isin(REFERENCE_CELL_TYPES).values
    n_ref = ref_mask.sum()
    if n_ref < 100:
        log.warning(
            f"  Only {n_ref} reference cells found. CNV normalisation\n"
            "  may be unreliable. Consider using more reference types."
        )
    else:
        log.info(f"  Using {n_ref:,} reference cells")

    # ---- 4. Compute CNV matrix ------------------------------------------
    log.info("\nSTEP 4 — Computing CNV matrix (sliding chromosomal windows)")
    log.info(
        f"  Window size: {WINDOW_SIZE} genes\n"
        "  This is the standard inferCNV window (Tirosh et al. 2016).\n"
        "  Larger windows increase signal-to-noise but reduce resolution.\n"
        "  Each window represents ~1-3 Mb depending on gene density."
    )
    cnv_matrix, pos_ordered, gene_order, n_windows = compute_cnv_matrix(
        adata, positions, ref_mask
    )

    # ---- 5. Per-cell CNV load score ------------------------------------
    log.info("\nSTEP 5 — Computing per-cell CNV load score")
    log.info(
        "  CNV load = mean |CNV score| across all windows.\n"
        "  Higher score = more chromosomal instability = likely tumour cell.\n"
        "  Threshold for tumour call: > mean_reference + 3 × SD_reference"
    )
    cnv_load = cnv_summary_score(cnv_matrix)
    adata.obs["cnv_load"] = cnv_load

    # Threshold based on reference cell distribution
    ref_load   = cnv_load[ref_mask]
    threshold  = ref_load.mean() + 3 * ref_load.std()
    adata.obs["is_tumour_cnv"] = (cnv_load > threshold).astype(str)
    adata.obs["is_tumour_cnv"] = adata.obs["is_tumour_cnv"].map(
        {"True": "Tumour", "False": "Normal/uncertain"}
    )
    n_tumour = (adata.obs["is_tumour_cnv"] == "Tumour").sum()
    log.info(
        f"  CNV threshold: {threshold:.4f}\n"
        f"  Cells called as tumour: {n_tumour:,} / {adata.n_obs:,} "
        f"({100 * n_tumour / adata.n_obs:.1f}%)"
    )

    # Store CNV matrix in obsm
    adata.obsm["X_cnv"] = cnv_matrix.astype(np.float32)

    # ---- 6. Region-specific CNV scores ----------------------------------
    log.info("\nSTEP 6 — Computing region-specific RB SCNA scores")
    log.info(
        "  Scoring key RB chromosomal regions:\n"
        "    chr2p  : MYCN amplification (Subtype 2)\n"
        "    chr6p  : E2F3/DEK amplification (aggressive)\n"
        "    chr13q : RB1 deletion (initiating event)\n"
        "    chr16q : loss (metastatic)\n"
        "    chr1q  : gain (MDM4 / p53 pathway)"
    )
    # Map windows back to positions for region scoring
    window_chroms = pos_ordered["chromosome"].values[:n_windows].astype(int)
    window_starts = pos_ordered["start"].values[:n_windows]

    for region_name, (chrom, start, end) in RB_SCNA_REGIONS.items():
        chrom_int = int(chrom)
        win_mask  = (
            (window_chroms == chrom_int)
            & (window_starts >= start)
            & (window_starts <  end)
        )
        if win_mask.sum() > 0:
            region_score = cnv_matrix[:, win_mask].mean(axis=1)
            adata.obs[f"cnv_{region_name}"] = region_score
            log.info(
                f"  {region_name}: {win_mask.sum()} windows, "
                f"mean score = {region_score.mean():.4f}"
            )
        else:
            log.warning(f"  No windows found for region {region_name}")

    # ---- 7. Plots -------------------------------------------------------
    log.info("\nSTEP 7 — Generating CNV visualizations")
    plot_cnv_heatmap(
        cnv_matrix, adata.obs["cell_type_broad"],
        pos_ordered,
        FIG_DIR / "cnv_heatmap_chromosome.pdf",
    )
    plot_cnv_umap(adata, FIG_DIR / "cnv_umap_score.pdf")

    # Save CNV scores table
    cnv_cols = ["cnv_load", "is_tumour_cnv"] + [
        f"cnv_{r}" for r in RB_SCNA_REGIONS if f"cnv_{r}" in adata.obs.columns
    ]
    adata.obs[cnv_cols].to_csv(TAB_DIR / "cnv_scores_per_cell.csv")
    log.info(f"  Saved CNV scores → {TAB_DIR / 'cnv_scores_per_cell.csv'}")

    # Summarize by cell type
    cnv_summary = adata.obs.groupby("cell_type_broad")["cnv_load"].agg(
        ["mean", "median", "std", "count"]
    ).round(4)
    cnv_summary.to_csv(TAB_DIR / "cnv_load_by_celltype.csv")
    log.info(f"  CNV load by cell type:\n{cnv_summary.to_string()}")

    # ---- 8. Save --------------------------------------------------------
    log.info(f"\nSTEP 8 — Saving CNV-annotated atlas → {OUT_H5AD.name}")
    log.info(
        "  Contents added:\n"
        "    .obs['cnv_load']       : per-cell CNV instability score\n"
        "    .obs['is_tumour_cnv']  : Tumour / Normal call\n"
        "    .obs['cnv_chr*']       : region-specific CNV scores\n"
        "    .obsm['X_cnv']         : full CNV matrix (cells × windows)"
    )
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    return adata


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_cnv()
