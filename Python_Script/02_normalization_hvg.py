"""
=============================================================================
Script 02 — Normalization, Log-Transformation, and HVG Selection
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
Raw UMI counts are dominated by library-size effects: a cell sequenced more
deeply will have higher counts for every gene, independent of true biology.
Normalization removes this technical covariate so that expression values
reflect actual gene activity rather than sequencing luck.

The three-step normalization here (normalize → log1p → scale) transforms
the data into a space where:
  - Shallow and deep cells contribute equally to downstream analyses
  - Log-transformation compresses extreme values and makes the data closer
    to normally distributed (required for PCA, which underpins UMAP and
    clustering)
  - Scaling (zero-mean, unit-variance) ensures that highly expressed "house-
    keeping" genes don't dominate the principal components

SIGNIFICANCE OF HVG SELECTION
------------------------------
The human genome encodes ~20,000 protein-coding genes.  For this ~100,000-cell
dataset, the full gene × cell matrix would require enormous memory and
introduce noise from genes that do not vary meaningfully across cell types.
Selecting the top 3,000 Highly Variable Genes (HVGs) captures genes whose
expression actually distinguishes cell types and states, while reducing
dimensionality ~7-fold.

In retinoblastoma specifically, HVG selection captures cone-specific markers
(ARR3, RXRG, THRB), proliferation markers (MKI67, TOP2A), and immune markers
(CD68, IBA1) that are central to the atlas cell-type annotation (Aim 2) and
the cone precursor trajectory (Aim 3).

NOTE ON BATCH-AWARE HVG SELECTION
----------------------------------
Standard HVG selection on the merged matrix can be dominated by between-
batch technical variability.  We use Seurat v3-flavoured HVG selection
(vstfeed with per-batch variance estimation via `batch_key='sample_id'`) so
that genes selected are variable *within* samples, not just between them.
This is the recommended approach before scVI integration (Lopez et al., 2018,
Nat Methods).

REFERENCES
----------
- Hafemeister C, Satija R. Normalization and variance stabilization of single-
  cell RNA-seq data using regularized negative binomial regression.
  Genome Biol. 2019;20:296. https://doi.org/10.1186/s13059-019-1874-1

- Lopez R et al. Deep generative modeling for single-cell transcriptomics.
  Nat Methods. 2018;15(12):1053-1058.
  https://doi.org/10.1038/s41592-018-0229-2

INPUT  : data/processed/02_qc_filtered.h5ad
OUTPUT :
  data/processed/03_normalized.h5ad
  results/figures/hvg_mean_variance.pdf
  results/figures/pca_variance_explained.pdf
  results/tables/highly_variable_genes.csv

USAGE  : pixi run python script/02_normalization_hvg.py
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
import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data" / "processed" / "02_qc_filtered.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "03_normalized.h5ad"
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
# Configuration
# ---------------------------------------------------------------------------
N_HVG        = 3_000      # number of highly variable genes to select
N_PCS        = 50         # PCs for pre-scVI neighbour graph (diagnostic)
TARGET_SUM   = 10_000     # library-size normalization target (CPM-like)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_hvg_mean_var(adata: sc.AnnData, out_path: Path) -> None:
    """Mean-variance plot coloured by HVG status.

    SIGNIFICANCE:
    The mean-variance relationship in raw count data follows an overdispersed
    Poisson (negative binomial) distribution. Seurat v3 HVG selection fits a
    regularized model to this relationship and selects genes that deviate
    upward (higher variance than expected for their mean expression level).
    These are the genes that are genuinely differentially expressed across
    cell types rather than stochastically noisy low-count genes.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    not_hvg = ~adata.var["highly_variable"]
    ax.scatter(
        adata.var.loc[not_hvg, "means"],
        adata.var.loc[not_hvg, "dispersions_norm"],
        s=2, alpha=0.3, color="#AAAAAA", label="Not HVG", rasterized=True,
    )
    ax.scatter(
        adata.var.loc[adata.var["highly_variable"], "means"],
        adata.var.loc[adata.var["highly_variable"], "dispersions_norm"],
        s=3, alpha=0.7, color="#E84C4C", label=f"HVG (n={adata.var['highly_variable'].sum():,})",
        rasterized=True,
    )
    ax.set_xscale("log")
    ax.set_xlabel("Mean expression", fontsize=11)
    ax.set_ylabel("Normalised dispersion", fontsize=11)
    ax.set_title(
        f"Highly Variable Gene selection\n"
        f"(Seurat v3, batch_key='sample_id', top {N_HVG:,} genes)",
        fontsize=10,
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_pca_variance(adata: sc.AnnData, out_path: Path) -> None:
    """Scree plot of PCA variance explained.

    SIGNIFICANCE:
    The elbow in the variance-explained curve is used to choose the number
    of principal components for the neighbour graph.  Including too many PCs
    adds noise from low-variance components; too few loses biological signal.
    For retinoblastoma, the first 20-30 PCs typically capture the major
    cell-type axes (tumour vs. immune vs. glial) while later PCs capture
    within-type variation relevant to trajectory analysis.

    NOTE: This diagnostic PCA is run on the HVG-scaled matrix BEFORE scVI
    batch correction and is used only for quality assessment, not for the
    final analysis (which uses the scVI latent space).
    """
    var_ratio = adata.uns["pca"]["variance_ratio"]
    cum_var   = np.cumsum(var_ratio)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].bar(range(1, len(var_ratio) + 1), var_ratio * 100,
                color="#4C9BE8", edgecolor="none")
    axes[0].set_xlabel("Principal component", fontsize=11)
    axes[0].set_ylabel("Variance explained (%)", fontsize=11)
    axes[0].set_title("Scree plot", fontsize=11)

    axes[1].plot(range(1, len(cum_var) + 1), cum_var * 100,
                 color="#E84C4C", lw=2)
    axes[1].axhline(80, color="grey", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Number of PCs", fontsize=11)
    axes[1].set_ylabel("Cumulative variance (%)", fontsize=11)
    axes[1].set_title("Cumulative variance explained", fontsize=11)
    pc80 = int(np.searchsorted(cum_var, 0.8)) + 1
    axes[1].annotate(f"80% at PC{pc80}", xy=(pc80, 80),
                     xytext=(pc80 + 2, 82), fontsize=9,
                     arrowprops=dict(arrowstyle="->"))

    plt.suptitle("PCA diagnostics (HVG-scaled matrix, pre-scVI)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_normalization() -> sc.AnnData:
    """Normalization and HVG selection pipeline.

    Steps
    -----
    1. Load QC-filtered atlas.
    2. Normalize to fixed library size (10,000 counts / cell).
       WHY: Makes expression values comparable across cells with different
       sequencing depths.  10,000 is convention (approximately CPM) and
       produces values in a convenient numeric range.
    3. Log1p transform.
       WHY: Compresses the dynamic range of expression values.
       log(1 + x) is preferred over log(x) to handle zero counts without
       undefined values.  The +1 pseudo-count also makes the transformed
       distribution approximately symmetric for downstream linear methods
       (PCA, regression).
    4. Batch-aware HVG selection.
       WHY: Selects genes that are variable within samples, not just between
       them, reducing the influence of batch-specific artefacts on the HVG
       set used for initial clustering (pre-scVI).
    5. Scale to unit variance (clip at max_value=10 to limit outlier influence).
       WHY: PCA is a variance-maximising projection; without scaling,
       highly expressed genes (e.g., ACTB, GAPDH) dominate the first
       principal components regardless of their biological relevance.
    6. Run diagnostic PCA on HVG-scaled matrix.
    7. Save processed AnnData.
    """

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading QC-filtered atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 2. Normalize ---------------------------------------------------
    log.info("\nSTEP 2 — Library-size normalization (target = 10,000 counts)")
    log.info(
        "  WHY: Cells with more total counts are deeper-sequenced, not\n"
        "  more transcriptionally active. Normalization removes this\n"
        "  library-size technical covariate."
    )
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    log.info(f"  Normalized. Median total counts = {adata.X.sum(axis=1).A1.median():.0f}" if hasattr(adata.X, 'A1') else "  Normalized.")

    # ---- 3. Log1p -------------------------------------------------------
    log.info("\nSTEP 3 — log1p transformation")
    log.info(
        "  WHY: scRNA-seq data is highly right-skewed (a few genes are\n"
        "  extremely highly expressed). log(1+x) compresses extreme values\n"
        "  and makes count distributions closer to Gaussian, which is\n"
        "  assumed by PCA and other downstream methods."
    )
    sc.pp.log1p(adata)
    # Store lognorm layer for reference
    adata.layers["lognorm"] = adata.X.copy()

    # ---- 4. Batch-aware HVG selection -----------------------------------
    log.info(f"\nSTEP 4 — Highly Variable Gene (HVG) selection")
    log.info(
        f"  Selecting top {N_HVG:,} HVGs using Seurat v3 method\n"
        "  with batch_key='sample_id' (variance estimated within each\n"
        "  sample to avoid batch-driven HVG inflation)."
    )
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_HVG,
        flavor="seurat_v3",
        layer="counts",            # use raw counts for dispersion fitting
        batch_key="sample_id",     # batch-aware selection
        subset=False,              # keep all genes; HVG flag stored in .var
    )
    n_hvg = adata.var["highly_variable"].sum()
    log.info(f"  Selected {n_hvg:,} HVGs")

    # Save HVG table
    hvg_df = adata.var[
        ["highly_variable", "means", "dispersions_norm", "highly_variable_rank"]
    ].query("highly_variable").sort_values("highly_variable_rank")
    hvg_df.to_csv(TAB_DIR / "highly_variable_genes.csv")
    log.info(f"  Saved HVG list → {TAB_DIR / 'highly_variable_genes.csv'}")
    plot_hvg_mean_var(adata, FIG_DIR / "hvg_mean_variance.pdf")

    # ---- 5. Scale (HVG-only copy for PCA) --------------------------------
    log.info("\nSTEP 5 — Scaling HVG matrix for PCA diagnostic")
    log.info(
        "  WHY: Zero-mean, unit-variance scaling ensures all HVGs\n"
        "  contribute equally to PCA regardless of their expression level.\n"
        "  max_value=10 caps outliers (e.g., extreme MALAT1 expression)\n"
        "  that would otherwise distort principal components."
    )
    # Work on HVG subset for PCA
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)

    # ---- 6. Diagnostic PCA ----------------------------------------------
    log.info(f"\nSTEP 6 — Diagnostic PCA ({N_PCS} PCs)")
    log.info(
        "  NOTE: This PCA is computed on the raw merged matrix (before scVI\n"
        "  batch correction) for diagnostic purposes only. The actual\n"
        "  dimensionality reduction used in the atlas will use the scVI\n"
        "  latent space (script 03)."
    )
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver="arpack")
    # Copy PCA back into main adata for reference
    adata.obsm["X_pca_pre_scvi"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"]              = adata_hvg.uns["pca"]
    plot_pca_variance(adata_hvg, FIG_DIR / "pca_variance_explained.pdf")

    # ---- 7. Save --------------------------------------------------------
    log.info(f"\nSTEP 7 — Saving normalized AnnData → {OUT_H5AD.name}")
    log.info(
        "  Contents:\n"
        "    .X              : log-normalized counts (float32)\n"
        "    .layers['counts']: raw UMI counts (int) — preserved for scVI\n"
        "    .layers['lognorm']: log-normalized copy\n"
        "    .var['highly_variable']: HVG flag\n"
        "    .obsm['X_pca_pre_scvi']: pre-correction PCA (diagnostic)"
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
    run_normalization()
