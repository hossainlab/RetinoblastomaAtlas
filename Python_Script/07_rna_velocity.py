"""
=============================================================================
Script 07 — RNA Velocity Analysis (scVelo Dynamical Model)
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
RNA velocity estimates the *future* transcriptional state of each cell by
analysing the ratio of unspliced (nascent) to spliced (mature) mRNA.
Unspliced pre-mRNA is rapidly produced after gene activation; spliced mRNA
accumulates more slowly and decays faster once a gene is silenced.  The
imbalance between the two forms reveals the *direction and speed* of
transcriptional change.

In the context of this project, RNA velocity allows us to:

  1. Determine the DIRECTIONALITY of the cone precursor trajectory:
     Are tumour cells transitioning TOWARD or AWAY FROM the mature cone
     state?  A velocity field pointing from early progenitor-like states
     toward mature cone states would imply differentiation; one pointing in
     the opposite direction would imply dedifferentiation — consistent with
     the Subtype 1 → Subtype 2 invasion model.

  2. Identify TRANSITIONAL STATES in the invasion trajectory:
     Cells at bifurcation points (where velocity vectors diverge) represent
     critical decision points where microenvironmental signals (TGF-β,
     macrophage interactions) may push cells toward the invasive fate.

  3. Estimate KINETIC PARAMETERS per gene:
     The dynamical model (Bergen et al. 2020) fits four parameters for each
     gene: transcription rate (α), splicing rate (β), degradation rate (γ),
     and switching time (t_s). Genes with high induction rate and late
     switch-off are candidate invasion driver genes.

WHY THE DYNAMICAL MODEL OVER THE STEADY-STATE MODEL?
------------------------------------------------------
The original steady-state RNA velocity model (La Manno et al. 2018) assumes
that all cells are at steady state (spliced/unspliced ratio is constant).
This assumption is violated in tumour evolution where cells are actively
transitioning.  The scVelo dynamical model (Bergen et al. 2020) relaxes
this assumption by inferring the full transcriptional kinetics for each
gene — making it the method of choice for developmental or oncological
trajectories.

INPUT REQUIREMENTS
------------------
RNA velocity requires BOTH spliced and unspliced count matrices.  These are
generated during alignment (STARsolo or Cell Ranger with --include-introns).
The matrices are stored in separate layers in the input AnnData:
  .layers['spliced']   — mature mRNA counts (aligned to exons)
  .layers['unspliced'] — nascent pre-mRNA counts (aligned to introns)

If these layers are absent, the pipeline will attempt to load them from
matching loom files produced by the Velocyto CLI.

REFERENCES
----------
- Bergen V et al. Generalizing RNA velocity to transient cell states through
  dynamical modeling. Nat Biotechnol. 2020;38(12):1408-1414.
  https://doi.org/10.1038/s41587-020-0591-3

- La Manno G et al. RNA velocity of single cells. Nature. 2018;560:494-498.
  https://doi.org/10.1038/s41586-018-0414-6

- Liu Y et al. Single-cell transcriptomic analysis reveals local extension
  mechanism of retinoblastoma. Cell Mol Life Sci. 2024;81:77.
  https://doi.org/10.1007/s00018-024-05102-9

INPUT  : data/processed/07_subtyped.h5ad
         (optional) data/raw/loom/*.loom — if spliced/unspliced absent
OUTPUT :
  data/processed/08_velocity.h5ad
  results/figures/velocity_stream_umap.pdf
  results/figures/velocity_confidence_umap.pdf
  results/figures/velocity_top_driver_genes.pdf
  results/tables/velocity_driver_genes.csv

USAGE  : pixi run python script/07_rna_velocity.py [--n-top-genes 2000]
=============================================================================
"""

from __future__ import annotations

import argparse
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
import scvelo as scv

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT      = Path(__file__).resolve().parents[1]
IN_H5AD   = ROOT / "data" / "processed" / "07_subtyped.h5ad"
OUT_H5AD  = ROOT / "data" / "processed" / "08_velocity.h5ad"
LOOM_DIR  = ROOT / "data" / "raw" / "loom"
FIG_DIR   = ROOT / "results" / "figures"
TAB_DIR   = ROOT / "results" / "tables"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)
warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
scv.settings.verbosity = 1
scv.settings.presenter_view = False


# ---------------------------------------------------------------------------
# Loom loading helper
# ---------------------------------------------------------------------------

def load_loom_velocyto(adata: sc.AnnData) -> sc.AnnData:
    """Load spliced/unspliced matrices from Velocyto loom files.

    SIGNIFICANCE:
    Velocyto (La Manno et al. 2018) processes BAM files to generate loom
    files with per-cell spliced and unspliced UMI counts.  When the atlas
    h5ad does not already contain these layers (common when the h5ad was
    created from Cell Ranger output without --include-introns), we merge
    in the loom data.

    Expected loom file names follow the pattern:
      <sample_id>.loom
    where sample_id values match adata.obs['sample_id'].
    """
    loom_files = list(LOOM_DIR.glob("*.loom"))
    if not loom_files:
        raise FileNotFoundError(
            f"No .loom files found in {LOOM_DIR}.\n"
            "  Please run Velocyto on the original BAM files:\n"
            "    velocyto run10x <cellranger_output_dir> <genome_gtf>\n"
            "  or use STARsolo with --soloFeatures Velocyto."
        )

    loom_adatas = []
    for loom_path in sorted(loom_files):
        loom_adata = scv.read(str(loom_path), cache=True)
        loom_adata.var_names_make_unique()
        loom_adatas.append(loom_adata)
        log.info(f"  Loaded loom: {loom_path.name} ({loom_adata.n_obs:,} cells)")

    # Merge loom data into main adata
    # scVelo's merge() aligns barcodes and transfers spliced/unspliced layers
    adata_merged = scv.utils.merge(adata, loom_adatas[0] if len(loom_adatas) == 1
                                   else sc.concat(loom_adatas))
    log.info(
        f"  After loom merge: {adata_merged.n_obs:,} cells "
        f"(some barcodes may be lost if loom/h5ad barcodes don't match)"
    )
    return adata_merged


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_velocity_stream(
    adata: sc.AnnData,
    color_keys: list[str],
    title: str,
    out_path: Path,
) -> None:
    """UMAP with RNA velocity stream plot."""
    ncols = len(color_keys)
    fig, axes = plt.subplots(1, ncols, figsize=(7 * ncols, 6))
    if ncols == 1:
        axes = [axes]
    for ax, key in zip(axes, color_keys):
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color=key,
            ax=ax, show=False, frameon=False,
            size=20, alpha=0.7,
            arrow_size=1.5, linewidth=0.8,
            title=key,
        )
    plt.suptitle(title, fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_velocity_confidence(adata: sc.AnnData, out_path: Path) -> None:
    """UMAP coloured by velocity confidence and length (speed).

    SIGNIFICANCE:
    Velocity confidence (coherence of velocity vectors in local neighbourhood)
    is a quality metric: cells in well-defined trajectories should have high
    confidence.  Low-confidence cells may be in ambiguous transcriptional
    states or represent noise.  Velocity length encodes transcriptional
    speed — faster-changing cells are often closer to fate decisions.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    scv.pl.umap(adata, color="velocity_confidence", ax=axes[0], show=False,
                frameon=False, cmap="RdYlGn", vmin=0, vmax=1,
                title="Velocity confidence")
    scv.pl.umap(adata, color="velocity_length", ax=axes[1], show=False,
                frameon=False, cmap="YlOrRd",
                title="Velocity length (transcriptional speed)")
    plt.suptitle("RNA velocity quality metrics", fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_top_driver_genes(adata: sc.AnnData, out_path: Path, n_top: int = 12) -> None:
    """Phase portraits of top velocity-driver genes.

    SIGNIFICANCE:
    Phase portraits show the spliced vs. unspliced expression of individual
    genes across cells.  In the dynamical model, the fitted kinetic curve
    (black line) represents the transcriptional cycle.  Genes that are
    induction drivers have most cells in the upper-left (unspliced rising)
    portion of the curve; repression drivers have cells in the lower-right
    (spliced decaying) portion.  Genes with clear phase portrait structure
    and high likelihood are the most reliable velocity genes.
    """
    # Get top likelihood genes
    if "fit_likelihood" not in adata.var.columns:
        log.warning("  'fit_likelihood' not in adata.var — run dynamical model first")
        return
    top_genes = (
        adata.var["fit_likelihood"]
        .dropna()
        .sort_values(ascending=False)
        .head(n_top)
        .index.tolist()
    )
    log.info(f"  Top {n_top} driver genes: {top_genes[:6]}...")

    ncols = 4
    nrows = int(np.ceil(n_top / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten()

    for ax, gene in zip(axes, top_genes):
        scv.pl.velocity(adata, var_names=[gene], ax=ax, show=False)
        ax.set_title(gene, fontsize=9)
    for ax in axes[n_top:]:
        ax.set_visible(False)

    plt.suptitle(
        f"Top {n_top} velocity driver genes (sorted by fit likelihood)",
        fontsize=11, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_velocity(
    n_top_genes: int = 2000,
    min_shared_counts: int = 30,
    subset_cell_type: str = "Cone_precursor",
) -> sc.AnnData:
    """RNA velocity pipeline using scVelo dynamical model.

    Steps
    -----
    1. Load subtyped atlas.
    2. Ensure spliced/unspliced layers are present (load from loom if needed).
    3. Filter and normalise for velocity analysis.
    4. Select top velocity genes.
    5. Compute moments (first- and second-order moments for EM algorithm).
    6. Fit dynamical model (recover full kinetics).
    7. Compute velocity vectors and project onto UMAP.
    8. Calculate velocity confidence and pseudotime.
    9. Identify top driver genes by fit likelihood.
    10. Optionally run velocity on Cone_precursor subset only.
    11. Save.
    """

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading subtyped atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 2. Ensure spliced/unspliced layers -----------------------------
    log.info("\nSTEP 2 — Checking for spliced/unspliced layers")
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        log.warning(
            "  Spliced/unspliced layers not found in h5ad.\n"
            "  Attempting to load from Velocyto loom files..."
        )
        adata = load_loom_velocyto(adata)
    else:
        n_s = (adata.layers["spliced"] > 0).sum()
        n_u = (adata.layers["unspliced"] > 0).sum()
        log.info(
            f"  Spliced layer: {n_s:,} non-zero entries\n"
            f"  Unspliced layer: {n_u:,} non-zero entries"
        )

    # ---- 3. Preprocessing for velocity ----------------------------------
    log.info("\nSTEP 3 — Preprocessing for velocity")
    log.info(
        "  scVelo requires both spliced and unspliced counts.\n"
        "  Cells with fewer than min_shared_counts are removed.\n"
        "  WHY min_shared_counts: Low-count cells have unreliable\n"
        "  spliced/unspliced ratios due to sampling noise."
    )
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=min_shared_counts,
        n_top_genes=n_top_genes,
        log=True,
    )
    log.info(f"  After filter_and_normalize: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 4. Moments ------------------------------------------------------
    log.info("\nSTEP 4 — Computing first/second order moments")
    log.info(
        "  Moments are computed in the scVI latent space (30-d KNN graph).\n"
        "  WHY: Using the scVI-corrected latent space rather than raw PCA\n"
        "  for KNN graph construction ensures that velocity neighbours are\n"
        "  biologically similar rather than batch-confounded."
    )
    use_rep = "X_scVI" if "X_scVI" in adata.obsm else "X_pca"
    scv.pp.moments(
        adata,
        n_pcs=30 if use_rep == "X_pca" else None,
        n_neighbors=30,
        use_rep=use_rep,
    )

    # ---- 5. Dynamical model (full kinetics) ----------------------------
    log.info(f"\nSTEP 5 — Fitting dynamical model (n_top_genes={n_top_genes})")
    log.info(
        "  The dynamical model (Bergen et al. 2020) estimates per-gene:\n"
        "    α : transcription rate\n"
        "    β : splicing rate\n"
        "    γ : degradation rate\n"
        "    t_s: switching time (induction → repression)\n"
        "  This is compute-intensive (~minutes for 100k cells).\n"
        "  Early stopping patience=100 avoids over-fitting noisy genes."
    )
    scv.tl.recover_dynamics(adata, n_jobs=4)
    scv.tl.velocity(adata, mode="dynamical")

    # ---- 6. Project velocity to UMAP -----------------------------------
    log.info("\nSTEP 6 — Projecting velocity onto UMAP")
    log.info(
        "  scv.tl.velocity_graph() builds a cell-transition probability\n"
        "  matrix based on the cosine similarity between each cell's\n"
        "  velocity vector and its neighbours' displacement vectors.\n"
        "  This graph is used for velocity embedding and pseudotime."
    )
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="umap")

    # ---- 7. Velocity confidence and length ----------------------------
    log.info("\nSTEP 7 — Computing velocity confidence and pseudotime")
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_pseudotime(adata)
    log.info(
        "  Velocity pseudotime assigns each cell a position along the\n"
        "  most probable RNA velocity trajectory (0 = root, 1 = tip).\n"
        "  This is used in CellRank (script 08) for fate mapping."
    )

    # ---- 8. Driver genes -----------------------------------------------
    log.info("\nSTEP 8 — Identifying top driver genes by fit likelihood")
    log.info(
        "  Fit likelihood quantifies how well the dynamical model fits\n"
        "  the phase portrait for each gene.  Genes with high likelihood\n"
        "  have coherent spliced/unspliced dynamics and are reliable\n"
        "  velocity drivers.  These are candidate genes regulating\n"
        "  the invasion transition."
    )
    if "fit_likelihood" in adata.var.columns:
        driver_df = (
            adata.var[["fit_likelihood", "fit_alpha", "fit_beta", "fit_gamma",
                        "fit_t_", "gene_count_corr"]]
            .dropna()
            .sort_values("fit_likelihood", ascending=False)
            .head(200)
        )
        driver_df.to_csv(TAB_DIR / "velocity_driver_genes.csv")
        log.info(
            f"  Top 5 driver genes: "
            f"{driver_df.head(5).index.tolist()}\n"
            f"  Saved driver gene table → {TAB_DIR / 'velocity_driver_genes.csv'}"
        )

    # ---- 9. Plots ------------------------------------------------------
    log.info("\nSTEP 9 — Generating velocity visualizations")
    color_keys = ["cell_type_broad", "rb_subtype", "velocity_pseudotime"]
    color_keys = [k for k in color_keys if k in adata.obs.columns]
    plot_velocity_stream(
        adata, color_keys,
        "RNA velocity stream (scVelo dynamical model)",
        FIG_DIR / "velocity_stream_umap.pdf",
    )
    plot_velocity_confidence(adata, FIG_DIR / "velocity_confidence_umap.pdf")
    plot_top_driver_genes(adata, FIG_DIR / "velocity_top_driver_genes.pdf")

    # ---- 10. Cone precursor subset velocity ---------------------------
    log.info(f"\nSTEP 10 — Running velocity on {subset_cell_type} subset only")
    log.info(
        "  Restricting to Cone_precursor cells for high-resolution\n"
        "  trajectory analysis.  In the full atlas, velocity arrows are\n"
        "  dominated by immune-to-tumour transitions; the subset view\n"
        "  reveals fine within-tumour directionality."
    )
    if subset_cell_type in adata.obs["cell_type_broad"].values:
        adata_sub = adata[adata.obs["cell_type_broad"] == subset_cell_type].copy()
        scv.pp.moments(adata_sub, n_pcs=20, n_neighbors=20)
        scv.tl.velocity(adata_sub, mode="dynamical")
        scv.tl.velocity_graph(adata_sub)
        scv.tl.velocity_embedding(adata_sub, basis="umap")
        plot_velocity_stream(
            adata_sub,
            ["cell_type_fine", "velocity_pseudotime"],
            f"RNA velocity — {subset_cell_type} subset",
            FIG_DIR / f"velocity_stream_{subset_cell_type.lower()}_subset.pdf",
        )
        # Store subset pseudotime back in full adata
        adata.obs.loc[adata_sub.obs_names, f"velocity_pseudotime_{subset_cell_type}"] = \
            adata_sub.obs["velocity_pseudotime"]
    else:
        log.warning(f"  No cells of type '{subset_cell_type}' found — skipping subset velocity")

    # ---- 11. Save -------------------------------------------------------
    log.info(f"\nSTEP 11 — Saving velocity-annotated atlas → {OUT_H5AD.name}")
    log.info(
        "  Contents added:\n"
        "    .layers['velocity']         : velocity vectors (spliced space)\n"
        "    .obs['velocity_pseudotime'] : RNA velocity pseudotime\n"
        "    .obs['velocity_confidence'] : per-cell velocity confidence\n"
        "    .obsm['velocity_umap']      : 2-d velocity embedding\n"
        "    .var['fit_likelihood']      : per-gene kinetic fit quality"
    )
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    return adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="RNA velocity with scVelo for RetinoblastomaAtlas",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--n-top-genes",     type=int, default=2000)
    p.add_argument("--min-shared-counts", type=int, default=30)
    p.add_argument("--subset-cell-type", type=str, default="Cone_precursor")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_velocity(
        n_top_genes=args.n_top_genes,
        min_shared_counts=args.min_shared_counts,
        subset_cell_type=args.subset_cell_type,
    )
