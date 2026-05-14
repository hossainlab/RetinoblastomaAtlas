"""
=============================================================================
Script 09 — Cell–Cell Communication (Ligand–Receptor Analysis)
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
Retinoblastoma progression from intraocular to extraocular disease is not
a cell-autonomous process: tumour cells communicate with the tumour
microenvironment (TME) through secreted ligands and surface receptors.
Identifying which cell pairs exchange the most biologically relevant signals
— and how this changes between intraocular and extraocular tumours —
provides mechanistic hypotheses for invasion-driving pathways.

Specific questions addressed:
  1. Do tumour-associated macrophages (TAMs) in extraocular samples send
     more TGF-β, VEGF, or immunosuppressive (IL-10, ARG1) signals to
     tumour cells than in intraocular samples?
  2. Do tumour cone precursor cells signal back to macrophages via
     colony-stimulating factors (CSF1, IL-34) to promote M2 polarisation
     in the extraocular TME?
  3. Which receptor-ligand interactions differ most between the Subtype 1
     (non-invasive) and Subtype 2 (invasive) tumour populations?

METHOD — LIANA (Consensus L-R Scoring)
----------------------------------------
We use LIANA (Dimitrov et al. 2022, Nature Communications), a consensus
method that aggregates scores from five L-R databases and methods:
  - CellChat (Jin et al. 2021)
  - NicheNet (Browaeys et al. 2020)
  - CellPhoneDB v2 (Efremova et al. 2020)
  - NATMI (Hou et al. 2020)
  - Connectome (Raredon et al. 2022)

Using a consensus avoids the database-specific biases of single-method
approaches.  Interactions are ranked by a combined specificity-sensitivity
score (aggregate_rank; Dimitrov et al. 2022).

COMPARISON STRATEGY
--------------------
We run L-R analysis separately on:
  a. Intraocular samples (GSE168434 + GSE249995 intraocular subset)
  b. Extraocular samples (GSE249995 extraocular subset)
And compare interaction strengths across conditions.

REFERENCES
----------
- Dimitrov D et al. Comparison of methods and resources for cell–cell
  communication inference from single-cell RNA-seq data.
  Nat Commun. 2022;13:3224.
  https://doi.org/10.1038/s41467-022-30755-0

- Jin S et al. Inference and analysis of cell-cell communication using
  CellChat. Nat Commun. 2021;12:1088.
  https://doi.org/10.1038/s41467-021-21246-9

- Browaeys R et al. NicheNet: modeling intercellular communication by
  linking ligands to target genes. Nat Methods. 2020;17:159-162.
  https://doi.org/10.1038/s41592-019-0667-5

- Wan W et al. Single-cell transcriptome landscape of intraocular and
  extraocular retinoblastoma. Ophthalmology. 2025.
  https://doi.org/10.1016/j.ophtha.2025.01.011

INPUT  : data/processed/09_cellrank.h5ad
OUTPUT :
  data/processed/10_communication.h5ad
  results/figures/liana_dotplot_intraocular.pdf
  results/figures/liana_dotplot_extraocular.pdf
  results/figures/liana_chord_intraocular.pdf
  results/figures/liana_chord_extraocular.pdf
  results/figures/liana_top_interactions_comparison.pdf
  results/tables/liana_interactions_intraocular.csv
  results/tables/liana_interactions_extraocular.csv

USAGE  : pixi run python script/09_ligand_receptor_communication.py
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
IN_H5AD  = ROOT / "data" / "processed" / "09_cellrank.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "10_communication.h5ad"
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
# Key interactions of interest
# ---------------------------------------------------------------------------

PATHWAYS_OF_INTEREST = [
    # TGF-β signaling (invasion-promoting; Wan et al. 2025)
    "TGFB1_TGFBR1", "TGFB2_TGFBR2", "TGFB1_TGFBR2",
    # Macrophage recruitment by tumour cells
    "CSF1_CSF1R", "IL34_CSF1R",
    # Immunosuppression
    "IL10_IL10RA", "IL10_IL10RB",
    # Invasion / EMT-promoting
    "FN1_ITGB1", "FN1_ITGA5", "MMP9_ITGAV",
    # VEGF / angiogenesis
    "VEGFA_FLT1", "VEGFA_KDR",
    # Optic nerve invasion signals
    "CXCL12_CXCR4", "CXCL12_CXCR7",
    # SOX4-related signaling
    "WNT5A_FZD1", "WNT5A_ROR2",
]

# Cell type pairs of interest for comparison
CELL_TYPE_PAIRS_OF_INTEREST = [
    ("Cone_precursor",    "Microglia_TAM"),
    ("Microglia_TAM",     "Cone_precursor"),
    ("Cone_precursor",    "Endothelial"),
    ("Microglia_TAM",     "Endothelial"),
    ("Cone_precursor",    "Muller_glia"),
    ("Muller_glia",       "Cone_precursor"),
]


# ---------------------------------------------------------------------------
# LIANA analysis wrapper
# ---------------------------------------------------------------------------

def run_liana_on_subset(
    adata: sc.AnnData,
    label: str,
    groupby: str = "cell_type_broad",
    min_cells: int = 20,
) -> pd.DataFrame | None:
    """Run LIANA consensus L-R analysis on a subset of cells.

    Parameters
    ----------
    adata    : AnnData subset (already filtered to condition)
    label    : Descriptive label for logging ('intraocular' or 'extraocular')
    groupby  : Column in adata.obs used as cell type labels
    min_cells: Minimum cells per cell type required for inclusion

    Returns
    -------
    DataFrame of L-R interactions with consensus scores, or None if LIANA
    is not installed or too few cell types are present.

    SIGNIFICANCE:
    LIANA runs five methods in one call and returns an aggregate_rank column
    (lower = more reliable interaction).  By running on each condition
    separately and then merging on (ligand, receptor, source, target), we
    can compute condition-specific log fold changes in interaction scores
    — a direct measure of rewired signaling during invasion.
    """
    try:
        import liana
        from liana.method import rank_aggregate
    except ImportError:
        log.warning(
            "  LIANA not installed — skipping L-R analysis.\n"
            "  Install with: pip install liana"
        )
        return None

    # Filter out cell types with < min_cells
    ct_counts = adata.obs[groupby].value_counts()
    valid_cts  = ct_counts[ct_counts >= min_cells].index.tolist()
    if len(valid_cts) < 2:
        log.warning(
            f"  {label}: fewer than 2 cell types with >= {min_cells} cells. "
            "Skipping."
        )
        return None

    adata_sub = adata[adata.obs[groupby].isin(valid_cts)].copy()
    log.info(
        f"  {label}: {adata_sub.n_obs:,} cells, "
        f"{len(valid_cts)} cell types: {valid_cts}"
    )

    # LIANA needs .X to contain log-normalized expression
    # Use scvi_normalized layer if available, otherwise log-normalized .X
    if "scvi_normalized" in adata_sub.layers:
        import numpy as np
        adata_sub.X = adata_sub.layers["scvi_normalized"]
        # LIANA works on log1p-normalized data
        import scipy.sparse as sp
        if sp.issparse(adata_sub.X):
            adata_sub.X = adata_sub.X.toarray()
        adata_sub.X = np.log1p(adata_sub.X)

    try:
        rank_aggregate(
            adata_sub,
            groupby=groupby,
            expr_prop=0.1,        # min fraction of cells expressing gene
            min_cells=min_cells,
            use_raw=False,
            verbose=True,
            resource_name="consensus",  # uses CellChatDB + CellPhoneDB + NATMI
        )
        results = adata_sub.uns.get("liana_res", None)
        if results is None:
            log.warning(f"  {label}: No LIANA results found in adata.uns")
            return None
        results["condition"] = label
        log.info(
            f"  {label}: {len(results):,} L-R interactions computed\n"
            f"  Top 5 by aggregate_rank:\n"
            f"{results.sort_values('aggregate_rank').head(5)[['source', 'target', 'ligand', 'receptor', 'aggregate_rank']].to_string()}"
        )
        return results
    except Exception as e:
        log.error(f"  LIANA failed for {label}: {e}")
        return None


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_liana_dotplot(
    results: pd.DataFrame,
    label: str,
    out_path: Path,
    n_top: int = 30,
) -> None:
    """Dot plot of top L-R interactions (source → target).

    SIGNIFICANCE:
    Dot size encodes interaction specificity; dot colour encodes expression
    magnitude.  This plot allows rapid identification of the most specific
    and highly expressed interactions in each condition.
    """
    try:
        import liana
        from liana.pl import dotplot
    except ImportError:
        return

    top = results.sort_values("aggregate_rank").head(n_top)

    fig, ax = plt.subplots(figsize=(14, max(8, n_top * 0.3)))
    try:
        dotplot(
            liana_res=top,
            colour="specificity_rank",
            size="magnitude_rank",
            inverse_size=True,
            inverse_colour=True,
            ax=ax,
            show=False,
        )
    except Exception:
        # Fallback: simple scatter plot if LIANA dotplot API changes
        ax.scatter(
            top["ligand"] + "_" + top["receptor"],
            top["source"] + "→" + top["target"],
            s=50, c=top.get("aggregate_rank", 0), cmap="viridis_r",
        )

    ax.set_title(
        f"Top {n_top} L-R interactions — {label}\n"
        "(lower aggregate_rank = more reliable)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_interaction_comparison(
    intra_results: pd.DataFrame,
    extra_results: pd.DataFrame,
    out_path: Path,
    pathways: list[str] | None = None,
) -> None:
    """Bar chart comparing interaction strengths between intraocular and
    extraocular conditions for pathways of interest.

    SIGNIFICANCE:
    Directly quantifies how the TME signaling landscape changes during
    local extension.  Pathways enriched in extraocular samples are candidate
    therapeutic targets for preventing optic nerve invasion.
    """
    if intra_results is None or extra_results is None:
        return

    pathways = pathways or PATHWAYS_OF_INTEREST[:12]

    def score_interaction(results: pd.DataFrame, pair: str) -> float:
        """Mean aggregate_rank for a given ligand_receptor pair (lower = stronger)."""
        lig, rec = pair.split("_", 1)
        mask = (results["ligand"] == lig) & (results["receptor"] == rec)
        sub  = results[mask]
        if sub.empty:
            return np.nan
        return 1 - sub["aggregate_rank"].min()  # invert so higher = stronger

    intra_scores = {p: score_interaction(intra_results, p) for p in pathways}
    extra_scores = {p: score_interaction(extra_results, p) for p in pathways}

    df = pd.DataFrame({
        "intraocular":  intra_scores,
        "extraocular":  extra_scores,
    }).dropna(how="all")

    if df.empty:
        log.warning("  No matching pathway scores for comparison plot")
        return

    x   = np.arange(len(df))
    w   = 0.35
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(x - w/2, df["intraocular"].fillna(0), width=w,
           label="Intraocular",  color="#4C9BE8", edgecolor="none")
    ax.bar(x + w/2, df["extraocular"].fillna(0), width=w,
           label="Extraocular",  color="#E84C4C", edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Interaction score (1 − aggregate_rank)", fontsize=10)
    ax.set_title(
        "L-R interaction strengths: intraocular vs. extraocular\n"
        "(invasion-associated pathways)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_communication() -> sc.AnnData:
    """Cell–cell communication inference pipeline.

    Steps
    -----
    1. Load CellRank-annotated atlas.
    2. Split into intraocular and extraocular subsets.
    3. Run LIANA consensus L-R analysis on each subset.
    4. Identify top interactions in each condition.
    5. Compare invasion-relevant pathway interactions.
    6. Identify condition-specific interactions (unique to extraocular).
    7. Visualize and save.
    """

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading CellRank-annotated atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    log.info(f"  Datasets: {adata.obs['dataset'].value_counts().to_dict()}")

    # ---- 2. Split by disease stage / dataset ----------------------------
    log.info("\nSTEP 2 — Splitting into intraocular and extraocular subsets")
    if "disease_stage" in adata.obs.columns:
        stage_col = "disease_stage"
        intra_mask = adata.obs[stage_col].str.contains("intraocular", case=False, na=False)
        extra_mask = adata.obs[stage_col].str.contains("extraocular", case=False, na=False)
    else:
        log.warning(
            "  'disease_stage' column not found — using 'dataset' as proxy.\n"
            "  GSE168434 = intraocular, GSE249995 = mixed.\n"
            "  For accurate comparison, add 'disease_stage' annotation."
        )
        intra_mask = adata.obs["dataset"] == "GSE168434"
        extra_mask = (adata.obs["dataset"] == "GSE249995") & \
                     adata.obs.get("sample_id", pd.Series("", index=adata.obs_names)).str.contains(
                         "ex|extraoc", case=False, na=False
                     )

    adata_intra = adata[intra_mask].copy()
    adata_extra = adata[extra_mask].copy()
    log.info(
        f"  Intraocular: {adata_intra.n_obs:,} cells\n"
        f"  Extraocular: {adata_extra.n_obs:,} cells"
    )

    # ---- 3. LIANA on each subset ----------------------------------------
    log.info("\nSTEP 3 — Running LIANA consensus L-R analysis")
    log.info(
        "  LIANA aggregates five L-R methods:\n"
        "    CellChat, NicheNet, CellPhoneDB v2, NATMI, Connectome\n"
        "  aggregate_rank: lower = more reliable interaction (consensus)\n"
        "  expr_prop=0.10: gene must be expressed in ≥10% of source/target cells"
    )
    intra_results = run_liana_on_subset(adata_intra, "intraocular")
    extra_results = run_liana_on_subset(adata_extra, "extraocular")

    # ---- 4. Save interaction tables ------------------------------------
    log.info("\nSTEP 4 — Saving interaction tables")
    if intra_results is not None:
        intra_results.to_csv(TAB_DIR / "liana_interactions_intraocular.csv", index=False)
        log.info(f"  Saved intraocular interactions → {TAB_DIR / 'liana_interactions_intraocular.csv'}")
    if extra_results is not None:
        extra_results.to_csv(TAB_DIR / "liana_interactions_extraocular.csv", index=False)
        log.info(f"  Saved extraocular interactions → {TAB_DIR / 'liana_interactions_extraocular.csv'}")

    # ---- 5. Identify extraocular-enriched interactions ------------------
    log.info("\nSTEP 5 — Identifying extraocular-enriched interactions")
    if intra_results is not None and extra_results is not None:
        merge_cols = ["source", "target", "ligand", "receptor"]
        merged = intra_results[merge_cols + ["aggregate_rank"]].rename(
            columns={"aggregate_rank": "rank_intra"}
        ).merge(
            extra_results[merge_cols + ["aggregate_rank"]].rename(
                columns={"aggregate_rank": "rank_extra"}
            ),
            on=merge_cols,
            how="outer",
        )
        merged["rank_intra"] = merged["rank_intra"].fillna(1.0)
        merged["rank_extra"] = merged["rank_extra"].fillna(1.0)
        merged["delta_rank"] = merged["rank_intra"] - merged["rank_extra"]
        # Positive delta = lower rank in extra = stronger in extraocular
        merged_sorted = merged.sort_values("delta_rank", ascending=False)
        merged_sorted.to_csv(
            TAB_DIR / "liana_interactions_comparison.csv", index=False
        )
        log.info(
            f"  Top 5 interactions enriched in extraocular:\n"
            f"{merged_sorted.head(5)[merge_cols + ['delta_rank']].to_string()}"
        )

    # ---- 6. Visualize ---------------------------------------------------
    log.info("\nSTEP 6 — Generating communication visualizations")
    if intra_results is not None:
        plot_liana_dotplot(
            intra_results, "Intraocular",
            FIG_DIR / "liana_dotplot_intraocular.pdf",
        )
    if extra_results is not None:
        plot_liana_dotplot(
            extra_results, "Extraocular",
            FIG_DIR / "liana_dotplot_extraocular.pdf",
        )
    plot_interaction_comparison(
        intra_results, extra_results,
        FIG_DIR / "liana_top_interactions_comparison.pdf",
    )

    # ---- 7. Store results in adata.uns ---------------------------------
    log.info("\nSTEP 7 — Storing results in adata.uns")
    if intra_results is not None:
        adata.uns["liana_intraocular"] = intra_results
    if extra_results is not None:
        adata.uns["liana_extraocular"] = extra_results

    # ---- 8. Save --------------------------------------------------------
    log.info(f"\nSTEP 8 — Saving communication-annotated atlas → {OUT_H5AD.name}")
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    return adata


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_communication()
