"""
=============================================================================
Script 08 — Cell Fate Mapping with CellRank
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
RNA velocity (script 07) provides a LOCAL measure of transcriptional
directionality — each cell has a vector pointing toward its next state.
CellRank (Lange et al. 2022, Nature Methods) extends this to a GLOBAL
probabilistic framework for cell fate prediction:

  1. It combines RNA velocity vectors with a cell–cell similarity kernel
     (computed in the scVI latent space) into a Markov chain transition
     matrix.  This produces robust fate predictions even for cells with
     noisy individual velocity vectors.

  2. It identifies INITIAL and TERMINAL STATES of the trajectory
     automatically using a directed random walk framework (GPCCA — General
     Perron Cluster Cluster Analysis).

  3. It computes FATE PROBABILITIES: for each cell, the probability of
     reaching each identified terminal state.  In the context of this
     project, these represent the probability of a cone precursor cell
     committing to:
       - Terminal state A: mature, non-invasive cone-like tumour (Subtype 1)
       - Terminal state B: dedifferentiated, invasive tumour (Subtype 2 /
         extraocular extension)

  4. It identifies DRIVER GENES: genes whose expression correlates with the
     probability of reaching each terminal state — candidates for
     microenvironmental regulators of the Subtype 1 → 2 transition.

SIGNIFICANCE IN THIS PROJECT
------------------------------
The CellRank fate mapping directly addresses Research Aim 3:
  "Reconstruct the cone precursor → invasive tumour cell transition
   trajectory and identify transcription factor drivers of dedifferentiation
   using RNA velocity and CellRank."

By comparing fate probabilities between intraocular (GSE168434) and
extraocular (GSE249995) samples, we can quantify whether extraocular tumour
cells are preferentially committed to the invasive terminal state —
providing mechanistic insight into why extraocular tumours are more
aggressive.

KERNEL COMBINATION STRATEGY
------------------------------
CellRank supports multiple kernels; we use a combined kernel:
  - VelocityKernel (RNA velocity direction, weight=0.8)
  - ConnectivityKernel (scVI-based KNN similarity, weight=0.2)
The ConnectivityKernel provides stability for cells with low velocity
confidence, preventing spurious fate assignments.

REFERENCES
----------
- Lange M et al. CellRank for directed single-cell fate mapping.
  Nat Methods. 2022;19:159-170.
  https://doi.org/10.1038/s41592-021-01346-6

- Reuter B et al. Generalized Perron Cluster Cluster Analysis (GPCCA):
  Identifying macrostates in non-reversible Markov chains.
  J Chem Theory Comput. 2018;14(7):3499-3510.
  https://doi.org/10.1021/acs.jctc.8b00079

- Bergen V et al. Generalizing RNA velocity to transient cell states through
  dynamical modeling. Nat Biotechnol. 2020;38(12):1408-1414.
  https://doi.org/10.1038/s41587-020-0591-3

INPUT  : data/processed/08_velocity.h5ad
OUTPUT :
  data/processed/09_cellrank.h5ad
  results/figures/cellrank_macrostates_umap.pdf
  results/figures/cellrank_fate_probabilities_umap.pdf
  results/figures/cellrank_driver_genes_heatmap.pdf
  results/figures/cellrank_lineage_drivers.pdf
  results/tables/cellrank_driver_genes.csv
  results/tables/cellrank_fate_probabilities.csv

USAGE  : pixi run python script/08_cellrank_fate_mapping.py [--n-states 4]
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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data" / "processed" / "08_velocity.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "09_cellrank.h5ad"
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
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_macrostates(adata: sc.AnnData, out_path: Path) -> None:
    """UMAP of CellRank macrostates (initial, terminal, transient)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, key in zip(axes, ["macrostates_fwd", "initial_states"]):
        if key in adata.obs.columns:
            sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                       size=3, alpha=0.8, title=key, legend_loc="right margin")
        else:
            ax.set_visible(False)
    plt.suptitle("CellRank macrostates", fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_fate_probabilities(adata: sc.AnnData, out_path: Path) -> None:
    """UMAP coloured by fate probability toward each terminal state."""
    fate_cols = [c for c in adata.obs.columns if c.startswith("fate_prob_")]
    if not fate_cols:
        log.warning("  No fate probability columns found — skipping fate prob plot")
        return

    ncols = min(len(fate_cols), 4)
    nrows = int(np.ceil(len(fate_cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes = np.array(axes).flatten() if nrows * ncols > 1 else [axes]

    for ax, key in zip(axes, fate_cols):
        sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                   size=2, cmap="viridis", vmin=0, vmax=1,
                   title=key.replace("fate_prob_", "→ "))
    for ax in axes[len(fate_cols):]:
        ax.set_visible(False)

    plt.suptitle(
        "Fate probabilities (CellRank GPCCA)\n"
        "Each panel shows P(cell reaches terminal state X)",
        fontsize=11, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_driver_gene_heatmap(
    driver_dict: dict[str, pd.DataFrame],
    adata: sc.AnnData,
    out_path: Path,
    n_genes: int = 20,
) -> None:
    """Heatmap of top driver genes per terminal state."""
    if not driver_dict:
        return
    all_drivers = pd.concat(driver_dict.values())
    top_genes   = all_drivers.index[:n_genes].tolist()
    present     = [g for g in top_genes if g in adata.var_names]
    if not present:
        return

    sc.pl.matrixplot(
        adata,
        var_names=present,
        groupby="cell_type_broad",
        standard_scale="var",
        cmap="RdBu_r",
        show=False,
        save=False,
    )
    plt.suptitle(
        f"Top {n_genes} CellRank driver genes\n(z-scored, grouped by cell type)",
        fontsize=11, fontweight="bold",
    )
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_cellrank(n_states: int = 4) -> sc.AnnData:
    """CellRank fate mapping pipeline.

    Steps
    -----
    1. Load velocity-annotated atlas.
    2. Build combined VelocityKernel + ConnectivityKernel.
    3. Compute cell–cell transition matrix.
    4. Run GPCCA to identify macrostates (initial, transient, terminal).
    5. Compute fate probabilities for each terminal state.
    6. Identify lineage driver genes.
    7. Visualize and save.
    """
    try:
        import cellrank as cr
        from cellrank.kernels import VelocityKernel, ConnectivityKernel
        from cellrank.estimators import GPCCA
    except ImportError:
        raise ImportError(
            "CellRank not installed. Install with:\n"
            "  pixi add cellrank\n"
            "or: pip install cellrank"
        )

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading velocity-annotated atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 2. Build combined kernel ---------------------------------------
    log.info("\nSTEP 2 — Building combined velocity + connectivity kernel")
    log.info(
        "  VelocityKernel:     uses RNA velocity cosine similarities\n"
        "  ConnectivityKernel: uses scVI KNN graph edge weights\n"
        "  Combined (0.8 velocity + 0.2 connectivity) for robustness\n"
        "  WHY combine: Pure velocity kernel is noisy; adding connectivity\n"
        "  stabilises transitions for low-confidence velocity cells."
    )
    vk = VelocityKernel(adata)
    vk.compute_transition_matrix()

    ck = ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    combined_kernel = 0.8 * vk + 0.2 * ck
    combined_kernel.write_to_adata()
    log.info("  Combined kernel computed and written to adata.obsp")

    # ---- 3. GPCCA estimator ---------------------------------------------
    log.info(f"\nSTEP 3 — Running GPCCA (n_states={n_states})")
    log.info(
        "  GPCCA partitions the transition matrix into macrostates\n"
        "  using Schur decomposition of the Markov chain.\n"
        "  Terminal states = absorbing states (low outflow probability)\n"
        "  Initial states  = high-outflow source states"
    )
    g = GPCCA(combined_kernel)
    g.compute_schur(n_components=n_states + 1)
    g.compute_macrostates(n_states=n_states, cluster_key="cell_type_broad")

    # Log identified states
    if hasattr(g, "macrostates"):
        log.info(f"  Macrostates identified: {g.macrostates.cat.categories.tolist()}")

    g.set_terminal_states_from_macrostates()
    log.info(f"  Terminal states set: {g.terminal_states.cat.categories.tolist()}")

    # ---- 4. Fate probabilities ------------------------------------------
    log.info("\nSTEP 4 — Computing fate probabilities")
    log.info(
        "  Fate probability P(cell → terminal state X) is computed via\n"
        "  absorption probability in the Markov chain.  Each cell gets a\n"
        "  probability vector over all terminal states summing to 1."
    )
    g.compute_fate_probabilities()

    # Store fate probabilities in adata.obs
    if hasattr(g, "fate_probabilities"):
        for state in g.terminal_states.cat.categories:
            col_name = f"fate_prob_{state.replace(' ', '_')}"
            adata.obs[col_name] = g.fate_probabilities[state].values
        log.info(f"  Fate probability columns: {[c for c in adata.obs.columns if c.startswith('fate_prob_')]}")

    # ---- 5. Macrostate annotation back to adata -------------------------
    adata.obs["macrostates_fwd"] = g.macrostates
    if hasattr(g, "initial_states"):
        adata.obs["initial_states"] = g.initial_states

    # ---- 6. Driver genes ------------------------------------------------
    log.info("\nSTEP 5 — Identifying lineage driver genes")
    log.info(
        "  Driver genes are those whose expression is most correlated\n"
        "  with the fate probability of each terminal state.\n"
        "  Positive correlation = induction in cells committing to that fate\n"
        "  Negative correlation = repression during commitment\n"
        "  These are candidate regulators of the invasion transition."
    )
    driver_dict = {}
    for state in g.terminal_states.cat.categories:
        col_name = f"fate_prob_{state.replace(' ', '_')}"
        if col_name not in adata.obs.columns:
            continue
        try:
            drivers = g.compute_lineage_drivers(
                lineages=state,
                cluster_key="cell_type_broad",
                use_raw=False,
                return_drivers=True,
            )
            if drivers is not None:
                driver_dict[state] = drivers
                log.info(
                    f"  Top 5 drivers → {state}: "
                    f"{drivers.head(5).index.tolist()}"
                )
        except Exception as e:
            log.warning(f"  Driver gene computation for {state} failed: {e}")

    # Save driver gene tables
    if driver_dict:
        all_drivers = pd.concat(
            [df.assign(terminal_state=s) for s, df in driver_dict.items()]
        )
        all_drivers.to_csv(TAB_DIR / "cellrank_driver_genes.csv")
        log.info(f"  Saved driver genes → {TAB_DIR / 'cellrank_driver_genes.csv'}")

    # Save fate probabilities table
    fate_cols = [c for c in adata.obs.columns if c.startswith("fate_prob_")]
    meta_cols = ["sample_id", "dataset", "cell_type_broad", "rb_subtype",
                 "velocity_pseudotime"]
    meta_cols = [c for c in meta_cols if c in adata.obs.columns]
    adata.obs[meta_cols + fate_cols].to_csv(
        TAB_DIR / "cellrank_fate_probabilities.csv"
    )
    log.info(f"  Saved fate probabilities → {TAB_DIR / 'cellrank_fate_probabilities.csv'}")

    # ---- 7. Plots -------------------------------------------------------
    log.info("\nSTEP 6 — Generating CellRank visualizations")
    plot_macrostates(adata, FIG_DIR / "cellrank_macrostates_umap.pdf")
    plot_fate_probabilities(adata, FIG_DIR / "cellrank_fate_probabilities_umap.pdf")
    plot_driver_gene_heatmap(
        driver_dict, adata,
        FIG_DIR / "cellrank_driver_genes_heatmap.pdf",
    )

    # Fate probability comparison: intraocular vs. extraocular
    if "disease_stage" in adata.obs.columns and fate_cols:
        log.info("  Plotting fate probability distribution by disease stage...")
        fig, axes = plt.subplots(1, len(fate_cols), figsize=(6 * len(fate_cols), 5))
        if len(fate_cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, fate_cols):
            sc.pl.violin(
                adata, keys=col, groupby="disease_stage",
                ax=ax, show=False, rotation=30,
                title=col.replace("fate_prob_", "P(→ "),
            )
        plt.suptitle(
            "Fate probabilities: intraocular vs. extraocular",
            fontsize=11, fontweight="bold",
        )
        plt.tight_layout()
        fig.savefig(
            FIG_DIR / "cellrank_fate_probs_by_disease_stage.pdf",
            dpi=200, bbox_inches="tight",
        )
        plt.close(fig)
        log.info("  Saved: cellrank_fate_probs_by_disease_stage.pdf")

    # ---- 8. Save --------------------------------------------------------
    log.info(f"\nSTEP 7 — Saving CellRank-annotated atlas → {OUT_H5AD.name}")
    log.info(
        "  Contents added:\n"
        "    .obs['macrostates_fwd']     : GPCCA macrostate labels\n"
        "    .obs['fate_prob_*']         : per-terminal-state fate probs\n"
        "    .obsp['T_fwd']              : forward transition matrix\n"
        "    .uns['cellrank_terminal_states']: terminal state names"
    )
    adata.uns["cellrank_terminal_states"] = \
        g.terminal_states.cat.categories.tolist()
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    return adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="CellRank fate mapping for RetinoblastomaAtlas",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--n-states", type=int, default=4,
                   help="Number of macrostates for GPCCA")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_cellrank(n_states=args.n_states)
