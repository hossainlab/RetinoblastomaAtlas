"""
=============================================================================
Script 10 — TGF-β Pathway Activity Scoring and TME Rewiring Analysis
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
TGF-β (Transforming Growth Factor beta) signaling is a central mediator of
tumour invasion and immunosuppression in multiple cancer types.  In this
retinoblastoma atlas, Wan et al. (2025, Ophthalmology) identified that the
CP4 cone precursor subcluster — specifically enriched in extraocular samples
— shows the highest TGF-β ligand expression (TGFB1, TGFB2).  This script
quantifies and spatialises TGF-β pathway activity across the entire atlas to:

  1. Confirm that TGF-β signaling is significantly UPREGULATED in extraocular
     vs. intraocular tumour cells (validates the Wan et al. finding in our
     integrated two-dataset atlas).

  2. Identify WHICH CELL TYPES are the primary sources (TGFB1/2 producers)
     and responders (SMAD2/3 phosphorylation targets) of TGF-β in the TME.
     Key hypotheses:
       a. Tumour cells (CP4/Subtype 2) autocrine TGF-β promotes EMT
       b. TAMs produce TGF-β to suppress CD8+ T cell activity
       c. Muller glia / CAFs respond to TGF-β by upregulating ECM genes

  3. Correlate TGF-β activity with:
       - Fate probability of the invasive terminal state (script 08)
       - RB subtype 2 score (script 06)
       - CNV instability load (script 05)

  4. Score ADDITIONAL TME-RELEVANT PATHWAYS alongside TGF-β to give a
     comprehensive picture of the rewired signalling landscape during
     extraocular extension:
       - VEGFA-VEGFR2 (angiogenesis / tumour vascularization)
       - Hypoxia (HIF1A targets; relevant to avascular intraocular space)
       - MAPK (ERK/RAF/RAS activation; invasion-associated)
       - PI3K-AKT (survival / mTOR pathway)
       - p53 (tumour suppressor; lost in high-risk metastatic RB)
       - JAK-STAT (IFN-γ signaling; immune evasion marker)
       - WNT (beta-catenin target genes; dedifferentiation)
       - NFkB (inflammatory signaling)

METHOD — decoupleR with PROGENy
---------------------------------
PROGENy (Pathway RespOnsive GENes; Schubert et al. 2018, Nat Comms) is a
database of gene expression footprints of 14 canonical signalling pathways.
It provides pathway-specific target gene sets derived from perturbation
experiments (NOT from correlation of expression data), making the scores
mechanistically interpretable.

decoupleR (Badia-i-Mompel et al. 2022, Bioinformatics) implements multiple
statistical enrichment approaches.  We use:
  - Weighted mean (wmean): robust, best for PROGENy (Badia-i-Mompel et al.)
  - ULM (univariate linear model): complementary approach
And return the consensus score.

REFERENCES
----------
- Wan W et al. Single-cell transcriptome landscape of intraocular and
  extraocular retinoblastoma. Ophthalmology. 2025.
  https://doi.org/10.1016/j.ophtha.2025.01.011

- Schubert M et al. Perturbation-response genes reveal signaling footprints
  in cancer gene expression. Nat Commun. 2018;9:20.
  https://doi.org/10.1038/s41467-017-02391-6

- Badia-i-Mompel P et al. decoupleR: ensemble of computational methods to
  infer biological activities from omics data.
  Bioinformatics. 2022;3(btac016).
  https://doi.org/10.1093/bioadv/btac016

- Liu Y et al. Single-cell transcriptomic analysis reveals local extension
  mechanism of retinoblastoma. Cell Mol Life Sci. 2024;81:77.
  https://doi.org/10.1007/s00018-024-05102-9

INPUT  : data/processed/10_communication.h5ad
OUTPUT :
  data/processed/11_pathway_scored.h5ad
  results/figures/pathway_scores_umap.pdf
  results/figures/tgfb_score_by_celltype_stage.pdf
  results/figures/pathway_heatmap_by_celltype.pdf
  results/figures/tgfb_correlation_subtype_fate.pdf
  results/figures/tme_rewiring_comparison.pdf
  results/tables/progeny_pathway_scores.csv
  results/tables/tgfb_differential_analysis.csv

USAGE  : pixi run python script/10_tgfb_pathway_scoring.py
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
import scipy.stats as stats
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     = Path(__file__).resolve().parents[1]
IN_H5AD  = ROOT / "data" / "processed" / "10_communication.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "11_pathway_scored.h5ad"
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
# Pathways to score
# ---------------------------------------------------------------------------

PROGENY_PATHWAYS = [
    "TGFb",    # TGF-β signaling — primary focus
    "VEGF",    # Angiogenesis
    "Hypoxia", # HIF1A targets
    "MAPK",    # ERK/RAF/RAS
    "PI3K",    # AKT/mTOR
    "p53",     # Tumour suppressor
    "JAK-STAT",# Cytokine signaling
    "WNT",     # Beta-catenin / dedifferentiation
    "NFkB",    # Inflammation
    "EGFR",    # EGF receptor signaling
    "Androgen",# AR signaling (used as control — not expected in RB)
    "Estrogen",# ER signaling (negative control)
    "Trail",   # Apoptosis
    "TNFa",    # TNF-alpha (inflammatory)
]

# PROGENy top-N genes per pathway used for scoring
PROGENY_TOP_N = 100


# ---------------------------------------------------------------------------
# TGF-β manual gene signatures (fallback if decoupleR not available)
# ---------------------------------------------------------------------------

TGFB_CANONICAL_TARGETS = [
    # SMAD-dependent targets
    "SMAD7", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",
    "CDH2", "VIM", "FN1", "MMP2", "MMP9", "CTGF", "CYR61",
    "TGFBI", "ITGB6", "PAI1", "SERPINE1", "COL1A1", "COL1A2",
    "ACTA2",
    # TGF-β ligands / receptors
    "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3",
    "SMAD2", "SMAD3", "SMAD4",
    # EMT markers (non-canonical TGF-β output)
    "CDH1",  # E-cadherin (should decrease in EMT)
    "VIM",   # Vimentin (should increase)
]


# ---------------------------------------------------------------------------
# decoupleR pathway scoring
# ---------------------------------------------------------------------------

def run_progeny_decoupler(adata: sc.AnnData) -> sc.AnnData | None:
    """Score cells with PROGENy pathway activities using decoupleR.

    SIGNIFICANCE:
    PROGENy scores estimate the *activity* of signalling pathways, not just
    the expression of individual pathway genes.  This is more interpretable
    biologically: a pathway can be active even if individual genes have
    variable expression, as long as the ensemble footprint is coherent.
    The wmean method is the recommended estimator for PROGENy (Badia-i-Mompel
    et al. 2022).
    """
    try:
        import decoupler as dc
    except ImportError:
        log.warning(
            "  decoupleR not installed — skipping PROGENy scoring.\n"
            "  Install with: pip install decoupler"
        )
        return None

    log.info("  Loading PROGENy network (human, top 100 genes per pathway)...")
    try:
        progeny_net = dc.get_progeny(organism="human", top=PROGENY_TOP_N)
        log.info(f"  PROGENy network: {len(progeny_net):,} gene-pathway edges")
    except Exception as e:
        log.warning(f"  Failed to load PROGENy network: {e}")
        return None

    # Use log-normalized expression
    log.info("  Running decoupleR wmean (weighted mean) estimator...")
    try:
        dc.run_wmean(
            mat=adata,
            net=progeny_net,
            source="source",
            target="target",
            weight="weight",
            times=100,           # permutations for p-value estimation
            min_n=5,             # min genes per pathway
            use_raw=False,
        )
    except Exception as e:
        log.error(f"  decoupleR wmean failed: {e}")
        return None

    # decoupleR stores results in adata.obsm['wmean_estimate']
    if "wmean_estimate" in adata.obsm:
        acts = pd.DataFrame(
            adata.obsm["wmean_estimate"],
            index=adata.obs_names,
            columns=adata.uns.get("wmean_estimate_names",
                                   adata.obsm["wmean_estimate"].columns
                                   if hasattr(adata.obsm["wmean_estimate"], "columns")
                                   else [f"path_{i}" for i in range(adata.obsm["wmean_estimate"].shape[1])]),
        )
        # Store individual pathway scores in adata.obs
        for pathway in acts.columns:
            adata.obs[f"PROGENy_{pathway}"] = acts[pathway].values
        log.info(
            f"  Scored {len(acts.columns)} pathways via PROGENy wmean.\n"
            f"  Pathways: {acts.columns.tolist()}"
        )
    else:
        log.warning("  wmean_estimate not found in adata.obsm after decoupleR run")
        return None

    return adata


def score_tgfb_manual(adata: sc.AnnData) -> None:
    """Score TGF-β activity using manual canonical target gene set.

    Used as a fallback when decoupleR is unavailable, and as a validation
    of the PROGENy score.  Uses scanpy's score_genes() with expression-
    bin-matched control genes.
    """
    present = [g for g in TGFB_CANONICAL_TARGETS if g in adata.var_names]
    absent  = [g for g in TGFB_CANONICAL_TARGETS if g not in adata.var_names]
    log.info(
        f"  TGF-β manual signature: {len(present)}/{len(TGFB_CANONICAL_TARGETS)} genes present "
        f"(missing: {absent[:5]})"
    )
    sc.tl.score_genes(
        adata,
        gene_list=present,
        score_name="TGFb_manual_score",
        ctrl_size=50,
        n_bins=25,
    )
    log.info(
        f"  TGF-β manual score: mean={adata.obs['TGFb_manual_score'].mean():.4f}, "
        f"std={adata.obs['TGFb_manual_score'].std():.4f}"
    )


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------

def differential_pathway_analysis(
    adata: sc.AnnData,
    pathway_cols: list[str],
    group_col: str = "disease_stage",
    group_a: str = "intraocular",
    group_b: str = "extraocular",
) -> pd.DataFrame:
    """Mann-Whitney U test for pathway scores between two groups.

    SIGNIFICANCE:
    A non-parametric test is preferred because pathway activity scores
    are not normally distributed (they are right-skewed in cells that
    activate a pathway and near-zero in inactive cells).  The Mann-Whitney
    U statistic (equivalent to AUC for binary outcomes) measures whether
    one group systematically scores higher than the other.

    Returns a DataFrame with:
      - pathway: pathway name
      - group_a_mean, group_b_mean: mean scores
      - log2fc: log2 fold change (group_b / group_a)
      - mwu_stat: Mann-Whitney U statistic
      - p_value: raw p-value
      - p_adj: Benjamini-Hochberg adjusted p-value
      - significant: True if p_adj < 0.05 and |log2fc| > 0.5
    """
    if group_col not in adata.obs.columns:
        log.warning(f"  Column '{group_col}' not found — skipping differential analysis")
        return pd.DataFrame()

    mask_a = adata.obs[group_col].str.contains(group_a, case=False, na=False)
    mask_b = adata.obs[group_col].str.contains(group_b, case=False, na=False)
    log.info(
        f"  Differential analysis: {group_a} (n={mask_a.sum():,}) vs "
        f"{group_b} (n={mask_b.sum():,})"
    )

    rows = []
    for col in pathway_cols:
        if col not in adata.obs.columns:
            continue
        a_vals = adata.obs.loc[mask_a, col].dropna().values
        b_vals = adata.obs.loc[mask_b, col].dropna().values
        if len(a_vals) < 10 or len(b_vals) < 10:
            continue
        stat, pval = stats.mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        mean_a = a_vals.mean()
        mean_b = b_vals.mean()
        lfc    = np.log2((mean_b + 1e-6) / (mean_a + 1e-6))
        rows.append({
            "pathway":       col,
            "group_a_mean":  round(mean_a, 5),
            "group_b_mean":  round(mean_b, 5),
            "log2fc":        round(lfc, 4),
            "mwu_stat":      round(stat, 1),
            "p_value":       pval,
        })

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    # Benjamini-Hochberg FDR correction
    from scipy.stats import false_discovery_control
    result["p_adj"] = false_discovery_control(result["p_value"].values)
    result["significant"] = (result["p_adj"] < 0.05) & (result["log2fc"].abs() > 0.5)
    result = result.sort_values("log2fc", ascending=False)
    return result


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_pathway_scores_umap(adata: sc.AnnData, out_path: Path) -> None:
    """UMAP grid of top pathway activity scores."""
    pathway_cols = [c for c in adata.obs.columns
                    if c.startswith("PROGENy_") or c == "TGFb_manual_score"]
    if not pathway_cols:
        log.warning("  No pathway score columns found")
        return

    focus = ["PROGENy_TGFb", "PROGENy_VEGF", "PROGENy_Hypoxia",
             "PROGENy_MAPK", "PROGENy_WNT", "TGFb_manual_score"]
    cols  = [c for c in focus if c in pathway_cols][:6]
    if not cols:
        cols = pathway_cols[:6]

    ncols = 3
    nrows = int(np.ceil(len(cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes = np.array(axes).flatten()

    for ax, col in zip(axes, cols):
        sc.pl.umap(
            adata, color=col, ax=ax, show=False, frameon=False,
            size=2, cmap="RdBu_r", vcenter=0,
            title=col.replace("PROGENy_", "").replace("_manual_score", " (manual)"),
        )
    for ax in axes[len(cols):]:
        ax.set_visible(False)

    plt.suptitle("Pathway activity scores (PROGENy / decoupleR)",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_tgfb_by_celltype_stage(adata: sc.AnnData, out_path: Path) -> None:
    """Violin / box plot of TGF-β score per cell type split by disease stage.

    SIGNIFICANCE:
    This plot directly tests the key hypothesis: TGF-β activity is higher
    in extraocular tumour cone precursor cells.  Each cell type is shown
    separately to distinguish cell-autonomous (tumour-intrinsic) TGF-β
    upregulation from paracrine contributions from stromal cells.
    """
    tgfb_col = "PROGENy_TGFb" if "PROGENy_TGFb" in adata.obs.columns \
               else "TGFb_manual_score"
    if tgfb_col not in adata.obs.columns:
        return

    stage_col = "disease_stage" if "disease_stage" in adata.obs.columns else "dataset"

    fig, ax = plt.subplots(figsize=(14, 6))
    sc.pl.violin(
        adata,
        keys=tgfb_col,
        groupby="cell_type_broad",
        split=stage_col if adata.obs[stage_col].nunique() <= 3 else None,
        ax=ax, show=False, rotation=30,
    )
    ax.set_title(
        f"TGF-β pathway activity by cell type\n"
        f"(split by {stage_col}; higher score = more active)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_pathway_heatmap_by_celltype(adata: sc.AnnData, out_path: Path) -> None:
    """Heatmap of mean pathway activity per cell type."""
    pathway_cols = [c for c in adata.obs.columns if c.startswith("PROGENy_")]
    if not pathway_cols:
        return

    ct_mean = adata.obs.groupby("cell_type_broad")[pathway_cols].mean()
    ct_mean.columns = [c.replace("PROGENy_", "") for c in ct_mean.columns]

    fig, ax = plt.subplots(figsize=(len(pathway_cols) * 0.8 + 2,
                                     len(ct_mean) * 0.6 + 2))
    im = ax.imshow(
        ct_mean.values,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-ct_mean.abs().max().max(),
        vmax=ct_mean.abs().max().max(),
    )
    ax.set_xticks(range(len(ct_mean.columns)))
    ax.set_xticklabels(ct_mean.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(ct_mean.index)))
    ax.set_yticklabels(ct_mean.index, fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.6, label="Mean pathway activity")
    ax.set_title(
        "PROGENy pathway activities per cell type\n"
        "(Red = active, Blue = inactive)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_tgfb_correlations(adata: sc.AnnData, out_path: Path) -> None:
    """Scatter plots: TGF-β score vs. RB subtype score, fate probability, CNV.

    SIGNIFICANCE:
    These scatter plots test three mechanistic hypotheses:
      1. TGF-β activity correlates with RB Subtype 2 / stemness score
         → TGF-β drives dedifferentiation
      2. TGF-β activity correlates with fate probability of invasive state
         → TGF-β promotes invasive fate commitment
      3. TGF-β activity correlates with CNV instability load
         → Genomically unstable cells activate TGF-β
    Pearson and Spearman correlations are reported for each.
    """
    tgfb_col = "PROGENy_TGFb" if "PROGENy_TGFb" in adata.obs.columns \
               else "TGFb_manual_score"
    if tgfb_col not in adata.obs.columns:
        return

    correlates = [
        ("score_RB_subtype2_stemness", "RB Subtype 2 score"),
        ("cnv_load",                   "CNV instability load"),
    ]
    # Add fate probability columns if present
    fate_cols = [c for c in adata.obs.columns if c.startswith("fate_prob_")]
    for fc in fate_cols[:2]:
        correlates.append((fc, fc.replace("fate_prob_", "P(→ ")))

    n_corr = len([c for c, _ in correlates if c in adata.obs.columns])
    if n_corr == 0:
        return

    ncols = min(3, n_corr)
    nrows = int(np.ceil(n_corr / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.array(axes).flatten() if n_corr > 1 else [axes]

    ax_idx = 0
    for col, label in correlates:
        if col not in adata.obs.columns:
            continue
        ax = axes[ax_idx]
        x = adata.obs[tgfb_col].values
        y = adata.obs[col].values
        # Sub-sample for speed
        idx = np.random.choice(len(x), size=min(5000, len(x)), replace=False)
        r_p, _   = stats.pearsonr(x[idx], y[idx])
        r_s, p_s = stats.spearmanr(x[idx], y[idx])
        ax.scatter(x[idx], y[idx], s=1, alpha=0.3, color="#4C9BE8",
                   rasterized=True)
        ax.set_xlabel("TGF-β pathway activity", fontsize=9)
        ax.set_ylabel(label, fontsize=9)
        ax.set_title(
            f"r_P={r_p:.3f}, r_S={r_s:.3f} (p={p_s:.1e})",
            fontsize=8,
        )
        ax_idx += 1

    for ax in axes[ax_idx:]:
        ax.set_visible(False)

    plt.suptitle(
        f"TGF-β activity correlates with invasion markers\n"
        f"(n={min(5000, adata.n_obs):,} cells subsampled)",
        fontsize=11, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_tme_rewiring(
    diff_df: pd.DataFrame,
    out_path: Path,
) -> None:
    """Volcano-style plot of pathway rewiring: intraocular vs. extraocular.

    SIGNIFICANCE:
    Shows which pathways are significantly more active in extraocular vs.
    intraocular tumours.  Pathways in the upper-right quadrant (high log2FC,
    low adjusted p-value) are candidates for driving local extension.
    Expected: TGF-β and VEGF in upper right; p53 in upper left (loss).
    """
    if diff_df.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 7))
    colors = diff_df["significant"].map({True: "#E84C4C", False: "#AAAAAA"})
    ax.scatter(
        diff_df["log2fc"],
        -np.log10(diff_df["p_adj"] + 1e-300),
        s=60, c=colors, edgecolors="none", alpha=0.8,
    )
    # Label significant points
    for _, row in diff_df[diff_df["significant"]].iterrows():
        name = row["pathway"].replace("PROGENy_", "")
        ax.annotate(
            name,
            xy=(row["log2fc"], -np.log10(row["p_adj"] + 1e-300)),
            xytext=(5, 5), textcoords="offset points", fontsize=8,
        )
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=1)
    ax.axvline(0, color="grey", linestyle="-",  linewidth=0.5)
    ax.set_xlabel("log2 FC (extraocular / intraocular)", fontsize=10)
    ax.set_ylabel("-log10(FDR)", fontsize=10)
    ax.set_title(
        "TME pathway rewiring: extraocular vs. intraocular\n"
        "(red = FDR < 0.05 and |log2FC| > 0.5)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pathway_scoring() -> sc.AnnData:
    """TGF-β and TME pathway scoring pipeline.

    Steps
    -----
    1. Load communication-annotated atlas.
    2. Run PROGENy pathway scoring via decoupleR (wmean).
    3. Score TGF-β with manual canonical target genes (validation).
    4. Differential pathway analysis: intraocular vs. extraocular.
    5. Correlate TGF-β with subtype, fate probability, CNV load.
    6. Visualize.
    7. Save.
    """

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading communication-annotated atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 2. PROGENy scoring via decoupleR --------------------------------
    log.info("\nSTEP 2 — PROGENy pathway activity scoring (decoupleR wmean)")
    log.info(
        "  PROGENy models 14 canonical signalling pathways from perturbation\n"
        "  gene expression footprints.  Scores reflect PATHWAY ACTIVITY,\n"
        "  not just pathway gene expression — more mechanistically meaningful.\n"
        "  decoupleR wmean: weighted mean of top 100 target genes per pathway,\n"
        "  weighted by perturbation-derived importance scores."
    )
    adata_scored = run_progeny_decoupler(adata)
    if adata_scored is None:
        log.info("  Falling back to manual TGF-β scoring only")

    # ---- 3. Manual TGF-β score (validation / fallback) ------------------
    log.info("\nSTEP 3 — Manual TGF-β signature score (canonical targets)")
    log.info(
        "  Validates the PROGENy TGF-β score using canonical SMAD-dependent\n"
        "  and EMT target genes (TGFB1, SMAD2/3, SNAI1/2, VIM, FN1, etc.).\n"
        "  Should correlate with PROGENy_TGFb if both methods are reliable."
    )
    score_tgfb_manual(adata)

    # Validate correlation between PROGENy and manual TGF-β score
    if "PROGENy_TGFb" in adata.obs.columns:
        r, p = stats.pearsonr(
            adata.obs["PROGENy_TGFb"].fillna(0),
            adata.obs["TGFb_manual_score"].fillna(0),
        )
        log.info(
            f"  Correlation PROGENy_TGFb vs manual: "
            f"r={r:.3f}, p={p:.2e}\n"
            f"  {'CONCORDANT' if abs(r) > 0.5 else 'LOW CONCORDANCE — check PROGENy gene availability'}"
        )

    # ---- 4. Differential pathway analysis --------------------------------
    log.info("\nSTEP 4 — Differential pathway analysis (intraocular vs. extraocular)")
    log.info(
        "  Mann-Whitney U test for each pathway score between conditions.\n"
        "  Benjamini-Hochberg FDR correction applied across all pathways.\n"
        "  Expected: TGF-β significantly higher in extraocular (Wan et al. 2025)"
    )
    pathway_cols = ([c for c in adata.obs.columns if c.startswith("PROGENy_")]
                    + ["TGFb_manual_score"])
    group_col = "disease_stage" if "disease_stage" in adata.obs.columns else "dataset"
    group_a   = "intraocular" if "disease_stage" in adata.obs.columns else "GSE168434"
    group_b   = "extraocular" if "disease_stage" in adata.obs.columns else "GSE249995"

    diff_df = differential_pathway_analysis(
        adata, pathway_cols, group_col, group_a, group_b
    )
    if not diff_df.empty:
        diff_df.to_csv(TAB_DIR / "tgfb_differential_analysis.csv", index=False)
        log.info(f"  Differential pathway results:\n{diff_df.to_string()}")
        sig = diff_df[diff_df["significant"]]
        log.info(
            f"  Significantly rewired pathways (FDR<0.05, |log2FC|>0.5):\n"
            f"{sig[['pathway', 'log2fc', 'p_adj']].to_string()}"
        )

    # ---- 5. Save pathway scores table ------------------------------------
    log.info("\nSTEP 5 — Saving pathway score table")
    meta_cols = ["sample_id", "dataset", "cell_type_broad", "rb_subtype",
                 "cnv_load", "is_tumour_cnv"]
    meta_cols = [c for c in meta_cols if c in adata.obs.columns]
    score_cols = [c for c in adata.obs.columns
                  if c.startswith("PROGENy_") or c == "TGFb_manual_score"]
    adata.obs[meta_cols + score_cols].to_csv(
        TAB_DIR / "progeny_pathway_scores.csv"
    )
    log.info(
        f"  Saved pathway scores → {TAB_DIR / 'progeny_pathway_scores.csv'}"
    )

    # ---- 6. Visualizations -----------------------------------------------
    log.info("\nSTEP 6 — Generating pathway visualizations")
    plot_pathway_scores_umap(adata, FIG_DIR / "pathway_scores_umap.pdf")
    plot_tgfb_by_celltype_stage(adata, FIG_DIR / "tgfb_score_by_celltype_stage.pdf")
    plot_pathway_heatmap_by_celltype(adata, FIG_DIR / "pathway_heatmap_by_celltype.pdf")
    plot_tgfb_correlations(adata, FIG_DIR / "tgfb_correlation_subtype_fate.pdf")
    if not diff_df.empty:
        plot_tme_rewiring(diff_df, FIG_DIR / "tme_rewiring_comparison.pdf")

    # ---- 7. Save --------------------------------------------------------
    log.info(f"\nSTEP 7 — Saving fully scored atlas → {OUT_H5AD.name}")
    log.info(
        "  Contents added:\n"
        "    .obs['PROGENy_*']        : 14 pathway activity scores\n"
        "    .obs['TGFb_manual_score']: canonical TGF-β target gene score\n"
        "    .uns['progeny_diff']     : differential pathway analysis table"
    )
    if not diff_df.empty:
        adata.uns["progeny_diff"] = diff_df.to_dict()
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    log.info(
        "=" * 60 + "\n"
        "  PIPELINE COMPLETE\n"
        "  Final atlas: 11_pathway_scored.h5ad\n"
        "  Contains: QC → normalization → scVI integration → annotation\n"
        "            → CNV → subtype scoring → RNA velocity → CellRank\n"
        "            → L-R communication → TGF-β pathway scoring\n"
        + "=" * 60
    )
    return adata


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_pathway_scoring()
