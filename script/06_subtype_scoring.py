"""
=============================================================================
Script 06 — RB Molecular Subtype Scoring
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
Liu et al. (2021, Nature Communications) used single-cell transcriptomics
to define two molecular subtypes of retinoblastoma tumour cells:

  SUBTYPE 1 — Mature cone-like
    • Expresses cone photoreceptor maturation markers: ARR3, RXRG, OPN1SW,
      THRB, GNAT2
    • Resembles the post-mitotic cone precursor cell of origin (Singh et al.
      2018, PNAS)
    • More common in heritable (germline RB1 mutation) cases
    • Better-differentiated; lower proliferative index
    • Lower tendency for vitreous seeding or optic nerve invasion

  SUBTYPE 2 — Dedifferentiated / stemness
    • Expresses mesenchymal/stemness/early retinal progenitor markers:
      VIM, NES, CD44, SOX2, PROM1, FN1
    • Associated with MYCN amplification (chr2p)
    • Higher MYCN/E2F target gene activity
    • Higher metastatic risk (extraocular extension)
    • Present in both germline and somatic cases but enriched in high-risk
      patients

SIGNIFICANCE IN THIS PROJECT
------------------------------
By scoring every cell in the GSE168434 (intraocular) and GSE249995
(intraocular + extraocular) atlas with Subtype 1/2 signatures, we can:
  1. Test whether extraocular tumour cells are enriched for Subtype 2 scores
     (supporting the intraocular → extraocular dedifferentiation model).
  2. Trace subtype composition along the cone precursor trajectory identified
     by RNA velocity (script 07) to determine where dedifferentiation occurs.
  3. Correlate subtype scores with TGF-β activity (script 10) to identify
     whether TGF-β signaling drives the Subtype 1 → Subtype 2 transition.

ADDITIONAL SIGNATURES
----------------------
We also score cells with:
  - Cell cycle (G1/S and G2/M) markers — proliferation index
  - Hypoxia (HIF1A targets) — relevant to avascular intraocular tumour
  - Cone photoreceptor identity (ARR3+, high-confidence progenitor cells)
  - Retinal progenitor cell (RPC) signature (FGF/NOTCH targets)
  - TAM M1/M2 polarisation (for macrophage scripts 09/10)

METHOD
------
Scanpy's `sc.tl.score_genes()` implements a modified z-score method:
  For each cell, it computes the mean expression of the signature genes
  and subtracts the mean expression of a random control gene set drawn
  from the same expression bins (Tirosh et al. 2016, Science).
This produces a score that is robust to technical variation in library size.

REFERENCES
----------
- Liu J et al. Genome and transcriptome sequencing of retinoblastoma
  reveals increased mutational burden in the unaffected eye.
  Nat Commun. 2021;12:5744.
  https://doi.org/10.1038/s41467-021-25595-x

- Singh HP et al. Developmental stage-specific proliferation and
  retinoblastoma genesis in RB-deficient human and mouse retinas.
  PNAS. 2018;115(44):E10337-E10346.
  https://doi.org/10.1073/pnas.1808903115

- Tirosh I et al. Dissecting the multicellular ecosystem of metastatic
  melanoma by single-cell RNA-seq. Science. 2016;352(6282):189-196.
  https://doi.org/10.1126/science.aad0501

- Wan W et al. Single-cell transcriptome landscape of intraocular and
  extraocular retinoblastoma. Ophthalmology. 2025.
  https://doi.org/10.1016/j.ophtha.2025.01.011

INPUT  : data/processed/06_cnv.h5ad
OUTPUT :
  data/processed/07_subtyped.h5ad
  results/figures/subtype_scores_umap.pdf
  results/figures/subtype_violin_by_dataset.pdf
  results/figures/subtype_composition_bar.pdf
  results/tables/subtype_scores_per_cell.csv

USAGE  : pixi run python script/06_subtype_scoring.py
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
IN_H5AD  = ROOT / "data" / "processed" / "06_cnv.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "07_subtyped.h5ad"
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
# Gene signatures
# ---------------------------------------------------------------------------

SIGNATURES: dict[str, list[str]] = {

    # ---- RB tumour subtypes (Liu et al. 2021, Singh et al. 2018) -----
    "RB_subtype1_cone": [
        "ARR3", "RXRG", "OPN1SW", "OPN1MW", "THRB", "GNAT2",
        "PDE6H", "CNGB3", "CNGA3", "KCNV2", "RCVRN", "CRABP1",
    ],
    "RB_subtype2_stemness": [
        "VIM", "NES", "CD44", "SOX2", "PROM1", "FN1", "SNAI2",
        "TWIST1", "ZEB1", "MYCN", "HMGA1", "ID1", "ID3",
    ],

    # ---- Cone precursor identity (Wan et al. 2025, CP subpopulation) --
    "cone_precursor_CP4_invasive": [
        "SOX4", "MMP2", "MMP9", "TGFB1", "TGFB2", "TGFBR1",
        "SMAD3", "CTGF", "FN1", "CDH2",
    ],

    # ---- Cell cycle signatures (Tirosh et al. 2016) -------------------
    "cell_cycle_S": [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1",
        "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1",
        "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
        "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
        "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN",
        "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
        "BRIP1", "E2F8",
    ],
    "cell_cycle_G2M": [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
        "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
        "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
        "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
        "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
        "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
        "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
        "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
    ],

    # ---- Hypoxia (HIF1A targets) ----------------------------------------
    "hypoxia_HIF1A": [
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1",
        "ALDOA", "CA9", "BNIP3", "BNIP3L", "LOX", "P4HA1",
    ],

    # ---- Retinal progenitor cells (RPC) -----------------------------------
    "retinal_progenitor_RPC": [
        "VSX2", "PAX6", "SOX2", "LHX2", "ASCL1", "HES1", "DLL3",
        "NOTCH1", "PCNA", "TOP2A",
    ],

    # ---- TAM M1 polarisation (pro-inflammatory) -------------------------
    "TAM_M1": [
        "CD86", "IL1B", "IL6", "TNF", "CXCL10", "NOS2", "IRF5",
        "STAT1", "HLA-DRA", "CD64",
    ],

    # ---- TAM M2 polarisation (immunosuppressive) ------------------------
    "TAM_M2": [
        "CD163", "MRC1", "ARG1", "IL10", "TGFB1", "CCL18",
        "FOLR2", "LYVE1", "STAB1", "PDCD1LG2",
    ],
}


# ---------------------------------------------------------------------------
# Scoring and plotting helpers
# ---------------------------------------------------------------------------

def score_all_signatures(adata: sc.AnnData) -> None:
    """Score cells with all predefined gene signatures.

    SIGNIFICANCE:
    Each signature score is stored as a column in adata.obs.  The scores
    enable:
      - Comparing subtype composition between intraocular and extraocular
        samples (test of the invasion model)
      - Continuous scoring rather than binary classification — cells near
        the decision boundary are biologically ambiguous and should not be
        over-interpreted
      - Cell cycle assignment for pseudotime and velocity analyses
        (cells in S/G2M phase should be accounted for in RNA velocity)
    """
    for sig_name, genes in SIGNATURES.items():
        # Keep only genes present in the dataset
        present = [g for g in genes if g in adata.var_names]
        absent  = [g for g in genes if g not in adata.var_names]
        if absent:
            log.info(f"  {sig_name}: {len(present)}/{len(genes)} genes found "
                     f"(missing: {absent[:5]}{'...' if len(absent) > 5 else ''})")
        if len(present) < 3:
            log.warning(f"  {sig_name}: fewer than 3 genes present — score unreliable")
            adata.obs[f"score_{sig_name}"] = np.nan
            continue
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=f"score_{sig_name}",
            ctrl_size=50,
            n_bins=25,
            use_raw=False,
        )
        score_col = adata.obs[f"score_{sig_name}"]
        log.info(
            f"  {sig_name}: genes={len(present)}, "
            f"score range [{score_col.min():.3f}, {score_col.max():.3f}], "
            f"mean={score_col.mean():.3f}"
        )


def assign_cell_cycle_phase(adata: sc.AnnData) -> None:
    """Assign S/G2M/G1 phase using Scanpy's built-in cell cycle scoring.

    SIGNIFICANCE:
    Cell cycle phase assignment is critical for RNA velocity: cycling cells
    produce unspliced mRNA that, without phase information, can be
    misinterpreted as directed transcriptional change.  Knowing which cells
    are in S or G2M also allows separate analysis of cycling vs. quiescent
    tumour populations, and supports identification of proliferating
    subclusters linked to MYCN amplification (Subtype 2).
    """
    s_genes  = [g for g in SIGNATURES["cell_cycle_S"]   if g in adata.var_names]
    g2m_genes= [g for g in SIGNATURES["cell_cycle_G2M"] if g in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    phase_counts = adata.obs["phase"].value_counts()
    log.info(f"  Cell cycle phase distribution:\n{phase_counts.to_string()}")


def assign_rb_subtype(adata: sc.AnnData) -> None:
    """Assign predominant RB subtype to each tumour cell.

    Logic
    -----
    For cells annotated as tumour (Cone_precursor, Tumour_proliferating,
    Cone_precursor_*), compute the difference:
      subtype_score = score_RB_subtype1_cone - score_RB_subtype2_stemness
    Positive → Subtype 1 (cone-like)
    Negative → Subtype 2 (stemness)
    Near zero (|diff| < 0.1) → Mixed / unresolved

    Non-tumour cells are labelled 'Non-tumour'.
    """
    tumour_types = {"Cone_precursor", "Tumour_proliferating"}
    tumour_mask  = adata.obs["cell_type_broad"].isin(tumour_types)

    # Also include fine subclusters starting with "CP_" prefix
    if "cell_type_fine" in adata.obs.columns:
        fine_tumour = adata.obs["cell_type_fine"].str.startswith(("CP_", "TP_"))
        tumour_mask = tumour_mask | fine_tumour

    subtype = pd.Series("Non-tumour", index=adata.obs_names)
    s1      = adata.obs.get("score_RB_subtype1_cone",     pd.Series(np.nan, index=adata.obs_names))
    s2      = adata.obs.get("score_RB_subtype2_stemness", pd.Series(np.nan, index=adata.obs_names))

    diff = s1 - s2
    subtype[tumour_mask & (diff >  0.1)] = "Subtype1_cone"
    subtype[tumour_mask & (diff < -0.1)] = "Subtype2_stemness"
    subtype[tumour_mask & (diff.abs() <= 0.1)] = "Mixed"

    adata.obs["rb_subtype"] = subtype.astype("category")
    counts = subtype[tumour_mask].value_counts()
    log.info(f"  RB subtype assignment (tumour cells only):\n{counts.to_string()}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_scores_umap(adata: sc.AnnData, out_path: Path) -> None:
    """UMAP grid coloured by key signature scores."""
    score_keys = [
        "score_RB_subtype1_cone",
        "score_RB_subtype2_stemness",
        "score_cone_precursor_CP4_invasive",
        "score_hypoxia_HIF1A",
        "rb_subtype",
        "phase",
    ]
    present = [k for k in score_keys if k in adata.obs.columns]
    ncols   = 3
    nrows   = int(np.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes = axes.flatten() if nrows > 1 else [axes] if ncols == 1 else axes.flatten()

    for ax, key in zip(axes, present):
        sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                   size=2, alpha=0.7)
        ax.set_title(key.replace("score_", ""), fontsize=9)

    for ax in axes[len(present):]:
        ax.set_visible(False)

    plt.suptitle("Molecular subtype and signature scores", fontsize=12,
                 fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_subtype_violin(adata: sc.AnnData, out_path: Path) -> None:
    """Violin plot of subtype scores split by dataset/disease_stage."""
    score_keys = ["score_RB_subtype1_cone", "score_RB_subtype2_stemness"]
    present    = [k for k in score_keys if k in adata.obs.columns]
    if not present:
        return

    groupby = "disease_stage" if "disease_stage" in adata.obs.columns else "dataset"
    fig, axes = plt.subplots(1, len(present), figsize=(7 * len(present), 5))
    if len(present) == 1:
        axes = [axes]

    for ax, key in zip(axes, present):
        sc.pl.violin(
            adata, keys=key, groupby=groupby,
            ax=ax, show=False, rotation=30,
        )
        ax.set_title(key.replace("score_", ""), fontsize=10)
        ax.set_xlabel("")

    plt.suptitle(
        f"Subtype scores by {groupby}\n"
        "(Subtype 2 enrichment in extraocular expected)",
        fontsize=11, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


def plot_subtype_composition(adata: sc.AnnData, out_path: Path) -> None:
    """Stacked bar chart of RB subtype composition per sample."""
    if "rb_subtype" not in adata.obs.columns:
        return

    tumour_types = {"Cone_precursor", "Tumour_proliferating"}
    tumour_adata = adata[adata.obs["cell_type_broad"].isin(tumour_types)]
    if tumour_adata.n_obs == 0:
        return

    comp = (
        tumour_adata.obs
        .groupby(["sample_id", "rb_subtype"])
        .size()
        .unstack(fill_value=0)
    )
    comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100

    colours = {
        "Subtype1_cone":     "#4C9BE8",
        "Subtype2_stemness": "#E84C4C",
        "Mixed":             "#F5A623",
        "Non-tumour":        "#AAAAAA",
    }
    cols = [c for c in colours if c in comp_pct.columns]
    fig, ax = plt.subplots(figsize=(max(8, len(comp_pct) * 0.8), 5))
    comp_pct[cols].plot(
        kind="bar", stacked=True, ax=ax,
        color=[colours[c] for c in cols],
        edgecolor="none",
    )
    ax.set_xlabel("Sample ID", fontsize=10)
    ax.set_ylabel("% tumour cells", fontsize=10)
    ax.set_title(
        "RB subtype composition per sample\n(tumour cells only)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_subtype_scoring() -> sc.AnnData:
    """Molecular subtype scoring pipeline.

    Steps
    -----
    1. Load CNV-annotated atlas.
    2. Score all gene signatures (subtype, cell cycle, hypoxia, RPC, TAM).
    3. Assign cell cycle phase.
    4. Assign RB subtype (Subtype 1 / 2 / Mixed) for tumour cells.
    5. Visualize scores on UMAP.
    6. Plot subtype distribution by sample and disease stage.
    7. Save.
    """

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading CNV-annotated atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 2. Score all signatures ----------------------------------------
    log.info("\nSTEP 2 — Scoring gene signatures")
    log.info(
        "  Using sc.tl.score_genes() with expression-bin-matched control\n"
        "  gene sets (Tirosh et al. 2016). Robust to library-size effects."
    )
    score_all_signatures(adata)

    # ---- 3. Cell cycle phase assignment ---------------------------------
    log.info("\nSTEP 3 — Cell cycle phase assignment (S / G2M / G1)")
    log.info(
        "  WHY: Cell cycle phase will be used in RNA velocity (script 07)\n"
        "  to avoid confounding velocity vectors with cell division events.\n"
        "  Proliferating Subtype 2 cells are expected to be enriched in S/G2M."
    )
    assign_cell_cycle_phase(adata)

    # ---- 4. RB subtype assignment ---------------------------------------
    log.info("\nSTEP 4 — RB subtype assignment")
    log.info(
        "  Subtype 1 (cone-like):   score_RB_subtype1 > score_RB_subtype2 + 0.1\n"
        "  Subtype 2 (stemness):    score_RB_subtype2 > score_RB_subtype1 + 0.1\n"
        "  Mixed:                   |score difference| <= 0.1"
    )
    assign_rb_subtype(adata)

    # ---- 5. Visualize ---------------------------------------------------
    log.info("\nSTEP 5 — Generating subtype visualizations")
    plot_scores_umap(adata, FIG_DIR / "subtype_scores_umap.pdf")
    plot_subtype_violin(adata, FIG_DIR / "subtype_violin_by_dataset.pdf")
    plot_subtype_composition(adata, FIG_DIR / "subtype_composition_bar.pdf")

    # ---- 6. Save score table --------------------------------------------
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    meta_cols  = ["sample_id", "dataset", "cell_type_broad", "rb_subtype",
                  "phase", "cnv_load", "is_tumour_cnv"]
    meta_cols  = [c for c in meta_cols if c in adata.obs.columns]
    adata.obs[meta_cols + score_cols].to_csv(
        TAB_DIR / "subtype_scores_per_cell.csv"
    )
    log.info(
        f"  Saved subtype scores table → "
        f"{TAB_DIR / 'subtype_scores_per_cell.csv'}"
    )

    # ---- 7. Save --------------------------------------------------------
    log.info(f"\nSTEP 7 — Saving subtyped atlas → {OUT_H5AD.name}")
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    return adata


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_subtype_scoring()
