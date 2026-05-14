"""
=============================================================================
Script 04 — Cell Type Annotation and Tumor Subcluster Identification
=============================================================================
Project : RetinoblastomaAtlas
Author  : Md. Jubayer Hossain
Date    : 2026-04-29

WHY THIS ANALYSIS?
------------------
The scVI-integrated atlas contains Leiden clusters that partition cells
by transcriptional similarity, but the biological identity of each cluster
must be assigned by domain knowledge.  Cell type annotation is the
*interpretive keystone* of the atlas: every downstream analysis — trajectory
inference, cell-cell communication, differential abundance — depends on
correct and reproducible cell type labels.

Retinoblastoma has a complex cellular composition:

  MALIGNANT compartment:
  ─ Cone precursor-like (CP): express ARR3, RXRG, THRB, PRDM1, CRX.
    These are the primary cells of interest — they represent the tumour
    cell of origin and the cells undergoing invasion (Aim 3).
  ─ Retinoma-like (RL): quiescent pre-malignant cells expressing p16/CDKN2A,
    p130/RBL2, low Ki67.
  ─ MKI67-high proliferating (MKI67PhrD): high MKI67, TOP2A, BIRC5, UBE2C.
    Elevated in extraocular samples (Liu et al. 2024, Commun Biol).

  MICROENVIRONMENT compartment:
  ─ Müller glia / Astrocyte-like: GFAP, VIM, SOX9, RLBP1
  ─ Microglia / TAM: IBA1, CD68, AIF1, CX3CR1
  ─ Vascular endothelium: PECAM1, CDH5, VWF
  ─ Pericyte: PDGFRB, ACTA2, RGS5
  ─ Fibroblast: DCN, COL1A1, LUM
  ─ T cell: CD3D, CD8A, CD4
  ─ NK / ILC: GNLY, NKG7, KLRD1
  ─ B cell / Plasma: CD79A, MS4A1, MZB1

SIGNIFICANCE OF FINE SUBCLUSTERING
------------------------------------
Liu et al. (2024, Commun Biol) resolved 10 cone precursor subpopulations
and identified MKI67PhrD cells as the primary extraocular-enriched
subtype driving SOX4-mediated local extension.  Wan et al. (2025, Sci Rep)
found that CP subcluster 4 shows elevated TGF-β signaling specifically in
invasive RB.  Reproducing and extending these subclusters in a larger
cohort is a primary deliverable of Aim 2.

ANNOTATION STRATEGY
--------------------
1. Score known marker genes (AUCell-like per-cell gene set scoring).
2. Assign broad cell types using sc.tl.score_genes().
3. Validate with dot plots and heatmaps.
4. Fine-subcluster the malignant and TAM compartments independently.
5. Compute differential proportions between intraocular and extraocular.

REFERENCES
----------
- Liu Y et al. Single-cell transcriptomics enable the characterization
  of local extension in retinoblastoma. Commun Biol. 2024;7:11.
  https://doi.org/10.1038/s42003-023-05732-y

- Yang J et al. Single-cell transcriptome profiling reveals intratumoural
  heterogeneity and malignant progression in retinoblastoma. Cell Death Dis.
  2021;12:1100. https://doi.org/10.1038/s41419-021-04390-4

- Wan W et al. Sci Rep. 2025;15:39954.
  https://doi.org/10.1038/s41598-025-23779-1

INPUT  : data/processed/04_scvi_integrated.h5ad
OUTPUT :
  data/processed/05_annotated.h5ad
  results/figures/umap_cell_types_broad.pdf
  results/figures/umap_cell_types_fine.pdf
  results/figures/dotplot_marker_genes.pdf
  results/figures/umap_disease_stage.pdf
  results/tables/cell_type_proportions.csv

USAGE  : pixi run python script/04_cell_type_annotation.py
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
IN_H5AD  = ROOT / "data" / "processed" / "04_scvi_integrated.h5ad"
OUT_H5AD = ROOT / "data" / "processed" / "05_annotated.h5ad"
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
# Marker gene dictionaries
# Compiled from Liu et al. 2024 (Commun Biol), Yang et al. 2021 (Cell Death Dis),
# and the Human Retina Cell Atlas.
# ---------------------------------------------------------------------------

BROAD_MARKERS: dict[str, list[str]] = {
    "Cone_precursor": ["ARR3", "RXRG", "THRB", "PRDM1", "CRX", "RCVRN",
                       "OPN1MW", "OPN1LW", "GNB3"],
    "Retinoma_like":  ["CDKN2A", "RBL2", "GADD45A"],
    "MKI67_high_RB":  ["MKI67", "TOP2A", "BIRC5", "UBE2C", "CENPF", "NUSAP1"],
    "Muller_glia":    ["GFAP", "VIM", "SOX9", "RLBP1", "CLU", "GLUL"],
    "Microglia_TAM":  ["AIF1", "CX3CR1", "PTPRC", "CD68", "IBA1", "CSF1R"],
    "Endothelium":    ["PECAM1", "CDH5", "VWF", "PLVAP", "KDR"],
    "Pericyte":       ["PDGFRB", "ACTA2", "RGS5", "CSPG4", "NOTCH3"],
    "Fibroblast":     ["DCN", "COL1A1", "LUM", "POSTN", "THY1"],
    "T_cell":         ["CD3D", "CD3E", "CD8A", "CD4", "TRAC"],
    "NK_ILC":         ["GNLY", "NKG7", "KLRD1", "NCR1", "FCGR3A"],
    "B_Plasma":       ["CD79A", "MS4A1", "MZB1", "JCHAIN", "IGKC"],
}

# Fine markers for TAM polarization (used in Aim 4)
TAM_FINE_MARKERS: dict[str, list[str]] = {
    "M1_like_TAM":  ["CXCL10", "CXCL9", "IDO1", "GBP1", "STAT1"],
    "M2_like_TAM":  ["CD163", "MRC1", "MIF", "CCL18", "FOLR2", "TREM2"],
    "Transiting_TAM": ["SPP1", "FN1", "MMP9", "HMOX1"],
}

# Cone precursor fine markers (Liu et al. 2024 Table 1)
CP_FINE_MARKERS: dict[str, list[str]] = {
    "CP_mature":       ["ARR3", "OPN1LW", "OPN1MW", "GUCA1C"],
    "CP_intermediate": ["RXRG", "THRB", "CRX", "PRDM1"],
    "CP_immature":     ["CCND1", "E2F1", "MYCN", "LIN28B"],
    "CP_SOX4_high":    ["SOX4", "TWIST1", "SNAI2"],   # extraocular enriched
    "CP_proliferating":["MKI67", "TOP2A", "UBE2C", "BIRC5"],
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def score_cell_types(adata: sc.AnnData, marker_dict: dict) -> sc.AnnData:
    """Score each cell for gene set activity using Scanpy's score_genes.

    WHY SCORE_GENES?
    ----------------
    sc.tl.score_genes computes a per-cell enrichment score for a gene list
    relative to a control set of randomly sampled genes with similar
    expression levels.  This is analogous to AUCell but computationally
    cheaper and sufficient for broad annotation.  The resulting scores are
    stored in .obs and used to assign cell type labels.
    """
    for ct, genes in marker_dict.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) < 2:
            log.warning(f"  {ct}: only {len(genes_present)} marker(s) found — skipping score")
            continue
        sc.tl.score_genes(adata, genes_present, score_name=f"score_{ct}")
        log.info(f"  Scored: {ct} ({len(genes_present)}/{len(genes)} markers found)")
    return adata


def assign_broad_cell_type(adata: sc.AnnData) -> sc.AnnData:
    """Assign broad cell type label based on highest score.

    The cell type with the maximum score_* value is assigned.  Cells where
    all scores are < 0.1 are labelled 'Unknown'.
    """
    score_cols = [f"score_{ct}" for ct in BROAD_MARKERS if f"score_{ct}" in adata.obs.columns]
    if not score_cols:
        log.error("No score columns found — run score_cell_types first")
        return adata

    score_matrix = adata.obs[score_cols].copy()
    score_matrix.columns = [c.replace("score_", "") for c in score_cols]

    max_scores = score_matrix.max(axis=1)
    max_types  = score_matrix.idxmax(axis=1)
    max_types[max_scores < 0.1] = "Unknown"

    adata.obs["cell_type_broad"] = pd.Categorical(max_types)
    log.info(f"  Broad cell type distribution:\n{adata.obs['cell_type_broad'].value_counts().to_string()}")
    return adata


def fine_subcluster(
    adata: sc.AnnData,
    cell_type: str,
    broad_key: str = "cell_type_broad",
    resolution: float = 1.5,
    key_added: str | None = None,
) -> sc.AnnData:
    """Sub-cluster a specific cell type compartment.

    WHY FINE SUBCLUSTERING?
    -----------------------
    Broad Leiden clustering at atlas level uses a resolution tuned for
    distinguishing major cell types.  Within-compartment biology requires
    higher resolution clustering on the subset, using the same scVI latent
    space but focusing neighbour distances on the relevant cells.
    """
    if key_added is None:
        key_added = f"leiden_fine_{cell_type.lower()}"

    mask = adata.obs[broad_key] == cell_type
    n_cells = mask.sum()
    log.info(f"  Fine subclustering {cell_type}: {n_cells:,} cells, resolution={resolution}")

    if n_cells < 50:
        log.warning(f"  Too few {cell_type} cells ({n_cells}) for subclustering — skipped")
        return adata

    adata_sub = adata[mask].copy()
    sc.pp.neighbors(adata_sub, use_rep="X_scVI", n_neighbors=20)
    sc.tl.leiden(adata_sub, resolution=resolution, key_added="leiden_fine")

    # Map labels back to main adata
    adata.obs[key_added] = "Other"
    adata.obs.loc[mask, key_added] = (
        f"{cell_type}_" + adata_sub.obs["leiden_fine"].astype(str)
    )
    adata.obs[key_added] = pd.Categorical(adata.obs[key_added])
    n_sub = adata_sub.obs["leiden_fine"].nunique()
    log.info(f"  → {n_sub} {cell_type} subclusters identified")
    return adata


def compute_cell_proportions(adata: sc.AnnData, out_path: Path) -> pd.DataFrame:
    """Compute cell type proportions per sample and disease stage.

    SIGNIFICANCE:
    Changes in cell type proportions between intraocular and extraocular
    tumours reflect the compositional remodelling of the tumour and its
    microenvironment during invasion.  For example, an increase in the
    proportion of MKI67-high proliferating cells or a decrease in mature
    cone-like cells in extraocular samples would support the dedifferentiation
    model of RB local extension.
    """
    if "disease_stage" not in adata.obs.columns:
        # Derive from sample IDs
        stage_map = {
            "S1_in1": "intraocular", "S2_in2": "intraocular",
            "S3_ex1": "extraocular", "S4_ex2": "extraocular",
        }
        adata.obs["disease_stage"] = adata.obs["sample_id"].map(
            lambda x: stage_map.get(x, "intraocular")  # GSE168434 are intraocular
        )

    prop = (
        adata.obs.groupby(["sample_id", "disease_stage", "cell_type_broad"])
        .size()
        .reset_index(name="n_cells")
    )
    total = prop.groupby("sample_id")["n_cells"].transform("sum")
    prop["proportion"] = prop["n_cells"] / total
    prop.to_csv(out_path, index=False)
    log.info(f"  Saved cell proportions → {out_path.name}")
    return prop


def plot_dotplot_markers(adata: sc.AnnData, out_path: Path) -> None:
    """Dot plot of canonical markers per broad cell type.

    SIGNIFICANCE:
    This is the primary evidence figure for cell type annotation — shows
    that each annotated cluster expresses its expected markers and not the
    markers of other cell types (specificity check).
    """
    selected_markers = {
        ct: genes[:4] for ct, genes in BROAD_MARKERS.items()
    }
    all_genes = [g for gs in selected_markers.values() for g in gs
                 if g in adata.var_names]
    if "cell_type_broad" not in adata.obs.columns:
        log.warning("cell_type_broad not found; skipping dot plot")
        return

    fig, ax = plt.subplots(figsize=(18, 8))
    sc.pl.dotplot(
        adata,
        var_names=all_genes,
        groupby="cell_type_broad",
        standard_scale="var",
        cmap="RdBu_r",
        dot_min=0.05,
        dot_max=1.0,
        ax=ax,
        show=False,
    )
    plt.suptitle("Canonical marker genes per cell type\n(dot size = % expressing, colour = scaled expression)",
                 fontsize=10, y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: {out_path.name}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_annotation() -> sc.AnnData:
    """Cell type annotation pipeline."""

    # ---- 1. Load --------------------------------------------------------
    log.info("=" * 60)
    log.info("STEP 1 — Loading scVI-integrated atlas")
    log.info("=" * 60)
    adata = sc.read_h5ad(IN_H5AD)
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 2. Score broad markers ----------------------------------------
    log.info("\nSTEP 2 — Scoring broad cell type markers")
    log.info(
        "  WHY: Gene set scoring provides a continuous, per-cell estimate\n"
        "  of each cell type's transcriptional programme. It is less\n"
        "  dependent on clustering resolution than looking at cluster-level\n"
        "  marker genes and allows cells to have partial identities."
    )
    adata = score_cell_types(adata, BROAD_MARKERS)

    # ---- 3. Assign broad cell types ------------------------------------
    log.info("\nSTEP 3 — Assigning broad cell type labels")
    adata = assign_broad_cell_type(adata)

    # UMAP coloured by cell type
    sc.pl.umap(adata, color="cell_type_broad", show=False,
               save=None, frameon=False, size=2)
    fig = plt.gcf()
    fig.savefig(FIG_DIR / "umap_cell_types_broad.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved: umap_cell_types_broad.pdf")

    # ---- 4. Fine subclustering: Cone precursors -------------------------
    log.info("\nSTEP 4 — Fine subclustering: Cone precursor-like cells")
    log.info(
        "  WHY: CP subclusters correspond to distinct maturation states\n"
        "  and invasion phenotypes (Liu et al. 2024 identified 10 CP\n"
        "  subpopulations). High-resolution CP subclustering is required\n"
        "  for the RNA velocity analysis (Aim 3)."
    )
    adata = fine_subcluster(adata, "Cone_precursor", resolution=1.5)
    adata = score_cell_types(adata, CP_FINE_MARKERS)

    # ---- 5. Fine subclustering: TAMs -----------------------------------
    log.info("\nSTEP 5 — Fine subclustering: Microglia/TAM cells")
    log.info(
        "  WHY: TAM subclusters separate M1-like (anti-tumour) from\n"
        "  M2-like (pro-tumour, immunosuppressive) and transitional\n"
        "  subpopulations. Zhang et al. (2025) showed TAM subcluster MG1\n"
        "  is specifically enriched in extraocular RB."
    )
    adata = fine_subcluster(adata, "Microglia_TAM", resolution=1.0)
    adata = score_cell_types(adata, TAM_FINE_MARKERS)

    # ---- 6. Disease stage label ----------------------------------------
    log.info("\nSTEP 6 — Adding disease stage metadata")
    stage_map = {
        "S1_in1": "intraocular", "S2_in2": "intraocular",
        "S3_ex1": "extraocular", "S4_ex2": "extraocular",
    }
    adata.obs["disease_stage"] = adata.obs["sample_id"].map(
        lambda x: stage_map.get(x, "intraocular")
    )
    log.info(f"  Stage distribution:\n{adata.obs['disease_stage'].value_counts().to_string()}")

    sc.pl.umap(adata, color="disease_stage", show=False, frameon=False, size=2)
    fig = plt.gcf()
    fig.savefig(FIG_DIR / "umap_disease_stage.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- 7. Marker gene dot plot ---------------------------------------
    log.info("\nSTEP 7 — Generating marker gene dot plot")
    plot_dotplot_markers(adata, FIG_DIR / "dotplot_marker_genes.pdf")

    # ---- 8. Cell type proportions per sample ---------------------------
    log.info("\nSTEP 8 — Computing cell type proportions")
    prop_df = compute_cell_proportions(adata, TAB_DIR / "cell_type_proportions.csv")

    # Bar chart of proportions by disease stage
    stage_prop = (
        prop_df.groupby(["disease_stage", "cell_type_broad"])["proportion"]
        .mean()
        .unstack(fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(12, 5))
    stage_prop.T.plot(kind="bar", ax=ax, colormap="tab20", edgecolor="none", width=0.8)
    ax.set_xlabel("Cell type", fontsize=11)
    ax.set_ylabel("Mean proportion", fontsize=11)
    ax.set_title("Cell type composition: intraocular vs. extraocular", fontsize=11)
    ax.legend(title="Stage", bbox_to_anchor=(1, 1), fontsize=9)
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "cell_type_proportions_bar.pdf", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- 9. Compute marker genes per cluster (for Supplementary) --------
    log.info("\nSTEP 9 — Computing Leiden cluster marker genes (Wilcoxon)")
    log.info(
        "  WHY WILCOXON: The Wilcoxon rank-sum test is the recommended\n"
        "  DE test for scRNA-seq data in Scanpy (Zimmermann et al. 2021)\n"
        "  because it is non-parametric and robust to the zero-inflation\n"
        "  of count data."
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type_broad",
        method="wilcoxon",
        n_genes=50,
        use_raw=False,
        layer="scvi_normalized" if "scvi_normalized" in adata.layers else None,
    )
    marker_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_genes_df.to_csv(TAB_DIR / "marker_genes_per_cell_type.csv", index=False)
    log.info(f"  Saved marker genes → {TAB_DIR / 'marker_genes_per_cell_type.csv'}")

    # ---- 10. Save -------------------------------------------------------
    log.info(f"\nSTEP 10 — Saving annotated atlas → {OUT_H5AD.name}")
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    log.info("  Done.\n")
    return adata


if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    run_annotation()
