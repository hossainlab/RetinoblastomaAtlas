"""
qc.py
=====
QC metrics, evidence-based filter decisions, cell/gene filtering,
sparse validation, export, and full pipeline.

Functions
---------
Metrics
    compute_qc_metrics_human   — MT/ribo/HB metrics for human data
    compute_qc_metrics_mouse   — MT/ribo/HB metrics for mouse data

Evidence decisions (per sample)
    decide_mt_filter           — should MT% filter be applied?
    decide_ribo_filter         — should ribo% filter be applied?
    decide_doublet_filter      — should doublets be removed?
    build_decision_table       — run all three tests, compile master table

Filtering
    check_sample_sizes         — find samples too small for MAD QC
    apply_mad_filters          — per-sample MAD-based cell filtering
    check_filter_impact        — warn if > 30% cells removed per sample
    filter_genes               — remove low-detection genes

Validation + export
    validate_sparse            — assert CSR format, report storage stats
    export_h5ad                — save filtered AnnData as gzip h5ad

Full pipeline
    run_pipeline               — end-to-end: extract → load → QC → export

Example
-------
    from utils.qc import (
        compute_qc_metrics_human,
        build_decision_table,
        apply_mad_filters,
        filter_genes,
        validate_sparse,
        export_h5ad,
    )
"""

import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.mixture import GaussianMixture

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0


# ─────────────────────────────────────────────────────────────────────────────
# PRIVATE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

_SPECIES_PATTERNS = {
    "human": {
        "mt":   {"startswith": "MT-"},
        "ribo": {"contains":   r"^RPS|^RPL"},
        "hb":   {"contains":   r"^HB[^BP]"},
    },
    "mouse": {
        "mt":   {"startswith": "mt-"},
        "ribo": {"contains":   r"^Rps|^Rpl"},
        "hb":   {"contains":   r"^Hb[^bp]"},
    },
}


def _mad(series: pd.Series) -> float:
    """Median absolute deviation."""
    return (series - series.median()).abs().median()


def _mad_thresholds(series: pd.Series, nmads: int) -> tuple:
    med = series.median()
    m   = _mad(series)
    return med - nmads * m, med + nmads * m


def _compute_qc_metrics_core(adata: ad.AnnData, species: str) -> None:
    """Shared QC computation. Operates in-place."""
    p = _SPECIES_PATTERNS[species]

    adata.var["mt"] = (
        adata.var_names.str.startswith(p["mt"]["startswith"])
        if "startswith" in p["mt"]
        else adata.var_names.str.contains(p["mt"]["contains"], regex=True)
    )
    adata.var["ribo"] = adata.var_names.str.contains(p["ribo"]["contains"], regex=True)
    adata.var["hb"]   = adata.var_names.str.contains(p["hb"]["contains"],   regex=True)

    print(f"  [{species}] MT genes:   {adata.var['mt'].sum()}")
    print(f"  [{species}] Ribo genes: {adata.var['ribo'].sum()}")
    print(f"  [{species}] HB genes:   {adata.var['hb'].sum()}")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"],
                               percent_top=None, log1p=False, inplace=True)

    adata.obs["log1p_total_counts"]      = np.log1p(adata.obs["total_counts"])
    adata.obs["log1p_n_genes_by_counts"] = np.log1p(adata.obs["n_genes_by_counts"])
    adata.obs["log1p_genes_per_count"]   = np.log1p(
        adata.obs["n_genes_by_counts"] / adata.obs["total_counts"].clip(lower=1)
    )
    print(f"[metrics] Done — {adata.n_obs:,} cells")
    _print_qc_summary(adata)


def _print_qc_summary(adata: ad.AnnData) -> None:
    cols = [c for c in ["total_counts","n_genes_by_counts","pct_counts_mt","pct_counts_ribo"]
            if c in adata.obs.columns]
    df   = adata.obs[cols].describe().loc[["mean","50%","min","max"]]
    df.index = ["mean","median","min","max"]
    print(df.round(2).to_string())


def _outlier_per_sample(adata, metric, nmads, direction, groupby) -> pd.Series:
    """Per-sample MAD outlier flag. direction: 'high' | 'low' | 'both'."""
    outlier = pd.Series(False, index=adata.obs_names)
    for smp in adata.obs[groupby].unique():
        mask  = adata.obs[groupby] == smp
        vals  = adata.obs.loc[mask, metric]
        lo, hi = _mad_thresholds(vals, nmads)
        if   direction == "high": outlier[mask] = vals > hi
        elif direction == "low":  outlier[mask] = vals < lo
        else:                     outlier[mask] = (vals < lo) | (vals > hi)
    return outlier


# ─────────────────────────────────────────────────────────────────────────────
# QC METRICS
# ─────────────────────────────────────────────────────────────────────────────

def compute_qc_metrics_human(adata: ad.AnnData) -> None:
    """
    Compute QC metrics for human scRNA-seq data (in-place).

    Gene patterns: MT- | RPS/RPL | HB[^BP]

    Adds to adata.var
    -----------------
    mt, ribo, hb  (bool flags)

    Adds to adata.obs
    -----------------
    total_counts, n_genes_by_counts,
    pct_counts_mt, pct_counts_ribo, pct_counts_hb,
    log1p_total_counts, log1p_n_genes_by_counts, log1p_genes_per_count

    Example
    -------
    >>> compute_qc_metrics_human(adata)
    # [human] MT genes: 13, Ribo genes: 87, HB genes: 6
    """
    print("[metrics] Species: HUMAN")
    _compute_qc_metrics_core(adata, "human")


def compute_qc_metrics_mouse(adata: ad.AnnData) -> None:
    """
    Compute QC metrics for mouse scRNA-seq data (in-place).

    Gene patterns: mt- | Rps/Rpl | Hb[^bp]

    Adds to adata.var
    -----------------
    mt, ribo, hb  (bool flags)

    Adds to adata.obs
    -----------------
    total_counts, n_genes_by_counts,
    pct_counts_mt, pct_counts_ribo, pct_counts_hb,
    log1p_total_counts, log1p_n_genes_by_counts, log1p_genes_per_count

    Example
    -------
    >>> compute_qc_metrics_mouse(adata)
    # [mouse] MT genes: 13, Ribo genes: 83, HB genes: 5
    """
    print("[metrics] Species: MOUSE")
    _compute_qc_metrics_core(adata, "mouse")


# ─────────────────────────────────────────────────────────────────────────────
# EVIDENCE-BASED FILTER DECISIONS
# ─────────────────────────────────────────────────────────────────────────────

def decide_mt_filter(adata: ad.AnnData,
                     groupby: str = "sample",
                     nmads:   int = 3) -> pd.DataFrame:
    """
    Per-sample evidence test: should MT% be used for filtering?

    FILTER_MT = True only when ALL 3 criteria hold:
      1. Right-skewed distribution   (skewness > 1)
      2. High-MT outlier subpop      (> 1% of cells beyond +nmads MAD)
      3. High-MT cells have fewer genes  (Spearman r < -0.2, p < 0.05)

    Parameters
    ----------
    adata   : AnnData with pct_counts_mt in obs.
    groupby : per-sample column (default 'sample').
    nmads   : MADs above median defining high-MT outlier (default 3).

    Returns
    -------
    pd.DataFrame  One row per sample. Key columns: FILTER_MT, reason.

    Example
    -------
    >>> mt_rep = decide_mt_filter(adata)
    >>> mt_rep[['FILTER_MT', 'reason']]
    """
    report: dict = {}
    for smp in adata.obs[groupby].unique():
        mask  = adata.obs[groupby] == smp
        mt    = adata.obs.loc[mask, "pct_counts_mt"]
        genes = adata.obs.loc[mask, "n_genes_by_counts"]

        skewness    = stats.skew(mt)
        _, hi       = _mad_thresholds(mt, nmads)
        outlier_pct = ((mt > hi).sum() / len(mt)) * 100
        corr, pval  = stats.spearmanr(mt, genes)

        verdict = (skewness > 1) and (outlier_pct > 1.0) and (corr < -0.2) and (pval < 0.05)
        report[smp] = {
            "mt_skewness":  round(skewness,    3),
            "high_mt_pct":  round(outlier_pct, 2),
            "mt_gene_corr": round(corr,        3),
            "corr_pval":    round(pval,        4),
            "FILTER_MT":    verdict,
            "reason": ("skewed + outliers + dying-cell pattern"
                       if verdict else "no evidence of dead-cell enrichment"),
        }
    return pd.DataFrame(report).T


def decide_ribo_filter(adata: ad.AnnData,
                        groupby: str = "sample",
                        nmads:   int = 5) -> pd.DataFrame:
    """
    Per-sample evidence test: should ribo% be used for filtering?

    FILTER_RIBO = True only when ALL 3 criteria hold:
      1. Bimodal distribution    (GMM mode gap > 15%)
      2. High-ribo cells have fewer counts  (Spearman r < -0.3)
      3. Outlier extent > 0.5%

    Bimodality check prevents filtering genuinely ribo-high cell types
    (hepatocytes, erythrocytes).

    Parameters
    ----------
    adata   : AnnData with pct_counts_ribo in obs.
    groupby : per-sample column (default 'sample').
    nmads   : MADs for ribo outlier definition (default 5).

    Returns
    -------
    pd.DataFrame  One row per sample. Key columns: FILTER_RIBO, reason.

    Example
    -------
    >>> ribo_rep = decide_ribo_filter(adata)
    >>> ribo_rep[['bimodal', 'FILTER_RIBO', 'reason']]
    """
    report: dict = {}
    for smp in adata.obs[groupby].unique():
        mask   = adata.obs[groupby] == smp
        ribo   = adata.obs.loc[mask, "pct_counts_ribo"]
        counts = adata.obs.loc[mask, "total_counts"]

        gm      = GaussianMixture(n_components=2, random_state=0).fit(ribo.values.reshape(-1, 1))
        means   = sorted(gm.means_.flatten())
        bimodal = (means[1] - means[0]) > 15

        corr, _     = stats.spearmanr(ribo, counts)
        _, hi       = _mad_thresholds(ribo, nmads)
        outlier_pct = ((ribo > hi).sum() / len(ribo)) * 100

        verdict = bimodal and (corr < -0.3) and (outlier_pct > 0.5)
        report[smp] = {
            "ribo_mode_low":  round(means[0],    2),
            "ribo_mode_high": round(means[1],    2),
            "bimodal":        bimodal,
            "ribo_corr":      round(corr,        3),
            "outlier_pct":    round(outlier_pct, 2),
            "FILTER_RIBO":    verdict,
            "reason": ("bimodal + low-count ribo-high cells"
                       if verdict else "no ribo subpopulation detected"),
        }
    return pd.DataFrame(report).T


def decide_doublet_filter(adata: ad.AnnData,
                           groupby: str = "sample") -> pd.DataFrame:
    """
    Per-sample Scrublet evidence test: should doublets be removed?

    FILTER_DOUBLETS = True only when ALL 3 criteria hold:
      1. Score separation > 0.2       (real signal)
      2. Observed rate < 2× expected  (plausible for 10x)
      3. Doublet count enrichment ≥ 1.5×  (hallmark)

    Side effect: adds 'doublet_score' and 'predicted_doublet' to adata.obs.

    Parameters
    ----------
    adata   : AnnData with layers['counts'] or .X as raw counts.
    groupby : per-sample column (default 'sample').

    Returns
    -------
    pd.DataFrame  One row per sample. Empty if scrublet not installed.

    Example
    -------
    >>> dbl_rep = decide_doublet_filter(adata)
    >>> dbl_rep[['n_cells', 'observed_dbl_rate', 'FILTER_DOUBLETS']]
    """
    try:
        import scrublet as scr
    except ImportError:
        print("[doublet] scrublet not installed — skipping. pip install scrublet")
        return pd.DataFrame()

    report: dict = {}
    for smp in adata.obs[groupby].unique():
        mask          = adata.obs[groupby] == smp
        sub           = adata[mask].copy()
        X             = sub.layers["counts"] if "counts" in sub.layers else sub.X
        scrub         = scr.Scrublet(X, expected_doublet_rate=0.06)
        scores, calls = scrub.scrub_doublets(min_counts=2, min_cells=3, verbose=False)

        adata.obs.loc[mask, "doublet_score"]     = scores
        adata.obs.loc[mask, "predicted_doublet"] = calls

        n_cells    = mask.sum()
        exp_rate   = min(0.008 * (n_cells / 1000), 0.20)
        obs_rate   = calls.mean()
        separation = scores[calls].mean() - scores[~calls].mean() if calls.any() else 0
        dbl_med    = adata.obs.loc[mask, "total_counts"][calls].median()
        sgl_med    = adata.obs.loc[mask, "total_counts"][~calls].median()
        enrich     = (dbl_med / sgl_med) if sgl_med > 0 else 0

        verdict = (separation > 0.2) and (obs_rate < exp_rate * 2) and (enrich > 1.5)
        report[smp] = {
            "n_cells":           n_cells,
            "expected_dbl_rate": round(exp_rate,   4),
            "observed_dbl_rate": round(obs_rate,   4),
            "score_separation":  round(separation, 3),
            "count_enrichment":  round(enrich,     3),
            "FILTER_DOUBLETS":   verdict,
            "reason": ("clear doublet signal + plausible rate + count enrichment"
                       if verdict else "weak doublet signal"),
        }
    return pd.DataFrame(report).T


def build_decision_table(adata: ad.AnnData,
                          groupby: str = "sample") -> tuple:
    """
    Run all three evidence tests and compile master per-sample decision table.

    Parameters
    ----------
    adata   : AnnData with QC metrics computed.
    groupby : per-sample column (default 'sample').

    Returns
    -------
    tuple : (decision, mt_report, ribo_report, dbl_report)
        decision — pd.DataFrame, one row per sample:
                   filter_mt, mt_reason, filter_ribo, ribo_reason,
                   filter_doublets, doublet_reason

    Example
    -------
    >>> decision, mt_rep, ribo_rep, dbl_rep = build_decision_table(adata)
    >>> decision[['filter_mt', 'filter_ribo', 'filter_doublets']]
    #                        filter_mt  filter_ribo  filter_doublets
    # GSM5139852_RB01_rep1      True         True             True
    # GSM5139853_RB01_rep2      True        False             True
    # GSM5139854_RB02_rep1     False        False            False
    """
    print("[decisions] Running per-sample evidence tests...")
    mt_rep   = decide_mt_filter(adata,      groupby=groupby)
    ribo_rep = decide_ribo_filter(adata,    groupby=groupby)
    dbl_rep  = decide_doublet_filter(adata, groupby=groupby)

    decision = pd.DataFrame({
        "filter_mt":       mt_rep["FILTER_MT"],
        "mt_reason":       mt_rep["reason"],
        "filter_ribo":     ribo_rep["FILTER_RIBO"]    if not ribo_rep.empty else False,
        "ribo_reason":     ribo_rep["reason"]          if not ribo_rep.empty else "skipped",
        "filter_doublets": dbl_rep["FILTER_DOUBLETS"]  if not dbl_rep.empty else False,
        "doublet_reason":  dbl_rep["reason"]           if not dbl_rep.empty else "skipped",
    })

    print("\n" + "=" * 70)
    print("FILTER DECISION TABLE")
    print("=" * 70)
    print(decision[["filter_mt", "filter_ribo", "filter_doublets"]].to_string())
    print("=" * 70)
    return decision, mt_rep, ribo_rep, dbl_rep


# ─────────────────────────────────────────────────────────────────────────────
# FILTERING
# ─────────────────────────────────────────────────────────────────────────────

def check_sample_sizes(adata: ad.AnnData,
                        groupby:   str = "sample",
                        min_cells: int = 50) -> tuple:
    """
    Identify samples too small for reliable MAD-based QC (< min_cells).

    MAD is unstable below ~50 cells. Small samples should be excluded
    from MAD filtering or merged before QC.

    Parameters
    ----------
    adata     : AnnData
    groupby   : per-sample column (default 'sample').
    min_cells : minimum cells required for MAD QC (default 50).

    Returns
    -------
    tuple : (ok_samples, small_samples)

    Example
    -------
    >>> ok, small = check_sample_sizes(adata)
    >>> if small:
    ...     adata = adata[adata.obs['sample'].isin(ok)].copy()
    """
    counts      = adata.obs.groupby(groupby).size()
    ok_samples  = counts[counts >= min_cells].index.tolist()
    small_samps = counts[counts <  min_cells].index.tolist()

    print(f"[sample_sizes] Threshold: {min_cells} cells")
    print(f"  OK:    {len(ok_samples)} samples")
    print(f"  Small: {len(small_samps)} samples")
    if small_samps:
        for s in small_samps:
            print(f"    {s}: {counts[s]} cells")
    return ok_samples, small_samps


def apply_mad_filters(adata: ad.AnnData,
                       decision: pd.DataFrame,
                       groupby:       str = "sample",
                       nmads_counts:  int = 5,
                       nmads_genes:   int = 5,
                       nmads_mt:      int = 3,
                       nmads_ribo:    int = 5) -> ad.AnnData:
    """
    Apply per-sample MAD-based cell filters guided by the decision table.

    Always filters
    --------------
    Low total counts (empty droplets), high total counts (multiplets),
    low gene count (poor-quality cells).

    Evidence-gated (only where decision table = True)
    -------------------------------------------------
    High MT%     → FILTER_MT = True per sample
    High ribo%   → FILTER_RIBO = True per sample
    Doublets     → FILTER_DOUBLETS = True per sample

    Parameters
    ----------
    adata         : AnnData with QC metrics computed.
    decision      : DataFrame from build_decision_table().
    groupby       : per-sample column.
    nmads_counts  : MADs for count thresholds (default 5).
    nmads_genes   : MADs for gene count threshold (default 5).
    nmads_mt      : MADs for MT% threshold (default 3, tighter).
    nmads_ribo    : MADs for ribo% threshold (default 5).

    Returns
    -------
    AnnData  Filtered copy.

    Example
    -------
    >>> adata_filt = apply_mad_filters(adata, decision)
    # Low counts removed:    312
    # High MT removed:       543
    # Retained: 70287 / 72340 cells (97.2%)
    """
    print("[filter] Applying MAD-based filters...")
    n_before = adata.n_obs
    keep     = pd.Series(True, index=adata.obs_names)

    low_counts = _outlier_per_sample(adata, "log1p_total_counts",      nmads_counts, "low",  groupby)
    low_genes  = _outlier_per_sample(adata, "log1p_n_genes_by_counts", nmads_genes,  "low",  groupby)
    hi_counts  = _outlier_per_sample(adata, "log1p_total_counts",      nmads_counts, "high", groupby)
    keep &= ~low_counts & ~low_genes & ~hi_counts
    print(f"  Low counts removed:   {low_counts.sum():>6}")
    print(f"  Low genes removed:    {low_genes.sum():>6}")
    print(f"  High counts removed:  {hi_counts.sum():>6}")

    mt_removed = 0
    for smp in adata.obs[groupby].unique():
        if smp not in decision.index or not decision.loc[smp, "filter_mt"]:
            continue
        mask   = adata.obs[groupby] == smp
        mt     = adata.obs.loc[mask, "pct_counts_mt"]
        _, hi  = _mad_thresholds(mt, nmads_mt)
        flag   = mask & (adata.obs["pct_counts_mt"] > hi)
        keep  &= ~flag
        mt_removed += flag.sum()
    print(f"  High MT removed:      {mt_removed:>6}")

    ribo_removed = 0
    for smp in adata.obs[groupby].unique():
        if smp not in decision.index or not decision.loc[smp, "filter_ribo"]:
            continue
        mask   = adata.obs[groupby] == smp
        ribo   = adata.obs.loc[mask, "pct_counts_ribo"]
        _, hi  = _mad_thresholds(ribo, nmads_ribo)
        flag   = mask & (adata.obs["pct_counts_ribo"] > hi)
        keep  &= ~flag
        ribo_removed += flag.sum()
    print(f"  High ribo removed:    {ribo_removed:>6}")

    dbl_removed = 0
    if "predicted_doublet" in adata.obs.columns:
        for smp in adata.obs[groupby].unique():
            if smp not in decision.index or not decision.loc[smp, "filter_doublets"]:
                continue
            flag       = (adata.obs[groupby] == smp) & adata.obs["predicted_doublet"]
            keep      &= ~flag
            dbl_removed += flag.sum()
    print(f"  Doublets removed:     {dbl_removed:>6}")

    adata_filt = adata[keep].copy()
    print(f"\n  Retained: {adata_filt.n_obs:,} / {n_before:,} cells "
          f"({adata_filt.n_obs / n_before * 100:.1f}%)")
    return adata_filt


def check_filter_impact(adata_before: ad.AnnData,
                         adata_after:  ad.AnnData,
                         groupby:      str   = "sample",
                         max_loss_pct: float = 30.0) -> None:
    """
    Warn if filtering removed more than max_loss_pct of cells per sample.

    Run immediately after apply_mad_filters to catch over-filtering
    that might remove rare cell populations.

    Parameters
    ----------
    adata_before  : AnnData before filtering.
    adata_after   : AnnData after filtering.
    groupby       : per-sample column.
    max_loss_pct  : warning threshold (default 30%).

    Example
    -------
    >>> adata_filt = apply_mad_filters(adata, decision)
    >>> check_filter_impact(adata, adata_filt)
    # All samples within threshold ✓
    # — or —
    # WARNING: GSM5139860 lost 47.2% — consider relaxing nmads thresholds
    """
    before  = adata_before.obs.groupby(groupby).size().rename("before")
    after   = adata_after.obs.groupby(groupby).size().reindex(before.index).fillna(0).rename("after")
    report  = pd.concat([before, after], axis=1)
    report["removed_pct"] = ((report["before"] - report["after"]) / report["before"] * 100).round(1)

    print("[filter_impact] Cells retained per sample:")
    print(report.to_string())

    over = report[report["removed_pct"] > max_loss_pct]
    if over.empty:
        print(f"\n  All samples within {max_loss_pct}% threshold ✓")
    else:
        print(f"\n  WARNING: {len(over)} sample(s) lost > {max_loss_pct}%:")
        for sid, row in over.iterrows():
            print(f"    {sid}: {row['removed_pct']:.1f}% removed "
                  f"({int(row['before'])} → {int(row['after'])} cells)")
        print("  → Relax nmads in apply_mad_filters() or inspect those samples")


def filter_genes(adata: ad.AnnData, min_cells: int = 10) -> ad.AnnData:
    """
    Remove genes detected in fewer than min_cells cells.

    Run after apply_mad_filters, before export.

    Parameters
    ----------
    adata     : AnnData.
    min_cells : minimum expressing cells (default 10).

    Returns
    -------
    AnnData  (scanpy filters in-place and returns same object)

    Example
    -------
    >>> adata_filt = filter_genes(adata_filt, min_cells=10)
    # [genes] Removed 4821 genes (< 10 cells). Remaining: 28717
    """
    n_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"[genes] Removed {n_before - adata.n_vars} genes "
          f"(< {min_cells} cells). Remaining: {adata.n_vars}")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION + EXPORT
# ─────────────────────────────────────────────────────────────────────────────

def validate_sparse(adata: ad.AnnData) -> None:
    """
    Assert adata.X is CSR sparse and report storage stats.
    Converts any dense layers to CSR in-place.

    Raises
    ------
    AssertionError if adata.X is not sparse CSR.

    Example
    -------
    >>> validate_sparse(adata_filt)
    # Format: CSR | Shape: (70287, 28717) | Sparsity: 90.7%
    # All checks passed ✓
    """
    print("[sparse] Validating...")
    assert sparse.issparse(adata.X), "adata.X is not sparse!"
    assert adata.X.format == "csr",  f"Expected CSR, got {adata.X.format}"

    sparsity = 1 - adata.X.nnz / (adata.n_obs * adata.n_vars)
    print(f"  Format:    {adata.X.format.upper()}")
    print(f"  Shape:     {adata.X.shape}")
    print(f"  Non-zeros: {adata.X.nnz:,}")
    print(f"  Sparsity:  {sparsity:.1%}")
    print(f"  RAM (data array): {adata.X.data.nbytes / 1e6:.1f} MB")

    for layer in adata.layers:
        lyr = adata.layers[layer]
        if sparse.issparse(lyr):
            print(f"  Layer '{layer}': {lyr.format.upper()} ✓")
        else:
            print(f"  Layer '{layer}': dense → converting to CSR")
            adata.layers[layer] = sparse.csr_matrix(lyr)

    print("[sparse] All checks passed")


def export_h5ad(adata: ad.AnnData,
                out_dir: str,
                prefix:  str = "qc_filtered") -> str:
    """
    Export filtered AnnData as gzip-compressed h5ad.
    Also writes a per-sample QC summary TSV.

    Parameters
    ----------
    adata   : Filtered AnnData.
    out_dir : Output directory (created if absent).
    prefix  : Filename prefix (default 'qc_filtered').

    Returns
    -------
    str  Full path to written .h5ad file.

    Output files
    ------------
    <out_dir>/<prefix>_filtered.h5ad
    <out_dir>/<prefix>_qc_summary.tsv

    Example
    -------
    >>> path = export_h5ad(adata_filt, "../results/", prefix="GSE168434_qc")
    # [export] h5ad: ../results/GSE168434_qc_filtered.h5ad  (892.3 MB)
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    adata.X = adata.X.tocsr()

    h5ad_path = out / f"{prefix}_filtered.h5ad"
    adata.write_h5ad(str(h5ad_path), compression="gzip")
    print(f"[export] h5ad: {h5ad_path}  ({h5ad_path.stat().st_size / 1e6:.1f} MB)")

    qc_cols = [c for c in adata.obs.columns
               if any(k in c for k in ["total_counts", "n_genes", "pct_counts"])]
    if qc_cols and "sample" in adata.obs.columns:
        summary  = adata.obs.groupby("sample")[qc_cols].median().round(2)
        tsv_path = out / f"{prefix}_qc_summary.tsv"
        summary.to_csv(str(tsv_path), sep="\t")
        print(f"[export] QC summary: {tsv_path}")

    return str(h5ad_path)


# ─────────────────────────────────────────────────────────────────────────────
# FULL PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(
    tar_path:  str | None = None,
    data_dir:  str | None = None,
    out_dir:   str        = "./scrna_output",
    species:   str        = "human",
    min_cells: int        = 10,
    groupby:   str        = "sample",
    prefix:    str        = "qc_filtered",
) -> ad.AnnData:
    """
    End-to-end pipeline: extract → load → QC → filter → export.

    Parameters
    ----------
    tar_path  : Path to .tar archive (mutually exclusive with data_dir).
    data_dir  : Path to already-extracted directory.
    out_dir   : Output directory.
    species   : 'human' or 'mouse'.
    min_cells : Minimum cells per gene (default 10).
    groupby   : per-sample column (default 'sample').
    prefix    : Output filename prefix.

    Returns
    -------
    AnnData  Filtered object (also written to disk).

    Example
    -------
    >>> adata = run_pipeline(
    ...     tar_path = "../data/GSE168434_RAW.tar",
    ...     species  = "human",
    ...     out_dir  = "../results/",
    ...     prefix   = "GSE168434_qc",
    ... )
    """
    from utils.data_exploration import extract_tar, load_all_samples

    if species not in ("human", "mouse"):
        raise ValueError(f"species must be 'human' or 'mouse', got '{species}'")
    if tar_path is None and data_dir is None:
        raise ValueError("Provide tar_path or data_dir")

    print("=" * 70)
    print("scRNA-seq QC PIPELINE") 
    print("=" * 70)

    if tar_path:
        data_dir = extract_tar(tar_path, str(Path(out_dir) / "raw"))

    print("\n[Step 2] Loading samples")
    adata = load_all_samples(data_dir)

    print("\n[Step 3] Computing QC metrics")
    if species == "mouse":
        compute_qc_metrics_mouse(adata)
    else:
        compute_qc_metrics_human(adata)

    print("\n[Step 4] Building filter decision table")
    decision, *_ = build_decision_table(adata, groupby=groupby)

    print("\n[Step 5] Filtering cells")
    adata = apply_mad_filters(adata, decision, groupby=groupby)

    print("\n[Step 6] Filtering genes")
    adata = filter_genes(adata, min_cells=min_cells)

    print("\n[Step 7] Validating sparse format")
    validate_sparse(adata)

    print("\n[Step 8] Exporting")
    out_path = export_h5ad(adata, out_dir, prefix=prefix)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print(f"  Cells:  {adata.n_obs:,}")
    print(f"  Genes:  {adata.n_vars:,}")
    print(f"  Output: {out_path}")
    print("=" * 70)
    return adata
