"""
scrna_qc_plots.py
=================
Publication-ready QC and doublet visualizations.
Formatted for high-impact journals (Nature Communications, Nature Methods).

Figures follow journal specifications:
  - Width:  89 mm (single column) or 183 mm (double column)
  - Font:   Arial/Helvetica, 7–8 pt labels, 8–9 pt titles
  - DPI:    300 (TIFF/PNG) + PDF vector copy
  - Colors: colorblind-safe palette (Wong 2011)
  - Style:  minimal axes, no top/right spines

Public functions
----------------
    plot_qc_overview(adata, ...)
        Violin plots of all QC metrics with MAD threshold indicators.
        → Figure 1 panel (per-sample QC distributions)

    plot_qc_scatter(adata, ...)
        Scatter: total counts vs genes, colored by MT%.
        → Shows cell quality landscape + filtering zones.

    plot_doublet_scores(adata, ...)
        Doublet score distributions + simulated vs observed histogram.
        → Doublet detection validation panel.

    plot_filter_summary(adata_before, adata_after, ...)
        Stacked bar: cells removed per category per sample.
        → Filter impact summary panel.

    plot_qc_heatmap(adata, ...)
        Per-sample QC metric heatmap (median values, z-scored).
        → Sample-level QC overview.

    plot_full_qc_panel(adata, adata_filt, decision, ...)
        Combined multi-panel figure (A–F) ready for submission.
        → Drop directly into manuscript.

Example
-------
    from utils.scrna_qc_plots import plot_full_qc_panel

    fig = plot_full_qc_panel(
        adata        = adata_raw,
        adata_filt   = adata_filtered,
        decision     = decision,
        groupby      = "sample",
        out_dir      = "../results/figures/",
        prefix       = "GSE168434_qc",
    )
"""

import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# JOURNAL STYLE CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

# Nature Communications figure widths (mm → inches)
_SINGLE_COL = 89  / 25.4   # 3.50 in
_DOUBLE_COL = 183 / 25.4   # 7.20 in
_DPI        = 300

# Wong (2011) colorblind-safe palette
_WONG = {
    "black":         "#000000",
    "orange":        "#E69F00",
    "sky_blue":      "#56B4E9",
    "green":         "#009E73",
    "yellow":        "#F0E442",
    "blue":          "#0072B2",
    "vermillion":    "#D55E00",
    "pink":          "#CC79A7",
}
_PALETTE   = list(_WONG.values())[1:]   # drop black — use as accent only
_MT_COLOR  = _WONG["vermillion"]
_DBL_COLOR = _WONG["blue"]
_PASS_COLOR = _WONG["green"]
_FAIL_COLOR = _WONG["vermillion"]

# Typography
_FONT_SM  = 6.5    # axis tick labels
_FONT_MD  = 7.5    # axis labels, legends
_FONT_LG  = 8.5    # panel titles
_FONT_PNL = 10     # panel letter (A, B, C...)
_LW       = 0.6    # default line width
_LW_THICK = 1.2    # threshold lines

def _journal_style():
    """Apply Nature-style rcParams globally."""
    mpl.rcParams.update({
        "font.family":          "sans-serif",
        "font.sans-serif":      ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":            _FONT_MD,
        "axes.titlesize":       _FONT_LG,
        "axes.labelsize":       _FONT_MD,
        "xtick.labelsize":      _FONT_SM,
        "ytick.labelsize":      _FONT_SM,
        "legend.fontsize":      _FONT_SM,
        "axes.linewidth":       _LW,
        "xtick.major.width":    _LW,
        "ytick.major.width":    _LW,
        "xtick.major.size":     2.5,
        "ytick.major.size":     2.5,
        "xtick.minor.size":     1.5,
        "ytick.minor.size":     1.5,
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "figure.dpi":           _DPI,
        "savefig.dpi":          _DPI,
        "savefig.bbox":         "tight",
        "savefig.pad_inches":   0.02,
        "pdf.fonttype":         42,    # editable text in Illustrator
        "ps.fonttype":          42,
    })

def _panel_label(ax, letter, x=-0.18, y=1.05):
    """Add bold panel letter (A, B, C...) to axis."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=_FONT_PNL, fontweight="bold",
            va="top", ha="left")

def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def _mad(series):
    med = series.median()
    return (series - med).abs().median()

def _mad_thresholds(series, nmads):
    med     = series.median()
    mad_val = _mad(series)
    return med - nmads * mad_val, med + nmads * mad_val

def _save(fig, out_dir, prefix, suffix):
    """Save figure as PNG + PDF."""
    if out_dir is None:
        return
    p = Path(out_dir)
    p.mkdir(parents=True, exist_ok=True)
    png = p / f"{prefix}_{suffix}.png"
    pdf = p / f"{prefix}_{suffix}.pdf"
    fig.savefig(str(png), dpi=_DPI, bbox_inches="tight")
    fig.savefig(str(pdf), bbox_inches="tight")
    print(f"[plot] Saved: {png}")
    print(f"[plot] Saved: {pdf}")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — QC METRICS OVERVIEW (violin + MAD thresholds)
# ─────────────────────────────────────────────────────────────────────────────

def plot_qc_overview(
    adata,
    groupby:    str   = "sample",
    nmads_counts: int = 5,
    nmads_genes:  int = 5,
    nmads_mt:     int = 3,
    nmads_ribo:   int = 5,
    out_dir:    str   = None,
    prefix:     str   = "qc",
    figsize:    tuple = (_DOUBLE_COL, 3.8),
):
    """
    Violin plots of four key QC metrics across all samples.
    MAD-based thresholds overlaid as dashed red lines with annotations.
    Each sample is a separate violin.

    Panels
    ------
    A: log10(total counts)  — low threshold (empty droplets)
    B: n genes detected     — low threshold (poor quality)
    C: MT %                 — high threshold (dying cells)
    D: Ribo %               — high threshold (stress / contamination)

    Parameters
    ----------
    adata         : AnnData with QC metrics computed.
    groupby       : obs column for per-sample grouping.
    nmads_counts  : MADs for counts threshold.
    nmads_genes   : MADs for gene count threshold.
    nmads_mt      : MADs for MT% threshold.
    nmads_ribo    : MADs for ribo% threshold.
    out_dir       : Directory to save PNG + PDF (None = don't save).
    prefix        : Filename prefix.
    figsize       : Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_qc_overview(adata, groupby="sample", out_dir="../results/figures/")
    """
    _journal_style()

    samples  = adata.obs[groupby].unique()
    n_samps  = len(samples)
    obs      = adata.obs

    metrics = [
        ("log1p_total_counts",      f"log$_{{10}}$(total counts + 1)", nmads_counts, "low",  "A"),
        ("log1p_n_genes_by_counts", "Genes detected (log$_{10}$+1)",   nmads_genes,  "low",  "B"),
        ("pct_counts_mt",           "Mitochondrial reads (%)",          nmads_mt,     "high", "C"),
        ("pct_counts_ribo",         "Ribosomal reads (%)",             nmads_ribo,   "high", "D"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=figsize)
    fig.subplots_adjust(wspace=0.45)

    sample_colors = {s: _PALETTE[i % len(_PALETTE)] for i, s in enumerate(samples)}

    for ax, (metric, ylabel, nmads, direction, letter) in zip(axes, metrics):
        if metric not in obs.columns:
            ax.set_visible(False)
            continue

        # ── Violins ─────────────────────────────────────────────────────
        data_per_sample = [obs.loc[obs[groupby] == s, metric].dropna().values
                           for s in samples]
        vp = ax.violinplot(data_per_sample,
                           positions=range(n_samps),
                           showmedians=True,
                           showextrema=False,
                           widths=0.7)

        for i, (body, smp) in enumerate(zip(vp["bodies"], samples)):
            body.set_facecolor(sample_colors[smp])
            body.set_alpha(0.75)
            body.set_linewidth(_LW)
            body.set_edgecolor("white")

        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(_LW_THICK)

        # ── Global MAD thresholds ────────────────────────────────────────
        all_vals  = obs[metric].dropna()
        lo, hi    = _mad_thresholds(all_vals, nmads)

        threshold = hi if direction == "high" else lo
        thr_label = (f"+{nmads} MAD" if direction == "high"
                     else f"−{nmads} MAD")

        ax.axhline(threshold, color=_MT_COLOR, lw=_LW_THICK,
                   linestyle="--", zorder=5, alpha=0.9)
        ax.text(n_samps - 0.4, threshold,
                thr_label, color=_MT_COLOR,
                fontsize=_FONT_SM - 0.5, va="bottom", ha="right",
                fontweight="bold")

        # Shade filtered region
        ylims = ax.get_ylim()
        if direction == "high":
            ax.axhspan(threshold, ylims[1] * 1.05,
                       alpha=0.06, color=_MT_COLOR, zorder=0)
        else:
            ax.axhspan(ylims[0] * 0.98, threshold,
                       alpha=0.06, color=_MT_COLOR, zorder=0)

        # ── Aesthetics ───────────────────────────────────────────────────
        ax.set_xticks(range(n_samps))
        ax.set_xticklabels(
            [s.split("_")[0] for s in samples],   # show GSM ID only
            rotation=90, ha="center", fontsize=_FONT_SM - 1,
        )
        ax.set_ylabel(ylabel, fontsize=_FONT_MD)
        ax.set_xlabel("Sample", fontsize=_FONT_MD)
        _despine(ax)
        _panel_label(ax, letter)

        # Cell count annotation on top
        n_cells = obs[groupby].value_counts()
        ax.set_title(f"n = {len(obs):,} cells", fontsize=_FONT_SM, pad=4, color="#555555")

    _save(fig, out_dir, prefix, "qc_overview")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — SCATTER: COUNTS vs GENES (colored by MT%)
# ─────────────────────────────────────────────────────────────────────────────

def plot_qc_scatter(
    adata,
    groupby:     str   = "sample",
    nmads_counts: int  = 5,
    nmads_genes:  int  = 5,
    nmads_mt:     int  = 3,
    max_samples:  int  = 9,
    out_dir:      str  = None,
    prefix:       str  = "qc",
    point_size:   float = 0.8,
):
    """
    Per-sample scatter: log(total counts) vs log(genes detected).
    Points colored by MT% (low=blue, high=red).
    MAD-based filter boundaries drawn as dashed lines.
    Filtered cells highlighted in the corner.

    Shows the full cell quality landscape in one view.

    Parameters
    ----------
    adata        : AnnData with QC metrics.
    groupby      : obs column for per-sample grouping.
    nmads_counts : MADs for count thresholds.
    nmads_genes  : MADs for gene count threshold.
    nmads_mt     : MADs for MT% threshold.
    max_samples  : Maximum samples to plot (first N, alphabetical).
    out_dir      : Save directory.
    prefix       : Filename prefix.
    point_size   : Scatter point size.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_qc_scatter(adata, max_samples=6, out_dir="../results/figures/")
    """
    _journal_style()

    samples  = sorted(adata.obs[groupby].unique())[:max_samples]
    n_samps  = len(samples)
    ncols    = min(3, n_samps)
    nrows    = int(np.ceil(n_samps / ncols))
    fw       = _DOUBLE_COL
    fh       = 2.2 * nrows + 0.3

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(fw, fh),
                             squeeze=False)
    fig.subplots_adjust(wspace=0.38, hspace=0.55)

    # Shared colormap: MT%
    mt_max = adata.obs["pct_counts_mt"].quantile(0.99) if "pct_counts_mt" in adata.obs else 30
    cmap   = mpl.cm.RdYlBu_r
    norm   = mpl.colors.Normalize(vmin=0, vmax=mt_max)

    for idx, smp in enumerate(samples):
        ax   = axes[idx // ncols][idx % ncols]
        mask = adata.obs[groupby] == smp
        sub  = adata.obs[mask]

        x = sub["log1p_total_counts"]
        y = sub["log1p_n_genes_by_counts"]
        c = sub["pct_counts_mt"] if "pct_counts_mt" in sub.columns else np.zeros(len(sub))

        sc = ax.scatter(x, y, c=c, cmap=cmap, norm=norm,
                        s=point_size, alpha=0.6, linewidths=0, rasterized=True)

        # ── MAD threshold lines ──────────────────────────────────────────
        x_lo, x_hi = _mad_thresholds(x, nmads_counts)
        y_lo, _    = _mad_thresholds(y, nmads_genes)
        mt_vals    = sub["pct_counts_mt"] if "pct_counts_mt" in sub.columns else pd.Series([0])
        _, mt_hi   = _mad_thresholds(mt_vals, nmads_mt)

        for val, orient, label_txt in [
            (x_lo, "v", f"−{nmads_counts}σ"),
            (x_hi, "v", f"+{nmads_counts}σ"),
            (y_lo, "h", f"−{nmads_genes}σ"),
        ]:
            if orient == "v":
                ax.axvline(val, color=_MT_COLOR, lw=_LW, linestyle="--", alpha=0.8)
            else:
                ax.axhline(val, color=_MT_COLOR, lw=_LW, linestyle="--", alpha=0.8)

        # Flag cells that would be removed
        removed = (
            (x < x_lo) | (x > x_hi) | (y < y_lo) |
            (c > mt_hi if "pct_counts_mt" in sub.columns else False)
        )
        n_rm = removed.sum()
        n_tot = len(sub)

        ax.set_title(
            f"{smp.split('_')[0]}\n"
            f"n={n_tot:,}  |  removed={n_rm} ({n_rm/n_tot*100:.1f}%)",
            fontsize=_FONT_SM, pad=3,
        )
        ax.set_xlabel("log$_{10}$(counts+1)", fontsize=_FONT_SM)
        ax.set_ylabel("log$_{10}$(genes+1)", fontsize=_FONT_SM)
        _despine(ax)

        if idx == 0:
            _panel_label(ax, "A")

    # Hide unused subplots
    for idx in range(n_samps, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
    cb.set_label("MT%", fontsize=_FONT_MD)
    cb.ax.tick_params(labelsize=_FONT_SM)
    cb.outline.set_linewidth(_LW)

    # Legend for threshold lines
    legend_line = mpl.lines.Line2D([], [], color=_MT_COLOR, lw=_LW,
                                   linestyle="--", label="MAD threshold")
    fig.legend(handles=[legend_line], fontsize=_FONT_SM,
               loc="lower center", ncol=1,
               bbox_to_anchor=(0.46, -0.02),
               frameon=False)

    _save(fig, out_dir, prefix, "qc_scatter")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — DOUBLET SCORES
# ─────────────────────────────────────────────────────────────────────────────

def plot_doublet_scores(
    adata,
    groupby:     str   = "sample",
    max_samples: int   = 9,
    out_dir:     str   = None,
    prefix:      str   = "qc",
):
    """
    Two panels per figure:

    A: Violin plot of doublet scores per sample.
       Samples colored by predicted doublet rate.
       Threshold line at Scrublet's auto-detected cutoff.

    B: Histogram overlay — singlets vs doublets.
       Shows score separation quality.

    Parameters
    ----------
    adata        : AnnData with 'doublet_score' and 'predicted_doublet' in obs.
    groupby      : obs column for per-sample grouping.
    max_samples  : Max samples plotted in violin panel.
    out_dir      : Save directory.
    prefix       : Filename prefix.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_doublet_scores(adata, out_dir="../results/figures/")
    """
    if "doublet_score" not in adata.obs.columns:
        raise ValueError("Run decide_doublet_filter(adata) first to add doublet_score to obs")

    _journal_style()

    fig, axes = plt.subplots(1, 2, figsize=(_DOUBLE_COL, 3.0),
                             gridspec_kw={"width_ratios": [2, 1]})
    fig.subplots_adjust(wspace=0.4)

    obs     = adata.obs
    samples = sorted(obs[groupby].unique())[:max_samples]
    n_samps = len(samples)

    # ── Panel A: Violin per sample ────────────────────────────────────────
    ax = axes[0]

    # Compute doublet rate per sample for color
    dbl_rates = {
        s: obs.loc[obs[groupby] == s, "predicted_doublet"].mean()
        for s in samples
    }
    rate_norm = mpl.colors.Normalize(vmin=0, vmax=max(dbl_rates.values()) * 1.1)
    rate_cmap = mpl.cm.YlOrRd

    data_per_sample = [obs.loc[obs[groupby] == s, "doublet_score"].values
                       for s in samples]
    vp = ax.violinplot(data_per_sample, positions=range(n_samps),
                       showmedians=True, showextrema=False, widths=0.7)

    for i, (body, smp) in enumerate(zip(vp["bodies"], samples)):
        color = rate_cmap(rate_norm(dbl_rates[smp]))
        body.set_facecolor(color)
        body.set_alpha(0.80)
        body.set_linewidth(_LW)
        body.set_edgecolor("white")

    vp["cmedians"].set_color("black")
    vp["cmedians"].set_linewidth(_LW_THICK)

    # Threshold line (global median threshold)
    threshold_vals = [
        obs.loc[(obs[groupby] == s) & (obs["predicted_doublet"]), "doublet_score"].min()
        for s in samples
        if obs.loc[obs[groupby] == s, "predicted_doublet"].any()
    ]
    if threshold_vals:
        global_thr = np.median(threshold_vals)
        ax.axhline(global_thr, color=_DBL_COLOR, lw=_LW_THICK,
                   linestyle="--", zorder=5, alpha=0.9)
        ax.text(n_samps - 0.3, global_thr + 0.01,
                "Doublet threshold", color=_DBL_COLOR,
                fontsize=_FONT_SM - 0.5, va="bottom", ha="right", fontweight="bold")

    ax.set_xticks(range(n_samps))
    ax.set_xticklabels([s.split("_")[0] for s in samples],
                       rotation=90, fontsize=_FONT_SM - 1)
    ax.set_ylabel("Doublet score", fontsize=_FONT_MD)
    ax.set_xlabel("Sample", fontsize=_FONT_MD)
    ax.set_ylim(bottom=0)
    _despine(ax)
    _panel_label(ax, "A")

    # Colorbar — doublet rate
    sm = mpl.cm.ScalarMappable(norm=rate_norm, cmap=rate_cmap)
    cb = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("Predicted doublet rate", fontsize=_FONT_SM)
    cb.ax.tick_params(labelsize=_FONT_SM)
    cb.outline.set_linewidth(_LW)

    # ── Panel B: Singlet vs Doublet score histogram ───────────────────────
    ax2   = axes[1]
    sing  = obs.loc[~obs["predicted_doublet"], "doublet_score"]
    dbl   = obs.loc[obs["predicted_doublet"],  "doublet_score"]
    bins  = np.linspace(0, obs["doublet_score"].quantile(0.999), 50)

    ax2.hist(sing, bins=bins, density=True, alpha=0.65,
             color=_PASS_COLOR, label="Singlet", linewidth=0)
    ax2.hist(dbl,  bins=bins, density=True, alpha=0.65,
             color=_FAIL_COLOR, label="Doublet", linewidth=0)

    if threshold_vals:
        ax2.axvline(global_thr, color=_DBL_COLOR, lw=_LW_THICK,
                    linestyle="--", alpha=0.9)

    # Score separation annotation
    if len(dbl) > 0:
        sep = dbl.mean() - sing.mean()
        ax2.text(0.97, 0.96,
                 f"Separation\n= {sep:.3f}",
                 transform=ax2.transAxes,
                 fontsize=_FONT_SM, va="top", ha="right",
                 color="#333333")

    leg = ax2.legend(fontsize=_FONT_SM, frameon=False, loc="upper left")
    ax2.set_xlabel("Doublet score", fontsize=_FONT_MD)
    ax2.set_ylabel("Density", fontsize=_FONT_MD)

    n_dbl  = obs["predicted_doublet"].sum()
    n_sing = (~obs["predicted_doublet"]).sum()
    ax2.set_title(
        f"Singlets: {n_sing:,}  |  Doublets: {n_dbl:,}\n"
        f"({n_dbl / len(obs) * 100:.1f}% removed)",
        fontsize=_FONT_SM, pad=4, color="#444444",
    )
    _despine(ax2)
    _panel_label(ax2, "B", x=-0.22)

    _save(fig, out_dir, prefix, "doublet_scores")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — FILTER SUMMARY (stacked bar)
# ─────────────────────────────────────────────────────────────────────────────

def plot_filter_summary(
    adata_before,
    adata_after,
    groupby:  str  = "sample",
    out_dir:  str  = None,
    prefix:   str  = "qc",
    figsize:  tuple = (_DOUBLE_COL, 2.8),
):
    """
    Stacked horizontal bar chart: cells retained and removed per sample.
    Removed cells broken down by filter category.

    Parameters
    ----------
    adata_before : AnnData before filtering (must have filter flag columns).
    adata_after  : AnnData after filtering.
    groupby      : obs column for per-sample grouping.
    out_dir      : Save directory.
    prefix       : Filename prefix.
    figsize      : Figure size.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_filter_summary(adata_raw, adata_filtered, out_dir="../results/figures/")
    """
    _journal_style()

    samples = sorted(adata_before.obs[groupby].unique())
    before  = adata_before.obs.groupby(groupby).size().reindex(samples)
    after   = adata_after.obs.groupby(groupby).size().reindex(samples).fillna(0)
    removed = before - after

    # Compute per-category removed counts where flags exist
    categories = {}
    obs = adata_before.obs

    flag_cols = {
        "Low counts":    obs.get("low_counts",   pd.Series(False, index=obs.index)),
        "High MT%":      obs.get("high_mt",      pd.Series(False, index=obs.index)),
        "Doublets":      obs.get("predicted_doublet", pd.Series(False, index=obs.index)),
    }

    # Build stacked bars
    fig, ax = plt.subplots(figsize=figsize)

    bar_colors  = [_PASS_COLOR, _WONG["sky_blue"], _MT_COLOR, _WONG["orange"], _DBL_COLOR]
    bar_labels  = ["Retained", "Low counts/genes", "High MT%", "High ribo%", "Doublets"]

    # Simple version: retained + total removed (if detailed flags missing)
    retained = after.values
    removed_v = removed.values
    y_pos = np.arange(len(samples))

    # Stacked horizontal bars
    ax.barh(y_pos, retained, height=0.6,
            color=_PASS_COLOR, label="Retained", linewidth=0)
    ax.barh(y_pos, removed_v, height=0.6,
            left=retained, color=_FAIL_COLOR, alpha=0.75,
            label="Removed", linewidth=0)

    # Percentage annotations
    for i, (ret, rem) in enumerate(zip(retained, removed_v)):
        total = ret + rem
        pct   = ret / total * 100 if total > 0 else 0
        ax.text(total + total * 0.01, i,
                f"{pct:.0f}%", va="center",
                fontsize=_FONT_SM, color="#333333")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([s.split("_")[0] for s in samples], fontsize=_FONT_SM)
    ax.set_xlabel("Number of cells", fontsize=_FONT_MD)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    _despine(ax)
    _panel_label(ax, "A")

    # Total cells annotation in title
    tot_before = int(before.sum())
    tot_after  = int(after.sum())
    ax.set_title(
        f"Total: {tot_before:,} → {tot_after:,} cells retained "
        f"({tot_after/tot_before*100:.1f}%)",
        fontsize=_FONT_SM, pad=5, color="#444444",
    )

    ax.legend(fontsize=_FONT_SM, frameon=False,
              loc="lower right", ncol=2)

    _save(fig, out_dir, prefix, "filter_summary")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5 — QC METRIC HEATMAP (per sample, z-scored)
# ─────────────────────────────────────────────────────────────────────────────

def plot_qc_heatmap(
    adata,
    groupby: str  = "sample",
    out_dir: str  = None,
    prefix:  str  = "qc",
    figsize: tuple = (_SINGLE_COL + 0.5, 3.2),
):
    """
    Per-sample QC metric heatmap using median values, z-scored across samples.

    Rows   = QC metrics
    Columns = samples
    Color  = z-score of median across samples (diverging: blue=low, red=high)

    Allows rapid identification of outlier samples across all metrics.

    Parameters
    ----------
    adata   : AnnData with QC metrics.
    groupby : obs column for per-sample grouping.
    out_dir : Save directory.
    prefix  : Filename prefix.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_qc_heatmap(adata, out_dir="../results/figures/")
    """
    _journal_style()

    metrics = {
        "log1p_total_counts":      "log(counts)",
        "log1p_n_genes_by_counts": "log(genes)",
        "pct_counts_mt":           "MT%",
        "pct_counts_ribo":         "Ribo%",
        "pct_counts_hb":           "HB%",
    }
    avail  = {k: v for k, v in metrics.items() if k in adata.obs.columns}
    samples = sorted(adata.obs[groupby].unique())

    # Compute per-sample medians
    medians = adata.obs.groupby(groupby)[list(avail.keys())].median()
    medians = medians.reindex(samples)

    # Z-score across samples per metric
    z = (medians - medians.mean()) / (medians.std() + 1e-9)
    z = z.T.rename(index=avail)

    fig, ax = plt.subplots(figsize=figsize)
    n_rows, n_cols = z.shape

    vmax = max(abs(z.values.min()), abs(z.values.max()), 1.5)
    im = ax.imshow(z.values, cmap="RdBu_r", aspect="auto",
                   vmin=-vmax, vmax=vmax)

    # Grid lines
    ax.set_xticks(np.arange(n_cols + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n_rows + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", size=0)

    # Tick labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([s.split("_")[0] for s in samples],
                       rotation=90, fontsize=_FONT_SM - 0.5, ha="center")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(list(avail.values()), fontsize=_FONT_SM)

    # Cell annotations (z-score values)
    for r in range(n_rows):
        for c in range(n_cols):
            val = z.values[r, c]
            color = "white" if abs(val) > vmax * 0.6 else "black"
            ax.text(c, r, f"{val:.1f}", ha="center", va="center",
                    fontsize=_FONT_SM - 1.5, color=color)

    # Colorbar
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label("Z-score (median)", fontsize=_FONT_SM)
    cb.ax.tick_params(labelsize=_FONT_SM)
    cb.outline.set_linewidth(_LW)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax.set_title("Per-sample QC metric summary", fontsize=_FONT_LG, pad=6)
    _panel_label(ax, "A")

    _save(fig, out_dir, prefix, "qc_heatmap")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 6 — BEFORE / AFTER COMPARISON
# ─────────────────────────────────────────────────────────────────────────────

def plot_before_after(
    adata_before,
    adata_after,
    groupby: str   = "sample",
    out_dir: str   = None,
    prefix:  str   = "qc",
    figsize: tuple = (_DOUBLE_COL, 3.0),
):
    """
    Side-by-side violin plots of key QC metrics before and after filtering.

    Panels
    ------
    A: log(total counts) before vs after
    B: Genes detected before vs after
    C: MT% before vs after
    D: Ribo% before vs after

    Clearly shows filtering effect for Methods section figure.

    Parameters
    ----------
    adata_before : AnnData before filtering.
    adata_after  : AnnData after filtering.
    groupby      : obs column for per-sample grouping.
    out_dir      : Save directory.
    prefix       : Filename prefix.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_before_after(adata_raw, adata_filtered, out_dir="../results/figures/")
    """
    _journal_style()

    metrics = [
        ("log1p_total_counts",      "log$_{10}$(counts+1)", "A"),
        ("log1p_n_genes_by_counts", "Genes detected (log+1)", "B"),
        ("pct_counts_mt",           "MT%",    "C"),
        ("pct_counts_ribo",         "Ribo%",  "D"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=figsize)
    fig.subplots_adjust(wspace=0.48)

    colors = {
        "Before": "#aaaaaa",
        "After":  _PASS_COLOR,
    }

    for ax, (metric, ylabel, letter) in zip(axes, metrics):
        if metric not in adata_before.obs.columns:
            ax.set_visible(False)
            continue

        data = {
            "Before": adata_before.obs[metric].dropna().values,
            "After":  adata_after.obs[metric].dropna().values,
        }

        vp = ax.violinplot(
            list(data.values()),
            positions=[0, 1],
            showmedians=True,
            showextrema=False,
            widths=0.65,
        )

        for body, (label, _) in zip(vp["bodies"], data.items()):
            body.set_facecolor(colors[label])
            body.set_alpha(0.75)
            body.set_linewidth(_LW)
            body.set_edgecolor("white")

        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(_LW_THICK)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Before", "After"], fontsize=_FONT_MD)
        ax.set_ylabel(ylabel, fontsize=_FONT_MD)
        _despine(ax)
        _panel_label(ax, letter)

        # Cell count annotations
        n_b = len(data["Before"])
        n_a = len(data["After"])
        ax.set_title(f"n={n_b:,} → {n_a:,}", fontsize=_FONT_SM, pad=3, color="#555555")

    _save(fig, out_dir, prefix, "before_after")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# COMBINED PUBLICATION FIGURE
# ─────────────────────────────────────────────────────────────────────────────

def plot_full_qc_panel(
    adata,
    adata_filt,
    decision         = None,
    groupby:  str    = "sample",
    nmads_counts: int = 5,
    nmads_genes:  int = 5,
    nmads_mt:     int = 3,
    nmads_ribo:   int = 5,
    out_dir:  str    = None,
    prefix:   str    = "qc",
):
    """
    Full publication-ready QC figure combining all panels.

    Layout (6 panels, 2 rows)
    -------------------------
    Row 1:  A) QC violins (4 metrics)   B) Scatter: counts vs genes
    Row 2:  C) Before/after comparison  D) Doublet scores
            E) Filter summary bar       F) QC heatmap

    Follows Nature Communications figure guidelines:
      - 183 mm wide, Arial 7–8 pt, 300 DPI
      - Colorblind-safe palette (Wong 2011)
      - All threshold lines annotated
      - Panel letters A–F

    Parameters
    ----------
    adata        : AnnData with QC metrics (pre-filter).
    adata_filt   : AnnData after filtering.
    decision     : DataFrame from build_decision_table() — optional.
    groupby      : obs column for per-sample grouping.
    nmads_counts : MADs for count thresholds.
    nmads_genes  : MADs for gene count threshold.
    nmads_mt     : MADs for MT% threshold.
    nmads_ribo   : MADs for ribo% threshold.
    out_dir      : Directory to save PNG + PDF (None = don't save).
    prefix       : Filename prefix.

    Returns
    -------
    matplotlib.figure.Figure

    Example
    -------
    >>> fig = plot_full_qc_panel(
    ...     adata      = adata_raw,
    ...     adata_filt = adata_filtered,
    ...     decision   = decision,
    ...     out_dir    = "../results/figures/",
    ...     prefix     = "GSE168434_qc",
    ... )
    """
    _journal_style()

    fw = _DOUBLE_COL
    fh = 7.5

    fig = plt.figure(figsize=(fw, fh))
    fig.subplots_adjust(hspace=0.62, wspace=0.45)

    gs_top = gridspec.GridSpec(1, 5, figure=fig,
                               top=0.97, bottom=0.57,
                               wspace=0.50)
    gs_bot = gridspec.GridSpec(1, 3, figure=fig,
                               top=0.46, bottom=0.05,
                               wspace=0.45)

    samples  = adata.obs[groupby].unique()
    n_samps  = len(samples)
    obs      = adata.obs

    # ── Metrics config ────────────────────────────────────────────────────
    metric_cfg = [
        ("log1p_total_counts",      f"log$_{{10}}$(counts+1)",  nmads_counts, "low",  "A"),
        ("log1p_n_genes_by_counts", "Genes detected (log+1)",   nmads_genes,  "low",  "B"),
        ("pct_counts_mt",           "MT%",                      nmads_mt,     "high", "C"),
        ("pct_counts_ribo",         "Ribo%",                    nmads_ribo,   "high", "D"),
    ]
    sample_colors = {s: _PALETTE[i % len(_PALETTE)] for i, s in enumerate(samples)}

    # ── Panels A–D: QC violins ────────────────────────────────────────────
    for col, (metric, ylabel, nmads, direction, letter) in enumerate(metric_cfg):
        ax = fig.add_subplot(gs_top[0, col])
        if metric not in obs.columns:
            ax.set_visible(False)
            continue

        data_list = [obs.loc[obs[groupby] == s, metric].dropna().values
                     for s in samples]
        vp = ax.violinplot(data_list, positions=range(n_samps),
                           showmedians=True, showextrema=False, widths=0.7)
        for body, smp in zip(vp["bodies"], samples):
            body.set_facecolor(sample_colors[smp])
            body.set_alpha(0.72)
            body.set_linewidth(_LW)
            body.set_edgecolor("white")
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(_LW_THICK)

        all_vals = obs[metric].dropna()
        lo, hi   = _mad_thresholds(all_vals, nmads)
        threshold = hi if direction == "high" else lo

        ax.axhline(threshold, color=_MT_COLOR, lw=_LW_THICK,
                   linestyle="--", alpha=0.85, zorder=5)
        ax.text(n_samps - 0.3, threshold,
                f"±{nmads} MAD",
                color=_MT_COLOR, fontsize=_FONT_SM - 1,
                va="bottom", ha="right", fontweight="bold")

        ylims = ax.get_ylim()
        if direction == "high":
            ax.axhspan(threshold, ylims[1] * 1.05,
                       alpha=0.06, color=_MT_COLOR, zorder=0)
        else:
            ax.axhspan(max(ylims[0] * 0.98, 0), threshold,
                       alpha=0.06, color=_MT_COLOR, zorder=0)

        ax.set_xticks(range(n_samps))
        ax.set_xticklabels([s.split("_")[0] for s in samples],
                           rotation=90, fontsize=max(_FONT_SM - 2, 4))
        ax.set_ylabel(ylabel, fontsize=_FONT_MD - 0.5)
        ax.set_xlabel("Sample", fontsize=_FONT_SM)
        _despine(ax)
        _panel_label(ax, letter, x=-0.22)

    # ── Panel E: Doublet scores ───────────────────────────────────────────
    ax_e = fig.add_subplot(gs_top[0, 4])

    if "doublet_score" in obs.columns:
        sing  = obs.loc[~obs["predicted_doublet"], "doublet_score"]
        dbl   = obs.loc[obs["predicted_doublet"],  "doublet_score"]
        bins  = np.linspace(0, obs["doublet_score"].quantile(0.999), 45)

        ax_e.hist(sing, bins=bins, density=True, alpha=0.65,
                  color=_PASS_COLOR, label="Singlet", linewidth=0)
        ax_e.hist(dbl,  bins=bins, density=True, alpha=0.65,
                  color=_FAIL_COLOR, label="Doublet", linewidth=0)

        # Threshold line
        if len(dbl) > 0:
            thr_vals = [obs.loc[(obs[groupby] == s) & (obs["predicted_doublet"]),
                                "doublet_score"].min()
                        for s in samples
                        if obs.loc[obs[groupby] == s, "predicted_doublet"].any()]
            if thr_vals:
                thr = np.median(thr_vals)
                ax_e.axvline(thr, color=_DBL_COLOR, lw=_LW_THICK,
                             linestyle="--", alpha=0.9)
                sep = dbl.mean() - sing.mean()
                ax_e.text(0.97, 0.96,
                          f"sep={sep:.3f}",
                          transform=ax_e.transAxes,
                          fontsize=_FONT_SM - 0.5, va="top", ha="right",
                          color="#333333")

        n_dbl = obs["predicted_doublet"].sum()
        ax_e.set_title(f"Doublets: {n_dbl:,}\n({n_dbl/len(obs)*100:.1f}%)",
                       fontsize=_FONT_SM, pad=3, color="#444444")
        ax_e.legend(fontsize=_FONT_SM - 0.5, frameon=False, loc="upper left")
        ax_e.set_xlabel("Doublet score", fontsize=_FONT_MD - 0.5)
        ax_e.set_ylabel("Density", fontsize=_FONT_MD - 0.5)
        _despine(ax_e)
    else:
        ax_e.text(0.5, 0.5, "Run decide_doublet_filter()\nfirst",
                  ha="center", va="center", transform=ax_e.transAxes,
                  fontsize=_FONT_SM, color="#888888")

    _panel_label(ax_e, "E", x=-0.22)

    # ── Panel F: Before/After violin (counts) ─────────────────────────────
    ax_f = fig.add_subplot(gs_bot[0, 0])
    for pos, (adata_set, label, color) in enumerate([
        (adata,      "Before", "#aaaaaa"),
        (adata_filt, "After",  _PASS_COLOR),
    ]):
        vals = adata_set.obs["log1p_total_counts"].dropna().values
        vp   = ax_f.violinplot([vals], positions=[pos],
                               showmedians=True, showextrema=False, widths=0.65)
        vp["bodies"][0].set_facecolor(color)
        vp["bodies"][0].set_alpha(0.75)
        vp["bodies"][0].set_linewidth(_LW)
        vp["bodies"][0].set_edgecolor("white")
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(_LW_THICK)

    ax_f.set_xticks([0, 1])
    ax_f.set_xticklabels(["Before", "After"], fontsize=_FONT_MD)
    ax_f.set_ylabel("log$_{10}$(counts+1)", fontsize=_FONT_MD - 0.5)
    n_b = adata.n_obs
    n_a = adata_filt.n_obs
    ax_f.set_title(f"n: {n_b:,} → {n_a:,}\n({n_a/n_b*100:.1f}% retained)",
                   fontsize=_FONT_SM, pad=3, color="#444444")
    _despine(ax_f)
    _panel_label(ax_f, "F", x=-0.28)

    # ── Panel G: Filter summary bar ───────────────────────────────────────
    ax_g = fig.add_subplot(gs_bot[0, 1])
    before  = adata.obs.groupby(groupby).size()
    after   = adata_filt.obs.groupby(groupby).size().reindex(before.index).fillna(0)
    retained = after.values
    removed  = (before - after).values
    y_pos    = np.arange(len(before))

    ax_g.barh(y_pos, retained, height=0.6, color=_PASS_COLOR,
              label="Retained", linewidth=0)
    ax_g.barh(y_pos, removed, height=0.6, left=retained,
              color=_FAIL_COLOR, alpha=0.7, label="Removed", linewidth=0)

    for i, (ret, rem) in enumerate(zip(retained, removed)):
        total = ret + rem
        pct   = ret / total * 100 if total > 0 else 0
        ax_g.text(total * 1.01, i, f"{pct:.0f}%",
                  va="center", fontsize=max(_FONT_SM - 1.5, 4.5), color="#333333")

    ax_g.set_yticks(y_pos)
    ax_g.set_yticklabels([s.split("_")[0] for s in before.index],
                         fontsize=max(_FONT_SM - 2, 4))
    ax_g.set_xlabel("Cells", fontsize=_FONT_MD - 0.5)
    ax_g.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax_g.legend(fontsize=_FONT_SM - 0.5, frameon=False, loc="lower right")
    _despine(ax_g)
    _panel_label(ax_g, "G", x=-0.18)

    # ── Panel H: QC heatmap ───────────────────────────────────────────────
    ax_h = fig.add_subplot(gs_bot[0, 2])
    hmap_metrics = {
        "log1p_total_counts":      "log(counts)",
        "log1p_n_genes_by_counts": "log(genes)",
        "pct_counts_mt":           "MT%",
        "pct_counts_ribo":         "Ribo%",
    }
    avail  = {k: v for k, v in hmap_metrics.items() if k in obs.columns}
    smp_ord = sorted(obs[groupby].unique())
    meds   = obs.groupby(groupby)[list(avail.keys())].median().reindex(smp_ord)
    z      = ((meds - meds.mean()) / (meds.std() + 1e-9)).T.rename(index=avail)

    vmax = max(abs(z.values.min()), abs(z.values.max()), 1.5)
    im   = ax_h.imshow(z.values, cmap="RdBu_r", aspect="auto",
                       vmin=-vmax, vmax=vmax)
    n_rows, n_cols = z.shape
    ax_h.set_xticks(range(n_cols))
    ax_h.set_xticklabels([s.split("_")[0] for s in smp_ord],
                         rotation=90, fontsize=max(_FONT_SM - 2, 4))
    ax_h.set_yticks(range(n_rows))
    ax_h.set_yticklabels(list(avail.values()), fontsize=_FONT_SM - 0.5)

    for r in range(n_rows):
        for c in range(n_cols):
            val   = z.values[r, c]
            color = "white" if abs(val) > vmax * 0.6 else "black"
            ax_h.text(c, r, f"{val:.1f}", ha="center", va="center",
                      fontsize=max(_FONT_SM - 2, 4), color=color)

    cb = fig.colorbar(im, ax=ax_h, shrink=0.65, pad=0.02)
    cb.set_label("Z-score", fontsize=_FONT_SM)
    cb.ax.tick_params(labelsize=_FONT_SM - 0.5)
    cb.outline.set_linewidth(_LW)
    ax_h.spines[:].set_visible(False)
    _panel_label(ax_h, "H", x=-0.20)

    # ── Figure caption note ───────────────────────────────────────────────
    fig.text(0.5, 0.01,
             "Dashed lines indicate MAD-based filtering thresholds. "
             "Shaded regions indicate removed cells. "
             "Palette: Wong (2011) colorblind-safe.",
             ha="center", fontsize=_FONT_SM - 0.5, color="#666666",
             style="italic")

    _save(fig, out_dir, prefix, "full_qc_panel")
    return fig
