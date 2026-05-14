"""
utils/
======
Single-cell RNA-seq utility package organized into 3 modules.

    data_exploration  — extraction, loading, data integrity
    qc                — metrics, filtering, export, pipeline
    visualization     — publication-ready figures

─────────────────────────────────────────────────────────────
ABOUT _function NAMING (underscore prefix)
─────────────────────────────────────────────────────────────
_name = PRIVATE — called internally, never by you directly
 name = PUBLIC  — everything below is yours to call

Private functions are implementation details inside each module:
  _load_mtx, _load_h5 ...          called inside load_sample()
  _compute_qc_metrics_core         called inside compute_qc_metrics_human/mouse()
  _mad, _outlier_per_sample        called inside apply_mad_filters()
  _strip_tar_ext                   called inside extract_tar()
  _LOADERS, _SPECIES_PATTERNS      internal lookup dicts
─────────────────────────────────────────────────────────────

Quick import (import everything at once)
----------------------------------------
    from utils import *

Selective import (recommended)
------------------------------
    from utils.data_exploration import inspect_tar, load_all_samples
    from utils.qc               import compute_qc_metrics_human, run_pipeline
    from utils.visualization    import plot_full_qc_panel
"""

from utils.data_exploration import (
    inspect_tar, extract_tar, extract_tar_if_needed, extract_multiple_tars,
    detect_format, discover_samples, load_sample, load_all_samples,
    check_raw_counts, detect_species, convert_ensembl_to_symbols,
    handle_duplicate_genes, inject_metadata, load_if_not_exists, diagnose_sample,
)

from utils.qc import (
    compute_qc_metrics_human, compute_qc_metrics_mouse,
    decide_mt_filter, decide_ribo_filter, decide_doublet_filter, build_decision_table,
    check_sample_sizes, apply_mad_filters, check_filter_impact,
    filter_genes, validate_sparse, export_h5ad, run_pipeline,
)

from utils.visualization import (
    plot_qc_overview, plot_qc_scatter, plot_doublet_scores,
    plot_filter_summary, plot_qc_heatmap, plot_before_after, plot_full_qc_panel,
)

__all__ = [
    "inspect_tar", "extract_tar", "extract_tar_if_needed", "extract_multiple_tars",
    "detect_format", "discover_samples", "load_sample", "load_all_samples",
    "check_raw_counts", "detect_species", "convert_ensembl_to_symbols",
    "handle_duplicate_genes", "inject_metadata", "load_if_not_exists", "diagnose_sample",
    "compute_qc_metrics_human", "compute_qc_metrics_mouse",
    "decide_mt_filter", "decide_ribo_filter", "decide_doublet_filter", "build_decision_table",
    "check_sample_sizes", "apply_mad_filters", "check_filter_impact",
    "filter_genes", "validate_sparse", "export_h5ad", "run_pipeline",
    "plot_qc_overview", "plot_qc_scatter", "plot_doublet_scores",
    "plot_filter_summary", "plot_qc_heatmap", "plot_before_after", "plot_full_qc_panel",
]
