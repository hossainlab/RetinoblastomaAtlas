# Research Plan: A Harmonized Single-Cell Atlas of Retinoblastoma Local Extension

**Project:** RetinoblastomaAtlas
**Data sources:** GSE168434 (~70,000 cells), GSE249995 (~30,000 cells)

> **Attribution:** All cited articles were retrieved from PubMed
> ([pubmed.ncbi.nlm.nih.gov](https://pubmed.ncbi.nlm.nih.gov)). DOI links are
> provided for every reference.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Background and Rationale](#2-background-and-rationale)
   - 2.1 Retinoblastoma — Epidemiology and Clinical Staging
   - 2.2 Cone Precursor Cell of Origin
   - 2.3 Intraocular vs. Extraocular Extension
   - 2.4 Tumor Microenvironment in RB
   - 2.5 The Case for a Harmonized scRNA-seq Atlas
   - 2.6 Computational Framework
3. [Specific Aims](#3-specific-aims)
4. [Methodology](#4-methodology)
   - Aim 1 — Quality Control, Normalization, Batch Correction
   - Aim 2 — Cell Type Annotation and Tumor Subtype Classification
   - Aim 3 — Cone Precursor Trajectory and RNA Velocity
   - Aim 4 — Cell-Cell Communication and TGF-β Rewiring
5. [Expected Outcomes and Significance](#5-expected-outcomes-and-significance)
6. [Timeline](#6-timeline)
7. [References](#7-references)

---

## 1. Executive Summary

Retinoblastoma (RB) is the most common primary intraocular malignancy of
childhood. Although most cases are curable when confined to the eye, extraocular
extension dramatically increases mortality. The molecular basis for the transition
from intraocular to extraocular disease remains incompletely understood, partly
because existing studies are limited by small sample sizes and lack of integrated
multi-dataset analyses. This project integrates two publicly available
single-cell RNA-seq (scRNA-seq) datasets comprising ~100,000 cells from 11 RB
patients spanning intraocular and extraocular phenotypes. Through a harmonized
analytical pipeline — batch correction via scVI, cell-type annotation, RNA
velocity with scVelo/CellRank, and ligand-receptor inference with CellChat — we
will produce a reference atlas, identify cone precursor evolutionary trajectories
toward invasive states, and decode how the tumor microenvironment (TME) rewires
TGF-β and other invasion-promoting signals during local extension.

---

## 2. Background and Rationale

### 2.1 Retinoblastoma — Epidemiology and Clinical Staging

Retinoblastoma is the prototypic hereditary cancer and the first human tumor
suppressor gene (*RB1*) was discovered in its context. It occurs at an incidence
of ~1 in 15,000–20,000 live births and almost exclusively affects children
under 5 years of age [[1]](#ref1). The presenting signs (leukocoria, strabismus)
are well characterized and the International Classification of Retinoblastoma
(ICRB, groups A–E) and the TNM staging system predict the likelihood of
high-risk pathologic features — including choroidal invasion, retrolaminar optic
nerve involvement, and anterior chamber seeding — that mandate adjuvant
chemotherapy [[2]](#ref2), [[3]](#ref3).

Most RB cases are curable, with survival exceeding 95% in high-income settings.
However, children presenting with extraocular disease (optic nerve transaction,
orbital invasion, or distant metastasis) face dramatically worse outcomes,
largely because RB cells that escape the eye can seed the CSF and disseminate
to the CNS or bone marrow [[4]](#ref4). High-risk pathologic features detected
at enucleation — retrolaminar and/or massive choroidal invasion — are the
strongest predictors of subsequent relapse, yet the single-cell mechanisms that
allow RB cells to breech ocular barriers remain unclear.

### 2.2 Cone Precursor Cell of Origin

Seminal work by Cobrinik and colleagues established that human RB originates from
maturing cone photoreceptor precursors following biallelic *RB1* inactivation.
They demonstrated that developmental stage-specific differences in MDM2 and MYCN
expression sensitize human (but not murine) maturing ARR3⁺ cone precursors to
RB1-loss-driven proliferation, explaining why mouse models fail to recapitulate
human RB [[5]](#ref5). Subsequently, they showed that MYCN-amplified RB
(representing ~2% of cases) originates from both immature (ARR3⁻) and maturing
(ARR3⁺) cone precursors, the latter undergoing rapid dedifferentiation with
loss of ARR3, L/M-opsin, and gain of progenitor markers (SOX9, RXRγ co-expressed
with NRL/BRN3) [[6]](#ref6). A very recent eLife study from the same group used
full-length scRNA-seq to map nascent, immature, and maturing photoreceptor states
and implicate LNCOC1/LNCOC2 lncRNA expression in cone-precursor-specific
proliferative responses to RB1 loss [[7]](#ref7).

Complementary multi-omics work at Institut Curie identified two RB molecular
subtypes [[8]](#ref8): Subtype 1 (earlier onset, mostly heritable, few genetic
alterations, maintains mature cone marker expression) and Subtype 2 (MYCN
amplification, dedifferentiated cone state co-expressing neuronal/ganglion cell
markers, stemness features including low immune/interferon response, E2F and
MYC activation, higher metastatic propensity). This molecular taxonomy directly
informs the present project's cell type annotation strategy. An hESC-based RB
organoid model confirmed the cone-precursor origin and provided a tractable
in vitro system [[9]](#ref9). Most recently, Gudenas et al. (2026) identified a
tumor-associated photoreceptor signature (TAPS) common to RB, pineoblastoma, and
Group 3 medulloblastoma, expanding the potential scope of cone-precursor-directed
therapy [[10]](#ref10).

At the scRNA-seq level, Yang et al. (2021) profiled 14,739 cells from two RB
samples and identified cone precursors as cells of origin, with UBE2C as a
candidate "switch gene" marking transition to malignant progression [[11]](#ref11).
Xu et al. (2022) analyzed publicly available scRNA-seq data and found that cone
precursors cycle abnormally: G1-phase cone precursors predominate early in disease,
shifting toward G2 phase with progression, while microglia exhibit HLA-family
immune response genes [[12]](#ref12).

### 2.3 Intraocular vs. Extraocular Extension

The mechanistic basis for local extension from intraocular to extraocular RB is
poorly characterized at single-cell resolution. Liu et al. (2024) performed the
most comprehensive scRNA-seq analysis of local extension to date: integrating
128,454 cells from 16 samples (their own plus public data) across intraocular
and extraocular phenotypes, they resolved 10 cone precursor-like subpopulations,
8 retinoma-like subpopulations, and 7 MKI67PhrD (proliferating photoreceptorness-
decreased) subpopulations [[13]](#ref13). SOX4 was identified as highly expressed
in extraocular samples — particularly in MKI67PhrD cells — and validated in
additional clinical samples, nominating SOX4 as a driver of local extension.
Copy number variation inference showed chromosome 6p amplification specific to
extraocular samples. Importantly, these data were generated as a prequel to the
present project's datasets and are directly compatible.

Zhang et al. (2025) characterized TME heterogeneity across non-invasive versus
invasive RB using integrated scRNA-seq, revealing that: (i) TAM subcluster MG1
and astrocyte-like subcluster AC1 are specifically enriched in extraocular tumors;
(ii) pseudotime trajectory analysis suggests invasion correlates with accumulation
of immature TAMs and depletion of terminally differentiated astrocyte-like cells;
and (iii) cell-cell communication networks are substantially rewired between
non-invasive and invasive RB [[14]](#ref14).

In the most recent study, Wan et al. (2025) analyzed 10 RB patients and showed
that cone precursor subcluster CP4 shows elevated TGF-β signaling specifically
in invasive RB, while fibroblast-CP interactions are increased in invasive tumors,
and DOK7 was identified as a key invasion-associated gene [[15]](#ref15). These
findings directly motivate our systematic TGF-β signaling axis investigation.

### 2.4 Tumor Microenvironment in RB

The RB TME is predominantly immunosuppressive. Cuadrado-Vilanova et al. (2022)
characterized 23 enucleated RB specimens and showed that CD163⁺ M2-polarized
TAMs are profuse, CD8⁺ TILs are only moderately abundant, and PD-L1 expression
is absent in ~95% of samples — concordant with FOXP3⁺ Treg infiltration
contributing to immune evasion [[16]](#ref16). They identified two key
immunosuppressive secreted factors — EMMPRIN and macrophage migration inhibitory
factor (MIF) — produced at high levels by RB primary cells and quantifiable in
aqueous humor as liquid biopsy biomarkers. Macrophage polarization toward M2
phenotype was directly induced by RB-conditioned medium or recombinant MIF,
establishing a mechanistic link between tumor secretome and microenvironmental
immunosuppression.

With respect to TGF-β specifically: although no scRNA-seq study has yet explicitly
mapped TGF-β ligand-receptor axes across the full RB TME at single-cell
resolution, Rb depletion has been shown to induce EMT-like transcriptional changes
(including E-cadherin downregulation, Slug and ZEB1 upregulation, and rho-GTPase
dysregulation) via a TGF-β-dependent mechanism [[17]](#ref17). Given that Wan
et al. found TGF-β pathway elevation in the invasive CP4 subcluster, and that
TGF-β is a central inducer of EMT and invasion in many solid tumors, systematic
characterization of TGF-β signaling within the RB TME is a high-priority aim.

### 2.5 The Case for a Harmonized scRNA-seq Atlas

Each prior scRNA-seq study of RB has been limited in one or more ways:
- **Small sample size** (≤5 patients in most single-institution studies)
- **Single phenotype** (either intraocular or extraocular, not both)
- **Shallow sequencing or 3′-only protocols** limiting splice-site information
  needed for RNA velocity
- **No systematic batch correction** across datasets
- **Incomplete TME profiling** (limited macrophage/glial subclustering)

The two GEO datasets integrated here (GSE168434: 7 patients, ~70,000 cells;
GSE249995: 4 samples including S1/S2 intraocular and S3/S4 extraocular) together
provide the largest harmonized RB scRNA-seq resource with matched intraocular
vs. extraocular representation. Harmonizing these datasets with rigorous batch
correction will enable statistical power unavailable to any single study.

### 2.6 Computational Framework

**Batch correction — scVI:** The scVI framework (single-cell variational
inference) uses a deep generative model (variational autoencoder) to jointly
learn a latent representation of gene expression that accounts for batch effects
while preserving biological variation [[18]](#ref18). Unlike linear methods
(ComBat, Harmony), scVI explicitly models the negative binomial distribution of
scRNA-seq counts and handles missing data, making it well suited for integrating
two datasets generated with different protocols and sequencing depths.

**Scalable single-cell analysis — Scanpy:** All preprocessing, clustering, and
visualization will use Scanpy [[19]](#ref19), a Python framework that handles
datasets of >1 million cells efficiently through graph-based data structures
and GPU acceleration.

**RNA velocity — scVelo:** Bergen et al. developed a generalized dynamical model
that solves the full transcriptional splicing kinetics (transcription, splicing,
degradation rates) per gene without requiring steady-state assumptions
[[20]](#ref20). Applied to RB cone precursors, scVelo will reveal which cells are
actively transitioning toward invasive states and identify putative driver genes.

**Trajectory and fate mapping — CellRank:** Lange et al. combined RNA velocity
directional information with Markov-chain-based trajectory inference to compute
cell fate probabilities in a statistically principled way [[21]](#ref21). CellRank
is particularly suited to disease contexts (as opposed to purely normal
development) because it does not require knowledge of a fixed trajectory
direction, making it appropriate for RB where cells may simultaneously progress
along invasive and quiescent trajectories.

**Cell-cell communication — CellChat:** CellChat infers and analyzes intercellular
signaling networks from scRNA-seq data using a curated ligand-receptor database
that includes multi-subunit molecular complexes and cofactor modulation
[[22]](#ref22). CellChat v2 adds expanded L-R pairs and comparative analysis
across conditions, enabling direct comparison of communication topologies between
intraocular (S1/S2) and extraocular (S3/S4) tumors [[23]](#ref23). Benchmarking
studies show CellChat among the top-performing tools when evaluated against
spatial transcriptomic ground truth [[24]](#ref24).

---

## 3. Specific Aims

### Aim 1 — Build a harmonized, batch-corrected RB single-cell atlas

**Hypothesis:** scVI-based integration of two independent 10x Genomics datasets
(GSE168434 + GSE249995) will resolve the major cell types and tumor cell
subpopulations present in retinoblastoma while removing technical batch effects.

**Deliverable:** `data/processed/atlas_scvi_integrated.h5ad` — a fully
annotated AnnData object with cell-type labels, UMAP coordinates, and quality
metrics.

---

### Aim 2 — Define the cellular ecosystem of intraocular vs. extraocular RB

**Hypothesis:** Cone precursor-like cells will show quantitatively distinct
subcluster proportions, copy number variation profiles, and cell cycle states
between intraocular and extraocular samples, mirroring the two RB subtypes
identified by Liu et al. (2021) and the SOX4/MKI67PhrD axis identified by Liu
et al. (2024).

**Deliverable:** High-resolution cell type taxonomy with ≥3 cone precursor
subclusters, ≥2 TAM subclusters, and ≥2 glial subclusters, compared across
disease stage.

---

### Aim 3 — Map the cone precursor invasion trajectory using RNA velocity

**Hypothesis:** Cone precursor-like cells undergoing transition toward the
invasive/extraocular state will show coherent velocity vectors and gene expression
kinetics (elevated SOX4, MYCN, cell cycle activators; reduced cone maturation
markers ARR3, L/M-opsin) consistent with active dedifferentiation-driven invasion.

**Deliverable:** scVelo velocity fields on cone precursor subclusters;
CellRank-derived fate probability maps; a list of ≥20 high-confidence transition
driver genes.

---

### Aim 4 — Decode TGF-β and invasion-promoting cell-cell communication

**Hypothesis:** The TME of extraocular RB is characterized by elevated TGF-β
signaling from fibroblasts and/or activated macrophages to cone precursor-like
cells, and by a reduction in interferon/immune surveillance signals, consistent
with an immunosuppressive, invasion-permissive niche.

**Deliverable:** CellChat differential interaction maps (intraocular vs.
extraocular); TGF-β pathway activity scored per cell type using decoupleR;
a mechanistic schematic of the invasion-promoting communication network.

---

## 4. Methodology

### Aim 1 — Quality Control, Normalization, and Batch Correction

#### 4.1.1 Quality Control

Per-cell QC will be performed using Scanpy's `sc.pp.calculate_qc_metrics()`.
Filtering thresholds (applied independently per sample to avoid batch-confounded
filtering):

| Metric | Lower bound | Upper bound |
|---|---|---|
| Genes detected (nFeature_RNA) | 300 | 99th percentile per sample |
| Total UMI counts (nCount_RNA) | 500 | 99th percentile per sample |
| Mitochondrial % | — | 20% |
| Log₁₀(nCount) / log₁₀(nGene) ratio | 1.1 | 1.8 (doublet proxy) |

Doublets will be detected independently per sample with Scrublet (threshold
0.25). Ambient RNA contamination will be corrected with SoupX using
sample-specific empty droplets from each raw 10x matrix.

#### 4.1.2 Normalization and Feature Selection

After merging filtered cells via `anndata.concat(join='outer')` (already
completed in the preprocessing notebook):

1. Library-size normalization to 10,000 counts per cell (`sc.pp.normalize_total`)
2. Log1p transformation (`sc.pp.log1p`)
3. Selection of 3,000 highly variable genes (HVGs) using `sc.pp.highly_variable_genes`
   with `flavor='seurat_v3'`, computed jointly across all samples but with
   per-sample variance estimation to avoid batch-dominated HVG selection
4. Regress out confounders: total counts and mitochondrial fraction

Raw counts will be preserved in `.layers['counts']` for scVI input.

#### 4.1.3 Batch Correction with scVI

scVI [[18]](#ref18) will be run with raw counts as input:

```python
import scvi

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="dataset",           # GSE168434 vs. GSE249995
    categorical_covariate_keys=["sample_id"],
    continuous_covariate_keys=["pct_counts_mt"]
)

model = scvi.model.SCVI(
    adata,
    n_latent=30,
    n_layers=2,
    gene_likelihood="nb",          # negative binomial
    dispersion="gene-batch"
)
model.train(max_epochs=400, early_stopping=True)
```

The 30-dimensional latent space will be used for:
- Batch-corrected neighbor graph (`sc.pp.neighbors(use_rep='X_scVI')`)
- UMAP visualization
- Leiden clustering (resolution 0.5, 1.0, 2.0 for hierarchical annotation)

Integration quality will be assessed with scib-metrics: batch correction
(iLISI, kBET, graph connectivity) vs. biological conservation (cLISI, NMI,
ARI) using an overall integration score.

---

### Aim 2 — Cell Type Annotation and Tumor Subtype Classification

#### 4.2.1 Reference-based Annotation

Cell types will be annotated in two passes:

**Pass 1 — Broad cell types** using marker genes from the Human Retina Cell
Atlas and prior RB scRNA-seq studies [[11]](#ref11), [[13]](#ref13):

| Cell Type | Key Markers |
|---|---|
| Cone precursor-like (CP) | ARR3, RXRG, THRB, PRDM1, CRX |
| Retinoma-like (RL) | p16 (CDKN2A), p130 (RBL2), low Ki67 |
| MKI67-high RB (MKI67PhrD) | MKI67, TOP2A, BIRC5, UBE2C |
| Müller glia / Astrocyte-like | GFAP, VIM, SOX9, RLBP1 |
| Microglia / TAM | IBA1, CD68, AIF1, CX3CR1 |
| Vascular endothelium | PECAM1, CDH5, VWF |
| Pericyte | PDGFRB, ACTA2, RGS5 |
| Fibroblast | DCN, COL1A1, LUM |
| T cell | CD3D, CD8A, CD4 |
| NK / ILC | GNLY, NKG7, KLRD1 |
| B cell / Plasma | CD79A, MS4A1, MZB1 |

**Pass 2 — Fine subclustering** within each major compartment at higher
resolution (leiden resolution 2.0–3.0), guided by literature-defined marker
panels for CP subclusters (CP1–CP10 as in [[13]](#ref13)), TAM subclusters
(M1-like: CXCL10, STAT1; M2-like: CD163, MRC1, MIF), and astrocyte subclusters.

#### 4.2.2 Copy Number Variation Inference

InferCNV will be run using normal retinal cells (Müller glia, endothelial cells)
as reference. Chromosome arm-level CNV profiles will be used to:
- Confirm malignant identity of tumor cells
- Detect chromosome 6p amplification specific to extraocular samples
  (as reported in [[13]](#ref13))
- Identify sample-specific clonal expansions

#### 4.2.3 Tumor Subtype Classification

Each tumor cell will be scored against Subtype 1 (mature cone-like) and
Subtype 2 (dedifferentiated/stemness) gene signatures from Liu et al. (2021)
[[8]](#ref8) using `sc.tl.score_genes()`. A third scoring set will capture the
tumor-associated photoreceptor signature (TAPS) identified by Gudenas et al.
(2026) [[10]](#ref10). Scores will be compared between GSE168434 patients
(baseline intraocular) and GSE249995 intraocular vs. extraocular samples.

---

### Aim 3 — Cone Precursor Invasion Trajectory and RNA Velocity

#### 4.3.1 Spliced/Unspliced Count Extraction

RNA velocity requires unspliced (nascent) and spliced (mature) mRNA counts.
Both GEO datasets were generated with 10x Chromium v2/v3 chemistry; pre-mRNA
alignment will be run using STARsolo with `--soloFeatures Gene Velocyto` on raw
FASTQs, or velocyto CLI on existing BAM files, to produce `.loom` files
per sample.

#### 4.3.2 scVelo Dynamical Model

```python
import scvelo as scv

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
```

Bergen et al.'s dynamical model [[20]](#ref20) will recover gene-specific
transcription (α), splicing (β), and degradation (γ) rates, enabling
identification of genes that are actively upregulated (α > γ: unspliced >
expected) or downregulated (α < γ: unspliced < expected) in each cell.

Analysis will be restricted to cone precursor-like cells and their nearest
neighbours (Müller glia, fibroblasts) to focus velocity interpretation.
Velocity confidence scores will be used to filter cells with low directional
consistency (threshold > 0.5).

#### 4.3.3 CellRank Fate Probability Mapping

```python
import cellrank as cr

vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck

estimator = cr.estimators.GPCCA(combined_kernel)
estimator.compute_macrostates(n_states=6, cluster_key="cell_type")
estimator.set_terminal_states_from_macrostates()
estimator.compute_fate_probabilities()
```

CellRank [[21]](#ref21) will compute the probability of each cone precursor cell
transitioning to: (a) the MKI67PhrD (highly proliferative) state, (b) the
invasive extraocular state (extraocular sample cells), or (c) quiescent
retinoma-like state. Driver gene identification will use CellRank's
`lineage_drivers()` function with Pearson correlation between gene expression
and fate probability.

---

### Aim 4 — Cell-Cell Communication and TGF-β Pathway Rewiring

#### 4.4.1 CellChat Inference

CellChat v2 [[23]](#ref23) will be run independently on:
- **Intraocular cohort:** S1_in1, S2_in2 (GSE249995) + intraocular samples
  from GSE168434 (RB01–RB07 if tumor is intraocular)
- **Extraocular cohort:** S3_ex1, S4_ex2 (GSE249995) + extraocular samples
  from GSE168434

```r
library(CellChat)
cellchat_io  <- createCellChat(adata_io,  group.by = "cell_type_fine")
cellchat_eo  <- createCellChat(adata_eo,  group.by = "cell_type_fine")
# Run both in parallel
for (cc in list(cellchat_io, cellchat_eo)) {
  cc <- setLigandReceptorInteraction(cc, use.mode = "max")
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
}
cellchat_merged <- mergeCellChat(list(cellchat_io, cellchat_eo),
                                  add.names = c("intraocular","extraocular"))
```

**Primary analysis targets:**
- TGF-β pathway: TGFB1/2/3 → TGFBR1/2 signaling between fibroblasts, TAMs,
  and cone precursor-like cells
- MIF pathway: TAM-to-CP MIF/CD44 and MIF/CXCR2 axes
- VEGF pathway: tumor-to-endothelium angiogenic signaling
- SPP1 (osteopontin): M2-TAM communication links

#### 4.4.2 TGF-β Pathway Activity Scoring

TGF-β pathway activity will be scored per cell using decoupleR with the
PROGENy footprint model:

```python
import decoupler as dc

progeny_model = dc.get_progeny(organism='human', top=500)
dc.run_mlm(
    mat=adata,
    net=progeny_model,
    source='source',
    target='target',
    weight='weight',
    obsm_key='progeny'
)
```

TGF-β activity scores will be compared between:
- Intraocular vs. extraocular cone precursor subclusters
- CP subtype 1 vs. subtype 2 cells
- TAM M1-like vs. M2-like subclusters

#### 4.4.3 Spatial Validation Strategy (Planned)

Where tissue availability permits, findings from the scRNA-seq atlas will be
orthogonally validated using published RB immunohistochemistry data (CD163 for
M2-TAMs, SOX4, TGFB1 IHC from Liu et al. [[13]](#ref13) and Wan et al.
[[15]](#ref15)). NicheNet prioritization will be used to identify top-ranked
upstream ligands driving the invasion-associated CP gene expression program,
providing mechanistic hypotheses testable in future wet-lab experiments.

---

## 5. Expected Outcomes and Significance

### Aim 1 Outcomes
- A fully harmonized AnnData object with >95,000 quality-filtered cells
- Demonstrated batch correction (scib overall score >0.7)
- Open-access resource deposited to Zenodo and linked to GEO records

### Aim 2 Outcomes
- First systematic comparison of 10+ RB cell type proportions across intraocular
  and extraocular samples at the single-cell level
- Validation and extension of the Subtype 1/2 taxonomy in a larger patient
  cohort than previously studied
- Identification of chromosome 6p amplification as a clonal marker of
  extraocular cells in the integrated dataset

### Aim 3 Outcomes
- First RNA velocity analysis of cone precursor invasion dynamics in RB
- Identification of putative invasion driver genes (predicted: SOX4, MKI67,
  BIRC5, SH3RF3, and novel candidates)
- CellRank-derived fate probability landscape separating invasive from
  quiescent cone precursor trajectories

### Aim 4 Outcomes
- First systematic quantitative comparison of cell-cell communication topologies
  in intraocular vs. extraocular RB
- Confirmation or refutation of the TGF-β/CP4 axis reported by Wan et al.
  [[15]](#ref15) in a larger, independently-analyzed cohort
- Identification of the dominant immunosuppressive communication axes (MIF,
  SPP1, TGF-β) and their rewiring patterns during local extension

### Broader Significance

This atlas will provide the field's most comprehensive single-cell view of RB
local extension, directly actionable for:
- **Biomarker discovery:** Invasion-specific surface antigens on CP subclusters
  are candidates for ophthalmic liquid biopsy (aqueous humor)
- **Therapeutic target nomination:** TME-rewiring signals (TGF-β, MIF)
  are druggable and may sensitize RB to immune checkpoint blockade
- **Cross-cancer insights:** The tumor-associated photoreceptor signature
  shared with pineoblastoma [[10]](#ref10) suggests shared dependencies that
  could be co-targeted

---

## 6. Timeline

| Phase | Duration | Activities |
|---|---|---|
| **Phase 1** | Months 1–2 | QC, doublet removal, SoupX; scVI batch correction; integration benchmarking (Aim 1) |
| **Phase 2** | Months 3–4 | Broad cell type annotation; InferCNV; tumor subtype scoring (Aim 2) |
| **Phase 3** | Months 5–6 | STAR/velocyto re-alignment for unspliced counts; scVelo dynamical model on CP subset (Aim 3) |
| **Phase 4** | Month 7 | CellRank fate probability mapping; driver gene identification (Aim 3 completion) |
| **Phase 5** | Months 8–9 | CellChat intraocular/extraocular inference; decoupleR TGF-β scoring; NicheNet prioritization (Aim 4) |
| **Phase 6** | Months 10–12 | Integration of all analyses; figure preparation; manuscript writing; data deposition |

**Notebooks roadmap:**

```
notebooks/
├── 01-data-preprocessing.ipynb        ✅ COMPLETE — merged_raw_atlas.h5ad
├── 02-qc-filtering.ipynb              🔲 Aim 1
├── 03-scvi-integration.ipynb          🔲 Aim 1
├── 04-cell-type-annotation.ipynb      🔲 Aim 2
├── 05-infercnv-cna.ipynb              🔲 Aim 2
├── 06-subtype-scoring.ipynb           🔲 Aim 2
├── 07-rna-velocity-scvelo.ipynb       🔲 Aim 3
├── 08-cellrank-trajectories.ipynb     🔲 Aim 3
├── 09-cellchat-communication.ipynb    🔲 Aim 4
└── 10-tgfb-decoupler-scoring.ipynb   🔲 Aim 4
```

---

## 7. References

> All articles below were retrieved from PubMed. DOI links open the primary
> publication record.

<a name="ref1"></a>
**[1]** Balmer A, Zografos L, Munier F. Diagnosis and current management of
retinoblastoma. *Oncogene.* 2006;25(38):5341–9.
[https://doi.org/10.1038/sj.onc.1209622](https://doi.org/10.1038/sj.onc.1209622)
(PMID: 16936756)

<a name="ref2"></a>
**[2]** Yousef YA, Al-Hussaini M, Mehyar M, et al. Predictive value of TNM
classification, International Classification, and Reese-Ellsworth staging of
retinoblastoma for the likelihood of high-risk pathologic features. *Retina.*
2015;35(9):1883–9.
[https://doi.org/10.1097/IAE.0000000000000547](https://doi.org/10.1097/IAE.0000000000000547)
(PMID: 25923953)

<a name="ref3"></a>
**[3]** Pérez V, Sampor C, Rey G, et al. Treatment of nonmetastatic unilateral
retinoblastoma in children. *JAMA Ophthalmol.* 2018;136(7):747–752.
[https://doi.org/10.1001/jamaophthalmol.2018.1501](https://doi.org/10.1001/jamaophthalmol.2018.1501)
(PMID: 29799944)

<a name="ref4"></a>
**[4]** Laurent VE, Sampor C, Solernou V, et al. Detection of minimally
disseminated disease in the cerebrospinal fluid of children with high-risk
retinoblastoma by reverse transcriptase-polymerase chain reaction for GD2
synthase mRNA. *Eur J Cancer.* 2013;49(13):2892–9.
[https://doi.org/10.1016/j.ejca.2013.04.021](https://doi.org/10.1016/j.ejca.2013.04.021)
(PMID: 23721779)

<a name="ref5"></a>
**[5]** Singh HP, Wang S, Stachelek K, et al. Developmental stage-specific
proliferation and retinoblastoma genesis in RB-deficient human but not mouse cone
precursors. *Proc Natl Acad Sci USA.* 2018;115(40):E9391–E9400.
[https://doi.org/10.1073/pnas.1808903115](https://doi.org/10.1073/pnas.1808903115)
(PMID: 30213853)

<a name="ref6"></a>
**[6]** Singh HP, Shayler DWH, Fernandez GE, et al. An immature, dedifferentiated,
and lineage-deconstrained cone precursor origin of N-Myc-initiated retinoblastoma.
*Proc Natl Acad Sci USA.* 2022;119(28):e2200721119.
[https://doi.org/10.1073/pnas.2200721119](https://doi.org/10.1073/pnas.2200721119)
(PMID: 35867756)

<a name="ref7"></a>
**[7]** Shayler DWH, Stachelek K, Cambier L, et al. Identification and
characterization of early human photoreceptor states and cell-state-specific
retinoblastoma-related features. *eLife.* 2025;13:e101918.
[https://doi.org/10.7554/eLife.101918](https://doi.org/10.7554/eLife.101918)
(PMID: 40767624)

<a name="ref8"></a>
**[8]** Liu J, Ottaviani D, Sefta M, et al. A high-risk retinoblastoma subtype
with stemness features, dedifferentiated cone states and neuronal/ganglion cell
gene expression. *Nat Commun.* 2021;12(1):5578.
[https://doi.org/10.1038/s41467-021-25792-0](https://doi.org/10.1038/s41467-021-25792-0)
(PMID: 34552068)

<a name="ref9"></a>
**[9]** Zhang X, Jin ZB. Reconstruct human retinoblastoma in vitro. *J Vis Exp.*
2022;(188).
[https://doi.org/10.3791/62629](https://doi.org/10.3791/62629)
(PMID: 36314812)

<a name="ref10"></a>
**[10]** Gudenas BL, Ahmad ST, Englinger B, et al. A tumor-associated photoreceptor
signature unifies distinct central nervous system malignancies. *Cancer Cell.*
2026;44(4):831–845.e10.
[https://doi.org/10.1016/j.ccell.2026.02.010](https://doi.org/10.1016/j.ccell.2026.02.010)
(PMID: 41791379)

<a name="ref11"></a>
**[11]** Yang J, Li Y, Han Y, et al. Single-cell transcriptome profiling reveals
intratumoural heterogeneity and malignant progression in retinoblastoma. *Cell
Death Dis.* 2021;12(12):1100.
[https://doi.org/10.1038/s41419-021-04390-4](https://doi.org/10.1038/s41419-021-04390-4)
(PMID: 34815392)

<a name="ref12"></a>
**[12]** Xu K, Nie W, Tong Q, et al. Analysis of progress characteristics of
retinoblastoma based on single cell transcriptome sequencing. *Sheng Wu Gong
Cheng Xue Bao.* 2022;38(10):3809–3824.
[https://doi.org/10.13345/j.cjb.220491](https://doi.org/10.13345/j.cjb.220491)
(PMID: 36305411)

<a name="ref13"></a>
**[13]** Liu Y, Hu W, Xie Y, et al. Single-cell transcriptomics enable the
characterization of local extension in retinoblastoma. *Commun Biol.*
2024;7(1):11.
[https://doi.org/10.1038/s42003-023-05732-y](https://doi.org/10.1038/s42003-023-05732-y)
(PMID: 38172218)

<a name="ref14"></a>
**[14]** Zhang X, Liu H, Lan Z, et al. Single-cell transcriptomic analyses provide
insights into the tumor microenvironment heterogeneity and invasion phenotype in
retinoblastoma. *Pathol Res Pract.* 2025;271:156009.
[https://doi.org/10.1016/j.prp.2025.156009](https://doi.org/10.1016/j.prp.2025.156009)
(PMID: 40378583)

<a name="ref15"></a>
**[15]** Wan W, Chen X, Liu H, Yang L, Huang P. Comprehensive analysis of
single-cell and bulk RNA sequencing uncover tumor microenvironment diversity in
invasive retinoblastoma. *Sci Rep.* 2025;15(1):39954.
[https://doi.org/10.1038/s41598-025-23779-1](https://doi.org/10.1038/s41598-025-23779-1)
(PMID: 41238845)

<a name="ref16"></a>
**[16]** Cuadrado-Vilanova M, Liu J, Paco S, et al. Identification of
immunosuppressive factors in retinoblastoma cell secretomes and aqueous humor
from patients. *J Pathol.* 2022;257(3):327–339.
[https://doi.org/10.1002/path.5893](https://doi.org/10.1002/path.5893)
(PMID: 35254670)

<a name="ref17"></a>
**[17]** Arima Y, Inoue Y, Shibata T, et al. Rb depletion results in deregulation
of E-cadherin and induction of cellular phenotypic changes that are characteristic
of the epithelial-to-mesenchymal transition. *Cancer Res.* 2008;68(13):5104–12.
[https://doi.org/10.1158/0008-5472.CAN-07-5680](https://doi.org/10.1158/0008-5472.CAN-07-5680)
(PMID: 18593909)

<a name="ref18"></a>
**[18]** Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling
for single-cell transcriptomics. *Nat Methods.* 2018;15(12):1053–1058.
[https://doi.org/10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2)
(PMID: 30504886)

<a name="ref19"></a>
**[19]** Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene
expression data analysis. *Genome Biol.* 2018;19(1):15.
[https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
(PMID: 29409532)

<a name="ref20"></a>
**[20]** Bergen V, Lange M, Peidli S, Wolf FA, Theis FJ. Generalizing RNA velocity
to transient cell states through dynamical modeling. *Nat Biotechnol.*
2020;38(12):1408–1414.
[https://doi.org/10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3)
(PMID: 32747759)

<a name="ref21"></a>
**[21]** Lange M, Bergen V, Klein M, et al. CellRank for directed single-cell
fate mapping. *Nat Methods.* 2022;19(2):159–170.
[https://doi.org/10.1038/s41592-021-01346-6](https://doi.org/10.1038/s41592-021-01346-6)
(PMID: 35027767)

<a name="ref22"></a>
**[22]** Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of
cell-cell communication using CellChat. *Nat Commun.* 2021;12(1):1088.
[https://doi.org/10.1038/s41467-021-21246-9](https://doi.org/10.1038/s41467-021-21246-9)
(PMID: 33597522)

<a name="ref23"></a>
**[23]** Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell
communication from single-cell transcriptomics. *Nat Protoc.*
2024;20(1):180–219.
[https://doi.org/10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)
(PMID: 39289562)

<a name="ref24"></a>
**[24]** Liu Z, Sun D, Wang C. Evaluation of cell-cell interaction methods by
integrating single-cell RNA sequencing data with spatial information. *Genome
Biol.* 2022;23(1):218.
[https://doi.org/10.1186/s13059-022-02783-y](https://doi.org/10.1186/s13059-022-02783-y)
(PMID: 36253792)

<a name="ref25"></a>
**[25]** Wang Y, Tang J, Liu Y, et al. Targeting ALDOA to modulate tumorigenesis
and energy metabolism in retinoblastoma. *iScience.* 2024;27(9):110725.
[https://doi.org/10.1016/j.isci.2024.110725](https://doi.org/10.1016/j.isci.2024.110725)
(PMID: 39262779)

<a name="ref26"></a>
**[26]** Zhang Z, Tang J, Liu Y, et al. The role of lactate metabolism in
retinoblastoma tumorigenesis and ferroptosis resistance. *Tissue Cell.*
2025;95:102893.
[https://doi.org/10.1016/j.tice.2025.102893](https://doi.org/10.1016/j.tice.2025.102893)
(PMID: 40188688)

<a name="ref27"></a>
**[27]** Zhang MG, Kuznetsoff JN, Owens DA, et al. Early mechanisms of
chemoresistance in retinoblastoma. *Cancers (Basel).* 2022;14(19):4966.
[https://doi.org/10.3390/cancers14194966](https://doi.org/10.3390/cancers14194966)
(PMID: 36230889)

<a name="ref28"></a>
**[28]** Berry JL, Kogachi K, Aziz HA, et al. Risk of metastasis and orbital
recurrence in advanced retinoblastoma eyes treated with systemic chemoreduction
versus primary enucleation. *Pediatr Blood Cancer.* 2017;64(4).
[https://doi.org/10.1002/pbc.26270](https://doi.org/10.1002/pbc.26270)
(PMID: 28221729)

---

*Document generated: 2026-04-29. All PubMed citations retrieved via the
PubMed MCP server.*
