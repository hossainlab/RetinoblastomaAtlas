# Project Instructions: RetinoblastomaAtlas (R Workflow)

## Core Mandates
- **Language:** R is the primary language for this project due to hardware constraints (32GB RAM, no GPU).
- **Framework:** Use **Seurat v5** for all single-cell analysis tasks.
- **Integration:** Prefer **Harmony** for batch correction across datasets (GSE168434, GSE249995) as it is memory-efficient and fast.
- **Coding Style:** Scripts must be clean, concise, and step-by-step. Minimal comments, focus on semantic code and progress messages.
- **Memory Management:** Always use `gc()` and delete large objects (matrices, lists) once they are merged or saved to RDS. Use sequential loading for large sample sets.

## Directory Structure
- `R_Script/`: Contains modular R scripts (`01_...`, `02_...`).
- `data/processed/`: Stores intermediate `.rds` objects.
- `results/figures/`: All plots must be publication-quality (high DPI, minimal clutter).
- `results/tables/`: All outputs (markers, scores) in CSV format.

## Standard Tools
- **QC:** Use `DoubletFinder` for doublets and MAD-based adaptive filtering.
- **Annotation:** Reference markers from `research_plan.md`.
- **Pathways:** Use `decoupleR` with PROGENy/DoRothEA for pathway activity scoring.
- **Communication:** Use `CellChat` for comparative ligand-receptor analysis.
- **Trajectory:** Use `Slingshot` for R-native trajectory inference.

## Verification
- After each script, verify the existence and size of the output RDS file.
- Ensure UMAPs show proper integration (mixing of datasets where expected) before proceeding to annotation.
