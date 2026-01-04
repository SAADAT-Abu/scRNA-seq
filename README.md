# Single-Cell RNA-Seq Analysis
A comprehensive single-cell RNA-seq analysis pipeline replicating the study by Ravikanth Nanduri et al. (2022) examining the role of HMGN proteins in white adipocyte browning using wild-type (WT) versus HMGN1/HMGN2 double knockout (DKO) mice.

## Overview

This analysis investigates adipocyte differentiation and browning mechanisms by comparing transcriptomic profiles between WT and DKO conditions across two developmental time points (Day 0: preadipocytes, Day 6: differentiated adipocytes). The workflow implements quality control with miQC, SCTransform normalization, dimensionality reduction, clustering, marker identification, and cell type annotation.

## Experimental Design

- **Conditions**: Wild-type (WT) vs HMGN1/HMGN2 double knockout (DKO)
- **Time Points**: Day 0 (preadipocytes) and Day 6 (differentiated adipocytes)
- **Replicates**: 2 biological replicates per condition per time point
- **Total Samples**: 8 (2 conditions × 2 time points × 2 replicates)
- **Data Format**: 10X Genomics HDF5 files

## Key Features

- **Quality Control**: Probabilistic identification of compromised cells using miQC
- **Normalization**: SCTransform with glmGamPoi backend for improved speed
- **Batch Correction**: Regression of mitochondrial percentage
- **Cell Type Identification**: 
  - Manual annotation using known adipocyte markers
  - Automated annotation using SingleR with MouseRNAseq reference
- **Differential Expression**: 
  - Single-cell level using Wilcoxon rank-sum test
  - Pseudo-bulk level using DESeq2
- **Visualization**: UMAP, heatmaps, feature plots, and dot plots

## Dependencies

### R Packages
```r
# Core single-cell analysis
Seurat (>= 4.0)
SingleCellExperiment
scater

# Quality control
miQC
scuttle

# Normalization and differential expression
glmGamPoi
presto

# Cell type annotation
SingleR
celldex

# Data manipulation and visualization
tidyverse
```

## Installation

```r
# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "scater", "miQC", 
                       "glmGamPoi", "SingleR", "celldex", "scuttle"))

# Install CRAN packages
install.packages(c("Seurat", "tidyverse", "presto"))
```

## Workflow

### 1. Data Loading
- Read 10X HDF5 files for all 8 samples
- Create individual Seurat objects with basic filtering (min.cells = 3, min.features = 200)
- Merge samples into a single object with preserved sample identity

### 2. Metadata Addition
Structured metadata includes:
- `condition`: WT or DKO
- `time_point`: Day 0 or Day 6
- `cond_tp`: Combined condition and time point

### 3. Quality Control (miQC)
- Identify mitochondrial genes (^mt-)
- Calculate per-cell QC metrics
- Fit mixture model to identify compromised cells
- Filter cells with probability of compromise < 0.75

### 4. Normalization
- SCTransform with:
  - glmGamPoi method for speed
  - v2 flavor for improved regularization
  - Regression of percent.mt
  - Memory conservation enabled

### 5. Dimensionality Reduction
- PCA on top variable features
- Visualize PC heatmaps and elbow plot
- Select first 20 PCs for downstream analysis

### 6. Clustering
- Construct k-nearest neighbor graph (dims = 1:20)
- Leiden clustering (algorithm = 4, resolution = 0.1)
- UMAP visualization

### 7. Marker Identification
- Identify cluster-specific markers using Wilcoxon test
- Compare time points (Day 0 vs Day 6)
- Pseudo-bulk differential expression with DESeq2

### 8. Cell Type Annotation
Manual annotation based on known markers:
- `Mmp3`: Preadipocytes
- `Mki67`: Proliferating cells
- `Fabp4`: Differentiating beige adipocytes
- `Scd1`, `Ucp1`, `Ppargc1a`, `Elovl3`, `Cidea`: Differentiated beige adipocytes

Automated annotation using SingleR with MouseRNAseq reference database.

## Project Structure

```
.
├── Script.qmd                 # Main analysis script
├── README.md                  # This file
├── data/                      # Input HDF5 files (not included)
│   └── *_filtered_feature_bc_matrix.h5
└── outputs/                   # Results and figures (generated)
```

## Usage

1. Place 10X HDF5 files in the `data/` directory
2. Update the working directory path in the setup chunk
3. Render the Quarto document:

```r
quarto::quarto_render("Script.qmd")
```

Or run interactively in RStudio.

## Key Findings

This analysis pipeline enables identification of:
- Cell type composition across conditions and time points
- Differentially expressed genes during adipocyte differentiation
- HMGN-dependent transcriptional changes
- Beige adipocyte marker expression patterns

## Output

The rendered HTML report includes:
- Interactive table of contents
- Embedded visualizations (UMAP, heatmaps, feature plots)
- Differential expression tables
- Quality control metrics

## Citation

If you use this analysis workflow, please cite:

Nanduri, R., et al. (2022). "HMGN proteins in white adipocyte browning."

---

**Keywords**: single-cell RNA-seq, Seurat, adipocyte differentiation, SCTransform, miQC, SingleR, pseudobulk analysis
