# pHGG scRNA-Seq Analysis

## Overview
This repository contains the code and processed outputs for single-cell RNA-seq (scRNA-seq) analysis of pediatric high-grade glioma (pHGG) cell lines and xenografts. The primary objective is to assess gene expression and cell state differences between radiation-resistant and wildtype cell lines transplanted into mice, with and without trametinib treatment.

**Note:** Raw sequencing data is NOT included and is stored externally.

---

## Project and Data Structure

### Data Organization

#### `data/Annotation_reference/`
- Contains all reference and metadata files needed for analysis
- Sample metadata mapping sample IDs to irradiation status and treatment groups

#### `data/R_Data_Objects/`
- Saved R data files (.RDS format)
- Ordination plots and metadata 

#### `data/Seurat_Outputs/`
- Excel files containing marker genes for each cluster
- Output from Seurat `FindAllMarkers()` analysis
- Gene expression summaries by cluster

#### `plots/`
- Visualization outputs and figures generated during downstream analysis
- Publication-ready and exploratory plots

### Important Note
Raw 10X sequencing data and FASTQ files are NOT included in this repository and are stored externally.

---

## Experimental Design

| Group                                 | Cell Line Status         | Mouse Treatment   | Sample IDs  | N |
|---------------------------------------|-------------------------|-------------------|-------------|---|
| Radiation-Resistant + Trametinib     | Radiation-resistant     | Trametinib        | 927, 929, 934  | 3 |
| Radiation-Resistant + Vehicle        | Radiation-resistant     | Vehicle           | 932           | 1 |
| Wildtype + Trametinib                | Wildtype (no radiation) | Trametinib        | 943, 946, 948, 949, 951 | 5 |
| Wildtype + Vehicle                   | Wildtype (no radiation) | Vehicle           | 940, 941, 942, 945, 947, 950, 952 | 7 |
| **TOTAL**                            |                         |                   | **16 samples** | **16** |

### Metadata Annotations
Each cell receives the following metadata:
- **Radiation:** "Radiation" or "No Radiation"
- **Treatment:** "Trametinib" or "Vehicle"
- **Runs:** Sample ID (927, 929, etc.)
- **seurat_clusters:** Cluster assignment (0-11)
- **Radiation_Resistance:** "Radiation Resistant" or "Radiation Sensitive"
- **Radiation_Trametenib_Resistance:** "Only Radiation Resistant" or "Radiation and Trametenib Resistant"

---

## Analysis Pipeline

### 1. Data Import & Merging
- Load 10X filtered feature-barcode count matrices for all 16 samples
- Create individual Seurat objects
- Merge all objects with appropriate cell IDs

### 2. Quality Control & Filtering
- Remove cells with:
  - High mitochondrial percentage (>25%)
  - Low feature count (<200 genes)
  - High feature count (>7,500 genes)

### 3. Preprocessing & Integration
- **Normalization:** Log normalization with `NormalizeData()`
- **Variable Features:** `FindVariableFeatures()`
- **Scaling:** `ScaleData()` with cell cycle regression (S.Score, G2M.Score)
- **Dimensionality Reduction:** PCA
- **Batch Correction:** Harmony integration by sample Run ID
- **UMAP:** Non-linear dimensionality reduction for visualization

### 4. Clustering
- **K-NN Graph:** `FindNeighbors()` with PCA dims 1:20
- **Leiden Clustering:** `FindClusters()` at resolution 0.6
- Results: 12 distinct clusters (0-11)

### 5. Cell Type Annotation
- **ScType:** Brain cell type classification (oligodendrocytes, astrocytes, microglia, cancer stem cells, etc.)
- **Azimuth:** Motor cortex reference mapping for cell state assignment
- **Marker Gene Expression:** DotPlots for visual confirmation of cell types

### 6. Differential Expression
Two main comparisons:

#### 6a. Radiation Resistance (Clusters 3, 6, 9 vs. 0, 1, 2, 4, 5, 7, 8, 10, 11)
- Identifies genes upregulated in radiation-resistant phenotypes
- GSEA analysis using Hallmark gene sets
- Results: `rad_DE` table with volcano plot

#### 6b. Radiation + Trametinib Resistance (Clusters 6, 9 vs. 3)
- Identifies genes conferring additional trametinib resistance
- GSEA analysis
- Results: `rad_trt_DE` table with volcano plot

### 7. Visualization
- UMAP by treatment, radiation status, and cluster
- Volcano plots showing significantly regulated genes
- DotPlots of cluster-specific marker genes
- Gene Set Enrichment plots (Hallmark pathways)

---

## Key Findings & Cluster Characterization

### Radiation Resistant Clusters: 3, 6, 9
- Enriched in proliferation pathways
- Upregulated cancer stem cell markers (SOX2, EGFR, TCF4, etc.)
- Higher expression of microglia/immune genes (P2RY12, SIGLEC8, etc.)

### Radiation Sensitive Clusters: 0, 1, 2, 4, 5, 7, 8, 10, 11
- Differentiated cell states
- Lower proliferation signatures
- Oligodendrocyte and astrocyte markers

### Additional Trametinib Resistance (Clusters 6, 9 vs. 3)
- Defines a subset of radiation-resistant cells that are also trametinib-resistant
- Specific gene signatures distinct from radiation resistance alone

---

## Dependencies

```r
library(Seurat)           # v4.0+
library(SeuratObject)
library(SeuratDisk)
library(harmony)          # Batch correction
library(ggplot2)          # Plotting
library(patchwork)        # Multi-panel plots
library(clusterProfiler)  # GSEA
library(msigdbr)          # Gene set databases
library(tidyverse)        # Data manipulation
library(ggrepel)          # Repelled labels for plots
library(azimuth)          # Cell type annotation
```

Install all dependencies:
```r
# CRAN packages
install.packages(c("ggplot2", "patchwork", "tidyverse", "ggrepel"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Seurat", "SeuratObject", "SeuratDisk", "clusterProfiler", "msigdbr"))

# Harmony (batch correction)
remotes::install_github("immunogenomics/harmony")

# Azimuth (cell type mapping)
remotes::install_github('satijalab/azimuth')
```

---

## How to Run

1. **Clone the repository**
    ```bash
    git clone https://github.com/jboconnor13/pHGG-scRNA-Analysis.git
    cd pHGG-scRNA-Analysis
    ```

2. **Install R packages** (see Dependencies section above)

3. **Open main analysis script in RStudio**
    - The workflow is structured as an R Markdown file with sequential chunks
    - Each chunk is labeled (Import, Processing, Differential Expression, Annotation)

4. **Update data paths** if needed
    - Paths to external 10X data directories
    - Update `file_dir` in markdown

5. **Run analysis**
    - Execute chunks sequentially or knit the entire document
    - Outputs include cluster assignments, marker genes, DE tables, and visualizations

---

## Output Files

| File | Description |
|------|-------------|
| `data/R_Data_Objects/` | Main Seurat objects including metadata ordination axes |
| `data/Seurat_Outputs/clustermarkers.xlsx` | Marker genes for each cluster (FindAllMarkers output) |
| `plots` | all visualizations |

---

## Contact & Attribution

**Author:** John O'Connor  
**Lab:** Green Lab, CU Anschutz  
**Date Updated:** February 2026

---

## License

This project is licensed under the MIT License - see LICENSE file for details.