
# RNS Gene Expression Analysis

This repository contains scripts used for gene expression processing,
integration, and visualization for the study:

Lee et al. (2026) – Impacts of gene duplication in the evolution of symbiotic root nodule symbiosis in legumes

Authors:
Hyun-oh Lee,  Andrew D. Farmer, Jamie A. O’Rourke, Jeffrey J. Doyle, Steven B. Cannon


The purpose of this repository is to provide a transparent and reproducible
record of the analytical steps used to generate expression summaries,
dimensionality‐reduction analyses, and visualization figures.

---

## Overview

The repository contains three major components:

1. RNA‑seq expression generation for selected species
2. Integration of external expression datasets
3. Downstream statistical analysis and visualization

Species included in the analysis:

- Chamaecrista fasciculata
- Medicago truncatula
- Glycine max
- Phaseolus vulgaris

---

## Repository Structure

scripts/

    expression_generation/
        chafa/
            Salmon_manual_chafa.md
            EdgeR_manual_chafa.md

        medtr/
            Salmon_manual_medtr.md
            EdgeR_manual_medtr.md

    downstream_analysis/
        gene_count_distribution/
            gene_count_distribution.R

        dimensionality_reduction/
            pca_umap_bias_plots.R

        heatmaps/
            heatmap_pdf.R
            heatmap_png.R

---

## Workflow Summary

### RNA‑seq Quantification

For *Chamaecrista* and *Medicago*:

1. RNA‑seq reads quantified using **Salmon**
2. Transcript abundance imported with **tximport**
3. Differential expression analysis performed using **edgeR**
4. Gene‑level expression summaries generated

### External Expression Integration

For *Glycine* and *Phaseolus*:

1. Public expression tables summarized by tissue
2. Relative expression calculated
3. Expression values merged into the gene‑family table

---

## Downstream Analyses

### Gene Count Distribution

Script:
scripts/downstream_analysis/gene_count_distribution/gene_count_distribution.R

Functions:

- Summarizes gene copy numbers across species
- Performs Kruskal–Wallis tests
- Generates boxplot visualizations

---

### PCA and UMAP Visualization

Script:
scripts/downstream_analysis/dimensionality_reduction/pca_umap_bias_plots.R

Functions:

- Principal component analysis (PCA)
- Uniform Manifold Approximation and Projection (UMAP)
- Visualization of expression bias patterns

---

### Heatmap Generation

Scripts:
scripts/downstream_analysis/heatmaps/

Functions:

- Generates relative expression heatmaps
- Species ordering based on phylogenetic relationships
- Visualization using ComplexHeatmap

---

## Software Requirements

R packages used include:

- ggplot2
- dplyr
- tidyr
- readxl
- ComplexHeatmap
- patchwork
- umap

Python scripts require:

- pandas
- numpy

---

## Notes

File paths inside scripts may need to be adjusted depending on the local
directory structure.

These scripts represent the analytical steps used in the manuscript and
are provided for reproducibility purposes.

---

## Citation

If you use these scripts, please cite:

Lee et al. (2026) – Impacts of gene duplication in the evolution of symbiotic root nodule symbiosis in legumes
