# WGCNA Co-expression Network Analysis Pipeline

## Overview
A comprehensive and automated pipeline for weighted gene co-expression network analysis (WGCNA) of transcriptomic data.

## Features
- **Complete WGCNA workflow**: From data preprocessing to network construction
- **Advanced filtering**: MAD/Variance-based gene selection  
- **Comprehensive visualization**: Module dendrograms, heatmaps, scatter plots
- **Cytoscape export**: Ready-to-use network files for visualization
- **Parameter tracking**: Automatic logging of analysis parameters
- **Modular output**: Organized results in dedicated subdirectories

## Requirements

### R Packages
- WGCNA
- ggplot2
- reshape2
- ggrepel
- dendextend
- gplots
- circlize
- ComplexHeatmap (Bioconductor)

Install all dependencies using:
```r
source("requirements.R")

Quick Start
1. Install Dependencies
r
source("requirements.R")
2. Run Analysis
r
source("WGCNA_Analysis.R")
3. Expected Output
The analysis will generate results in automatically created timestamped directories including:

Module identification and visualization

Module-trait relationships

Network files for Cytoscape

Comprehensive statistical results

Data Format
Expression Data Format (GENEFPKM.txt)
text
GeneID    Sample1    Sample2    Sample3
Gene_001    10.5      15.2       8.7
Gene_002    25.1      18.4      22.9
Trait Data Format (biaoxing.txt)
text
SampleID    Tissue    Compound content
Sample1     Leaf      25.1
Sample2     Leaf      18.5
Gene Annotation Format (gene_annotation.txt)
text
GeneID      GeneSymbol    Description
Gene_001    ACT1          肌动蛋白
Gene_002    TUB2          微管蛋白
File Structure
text
WGCNA-Pipeline/
├── README.md
├── WGCNA_Analysis.R
├── requirements.R
├── data/
│   ├── GENEFPKM.txt
│   ├── biaoxing.txt
│   └── gene_annotation.txt
└── results/
