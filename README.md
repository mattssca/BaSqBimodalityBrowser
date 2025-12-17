# BiModal BaSq Explorer

An interactive Shiny app for exploring bimodal gene expression patterns in Basal-Squamous (BaSq) bladder cancer samples.

![Shiny](https://img.shields.io/badge/Shiny-Apps-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)
![License](https://img.shields.io/badge/License-MIT-yellow)

## Overview

This application analyzes gene expression bimodality in BaSq bladder cancer samples from the IMvigor210 clinical trial. The tool identifies genes that exhibit dual expression states within BaSq tumors.

**Beyond Differential Expression**: While traditional differential expression analysis identifies genes with statistically different mean expression between clusters, this approach specifically discovers genes with **bimodal "on/off" switching patterns**. These binary-like expression states may represent distinct regulatory programs or cellular states that are missed by mean-based comparisons, potentially revealing novel therapeutic targets and biomarkers with discrete activation thresholds.

## Features

- 🔍 **Interactive Gene Explorer**: Visualize expression distributions for any gene
- 📊 **Bimodality Analysis**: Pre-computed bimodality metrics for ~15,000 genes  
- 🎯 **Smart Filtering**: Exclude proliferation/immune genes to find novel markers
- 📋 **Results Browser**: Sortable table of all bimodality analysis results
- 📁 **Export Functions**: Download plots and data for publication

## Quick Start

### Prerequisites

```r
# Required R packages
install.packages(c("shiny", "ggplot2", "patchwork", "dplyr", "DT"))
```

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/basq_bimodal_explorer.git
cd basq_bimodal_explorer

# Run the app
R -e "shiny::runApp('app.R')"
```

### Usage

1. **Gene Explorer Tab**: Enter any gene symbol to visualize its expression pattern
2. **Bimodality Results Tab**: Browse pre-computed results, filter by signature type
3. **Click any gene** in the results table to automatically visualize it

## Methodology

### Data Source
- **Dataset**: IMvigor210 clinical trial expression data
- **Samples**: BaSq subtype samples (n=90) classified using LundTaxR
- **Platform**: RNA-seq (log2 transformed)

### Analysis Pipeline

#### 1. Sample Clustering
- K-means clustering (k=2) applied to proliferation scores
- Identifies putative high/low proliferation subgroups within BaSq samples

#### 2. Bimodality Detection
**Kernel Density Estimation**: 
- Sheather-Jones bandwidth selection for optimal density estimation
- Generates smooth probability density function for each gene's expression

**Peak/Valley Identification**:
- Local maxima detection identifies expression peaks
- Local minima detection finds valleys between peaks
- Requires ≥2 peaks for bimodality consideration

**Bimodality Strength Calculation**:
```
Bimodality Score = Valley Depth × Peak Separation
Where:
- Valley Depth = (Min Peak Height - Valley Height) / Min Peak Height
- Peak Separation = |Expression difference between top 2 peaks|
```

#### 3. Cluster-Bimodality Alignment
**Valley-Based Classification**:
- Uses deepest valley between peaks as natural threshold
- Classifies samples as "high mode" (> valley) or "low mode" (≤ valley)

**Cluster Separation Quality**:
```
Separation Quality = |Proportion_Cluster1_High - Proportion_Cluster2_High|
Range: 0-1 (1 = perfect separation, 0 = no separation)
```

**Cluster Explains Bimodality Score**:
```
Explanation Score = Separation Quality × 4 × min(p1, 1-p1, p2, 1-p2)
Where:
- p1 = Proportion of cluster1 in high mode
- p2 = Proportion of cluster2 in high mode
- Factor of 4 normalizes the minimum function
- Capped at 1.0
```

#### 4. Quality Filters
- **Minimum samples**: ≥20 samples with valid expression
- **Zero-inflation filter**: Excludes genes with >30% near-zero expression (< 0.5)
- **Mode balance**: Both high/low modes must have ≥10 samples
- **Minimum bimodality**: Initial bimodality score ≥ 0.1

#### 5. Combined Scoring
```
Combined Score = Bimodality Strength × Cluster Explanation Score
```

This multiplicative approach ensures genes must exhibit both:
- Strong bimodal distribution patterns
- Clear alignment with proliferation-based clusters

### Key Metrics Explained

| Metric | Range | Interpretation |
|--------|--------|----------------|
| **Overall Bimodality** | 0-∞ | Strength of dual-peak distribution (higher = more bimodal) |
| **Cluster Separation Quality** | 0-1 | How distinctly clusters separate expression modes (1 = perfect) |
| **Valley Position** | Gene-specific | Expression threshold between low/high modes |
| **Cluster X in High Mode (%)** | 0-100% | Percentage of cluster samples above valley threshold |
| **Cluster Explains Bimodality** | 0-1 | How well clusters account for bimodal pattern (1 = perfect) |
| **Combined Score** | 0-∞ | Overall gene ranking (bimodality × cluster alignment) |

### Statistical Approach

**Advantages of this method**:
- Non-parametric (no distributional assumptions)
- Data-driven thresholding (valley detection vs. arbitrary cutoffs)  
- Integrates distribution shape with cluster structure
- Robust to outliers through kernel density estimation

**Biological Interpretation**:
- High-scoring genes represent discrete expression states within BaSq tumors
- Cluster alignment suggests proliferation-independent heterogeneity
- Non-blacklisted high scorers are candidate subtype-defining markers

### Validation Controls

**Positive Controls** (Expected high scores):
- Proliferation genes (e.g., MKI67, RFC4, TOP2A)
- Should show cluster1 high, cluster2 low pattern

**Negative Controls** (Expected low scores):
- Housekeeping genes with uniform expression
- Immune infiltration genes (filtered out as blacklisted)

**Novel Discoveries**:
- Non-blacklisted genes with high combined scores
- Represent proliferation-independent bimodal patterns

## File Structure

```
basq_bimodal_explorer/
├── app.R                          # Main Shiny application
├── README.md                      # This file
├── data/
│   └── app_data.RData             # Pre-computed analysis results
├── functions/
│   └── functions.R                # Visualization functions
├── screenshots/                   # App screenshots
└── docs/                          # Additional documentation
```

## Data Objects

The app loads several pre-computed data objects:

- `predicted_imvigor`: Expression matrix and sample classifications
- `basq_proliferations`: Sample clustering and proliferation scores  
- `binomial_results`: Complete bimodality analysis results (~15K genes)
- `proliferation_genes`: Cell cycle/proliferation signature genes
- `infiltration_genes`: Immune infiltration signature genes
- `blacklisted_genes`: Combined proliferation + infiltration genes

## Example Genes

**High Bimodality (Non-proliferation)**:
- `GATA3`: Transcription factor with distinct expression states
- `KRT5`: Epithelial marker showing dual patterns
- `TP63`: Tumor suppressor with bimodal distribution

**Proliferation Controls**:
- `MKI67`: Classic proliferation marker
- `RFC4`: DNA replication complex component
- `TOP2A`: Cell cycle checkpoint gene
