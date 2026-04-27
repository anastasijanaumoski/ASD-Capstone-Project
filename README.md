# ASD Capstone Project  
BrainSpan Gene Co-expression Network Analysis

## Overview
This project analyzes developmental brain gene expression data from the BrainSpan dataset to understand how gene networks are organized across brain regions and how they relate to autism spectrum disorder (ASD).

The goal of this project is to identify highly connected and biologically meaningful genes (network hubs), determine whether these patterns are non-random, and evaluate whether conserved hub genes overlap with known ASD-associated genes from the SFARI database.

## Methods

### Data Processing
- BrainSpan RNA-seq gene expression data was preprocessed and log-transformed
- Genes were filtered based on variability using a z-score threshold (z ≥ 2)
- Expression matrices were constructed for each brain region

### Network Construction
- Gene co-expression networks were built using Pearson correlation
- Networks were thresholded at multiple correlation cutoffs (0.7, 0.85, 0.9)
- Undirected graphs were generated for each brain region

### Network Analysis
- Degree centrality (hubs) and betweenness centrality (bridges) were calculated
- Network structure was compared to randomized networks using degree-preserving rewiring
- Clustering coefficient, path length, and assortativity were evaluated

### Sensitivity Analysis
- Multiple hub thresholds (1%, 2.5%, 5%, 10%) were tested
- Stability of hub selection was assessed using:
  - Jaccard similarity
  - Spearman correlation

### ASD Enrichment
- Hub genes were compared to the SFARI gene database
- Enrichment was evaluated using Fisher’s exact test
- High-confidence ASD genes (scores 1–2) were analyzed separately

## Key Findings
- Brain region networks show significantly higher clustering than randomized networks, indicating non-random structure
- Hub genes are stable across threshold choices, supporting the use of a top 5% cutoff
- Several conserved hub genes overlap with known ASD-associated genes
- These results suggest that network organization may play a role in ASD-related biology

## Repository Structure
scripts/
preprocess_brainspan.R
peoplexgene_hv_donor.R
genexgene_hv_donor.R
degree_and_betweenness.R
scaling.R
statistical_analysis_zscore.R
conserved_hubs_asd_overlap.R

## Requirements
- R (≥ 4.0)
- Packages:
  - dplyr
  - tidyr
  - ggplot2
  - igraph
  - purrr
  - poweRlaw

## Notes
Due to size restrictions, raw BrainSpan data is not included in this repository.  
Paths in scripts may need to be adjusted to local file locations.

## Author
Anastasija Naumoski  
M.S. Bioinformatics
