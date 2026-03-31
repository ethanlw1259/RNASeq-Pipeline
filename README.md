# RNA-Seq Analysis Pipeline

Bioinformatics pipeline developed during an internship at The Translational Genomics Research 
Institute (TGen) for large-scale differential expression analysis across hundreds of cases and 
controls (1B+ reads).

## Overview
Automates alignment, mapping, quantification, and differential expression analysis of RNA-Seq 
data on Linux HPC infrastructure, designed for reproducibility and batch processing.

## Components
- `pipeline.sh` — End-to-end shell pipeline for alignment, mapping, and quantification
- `diffExp.R` — Differential expression analysis and result reporting

## Tools & Languages
- Alignment: [e.g. STAR / HISAT2]
- Quantification: [e.g. featureCounts / Salmon]
- Differential expression: [e.g. DESeq2 / edgeR]
- R, Bash | Linux HPC environment
