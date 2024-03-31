# Thesis: Quality Control and Characterisation of High-Volume Gene Expression Data in Public Data Repositories

## Overview
This repository contains the R scripts developed for my thesis titled "Quality Control and Characterisation of High-Volume Gene Expression Data in Public Data Repositories." The research was conducted at the Agricultural University of Athens and the Biomedical Sciences Research Center "Alexander Fleming" (BSRC Fleming) in Greece. The scripts are entirely written in R.

## Objective
The primary aim of this work is to assess the quality of gene expression data stored in public repositories. We performed random sampling of human (Homo sapiens) and mouse (Mus musculus) data and subjected it to a series of quality control steps. Our goal was to evaluate the overall quality of data within high-volume gene expression data repositories. The analysis provides insights into data quality at both the base level and the sequence alignment level. Our findings may inform the development of bioinformatics tools, such as NGS simulators, and could inspire further research in RNA sequencing data quality.

## Methodology
- **Data Retrieval**: The analysis began with raw data from the GEO repository in SRR format, utilizing the SRA-toolkit for data conversion to FASTQ format.
- **Quality Control**: FASTQ files were analyzed using FastQC. Alignment was performed using HISAT2/BOWTIE2 aligners, resulting in BAM files.
- **Data Analysis**:
  - BAM statistics were gathered using the `bamstats.R` script from P. Moulos Lab at BSRC Fleming ([link](https://github.com/moulos-lab/genomics-facility-processes/blob/main/bamstats.R)). Metrics such as the number of aligned reads and uniquely aligned reads were extracted.
  - Density plots were generated using the `quality_graph.r` script, based on the extracted statistics.
  - The `multiqual.r` tool was developed to analyze the first 50 bases of 100,000 random reads, facilitating the creation of box plots with the `boxplot_homo_mus.r` script.
  - Additional density plots and histograms were derived from FastQC reports using the `quality_graph.r` script.

## Data Files
Scripts utilize large tables for generating graphs, specifically `Total_Table_graphs.xlsx` and `Total_Median_Quality.xlsx`.

## Tags
`#thesis` `#R` `#RNA-seq` `#bioinformatics`
