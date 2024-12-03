# README

This document outlines the data analysis performed in the paper "A comparison of genomic methods to assess DNA replication timing" 

DNA replication is a fundamental biological process of interest to a broad range of research areas, including biochemistry, cellular and developmental biology, and evolutionary and comparative biology. The timing of DNA replication is highly conserved across cell cycle divisions and certain commonalities in the replication timing program are observed from mammals to plants. For example, across many chromosomal genomes the distal end of chromosomes replicate earlier than centromeres, and the bulk of the genes replicate earlier than non-genic regions of the genome.

In the paper "A comparison of genomic methods to assess DNA replication timing" we characterized DNA replication timing in maize root tip tissue using both Repli-seq and DNA copy number, or S/G1, methods. The basic methods are previously published and commonly used in other systems, but we present a novel version of the S/G1 approach. Our improved S/G1 method adds a nascent DNA labeling step that enables bivariate flow cytometric sorting to avoid contamination between S and G1 nuclei. Without cross-contamination the comparison of S-phase and G1-phase DNA copy numbers is more accurate. 

Breifly, sequencing data from Repli-seq and two versions of S/G1 experiments were mapped and filtered, then the Repli-seq data was processed using the Repliscan program and the S/G1 ratio was calculated in non-overlapping 10 kb bins and haar wavelet smoothed. 


## Pipeline and software used
### Pipeline
Step  | File 
--- | --- 
Read pre-processing | 01_read_pre_processing
Map and filter trimmed reads | 02_mapping_and_filtering
High and low coverage droplists | 03_droplist
Repliscan | 04_repliscan
S/G1 ratio | 05_SG1_ratio
Haar wavelet smoothing | 06_haar_wavelet_smoothing

### Software required
trimmomatic v0.39, fastqc v0.12.1, bowtie2 v2.5.1, samtools v1.9, sambamba v1.0.0, DeepTools v3.5.4, bedtools v2.31, kentUtils v1.04

### Other files needed
- genome assembly
- chromosome length file - chromosomes only
- chromosome start/end file - chromosomes only
- 10 kb bin file - chromosomes only


