# methods_paper_test

This document outlines the data analysis pipelines used in the paper "_______________" 

DNA replication is a fundamental biological process. Assaying DNA replication timing is a critical thing to do......


Breifly, sequencing data from Repli-seq and two versions of S/G1 experiments were mapped and filtered, then the Repli-seq data was processed using the Repliscan program and the S/G1 data was ratioed and smoothed. 



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
trimomatic/x_____x, fastqc/_________, bowtie2/2.5.1, samtools/1.9, deeptools/3.5.4, bedtools/2.31, kentUtils/1.04



