# Pipeline
- Read pre-processing
- Map and filter trimmed reads
- High and low coverage droplists
- Normalize reads
- S/G1 ratio
- Haar Wavelet smoothing

# Software required
trimomatic/x_____x, fastqc/_________, bowtie2/2.5.1, samtools/1.9, deeptools/3.5.4, bedtools/2.31, kentUtils/1.04



# Read pre-processing

## Non-EdU sequencing files
### Get fastq files from NCBI SRA for primary roots
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_primary_G1_br1 | | BR1_Primary_G1_S26_L002_R1_001.fastq.gz | BR1_Primary_G1_S26_L002_R2_001.fastq.gz
B73_primary_S_br1 | | BR1_Primary_S_S27_L002_R1_001.fastq.gz | BR1_Primary_S_S27_L002_R2_001.fastq.gz
B73_primary_G1_br2 | | BR2_Primary_G1_S30_L002_R1_001.fastq.gz | BR2_Primary_G1_S30_L002_R2_001.fastq.gz
B73_primary_S_br2 | | BR2_Primary_S_S31_L002_R1_001.fastq.gz | BR2_Primary_S_S31_L002_R2_001.fastq.gz
B73_primary_G1_br3 | | BR3_Primary_G1_S33_L002_R1_001.fastq.gz | BR3_Primary_G1_S33_L002_R2_001.fastq.gz
B73_primary_S_br23 | | BR3_Primary_S_S34_L002_R1_001.fastq.gz | BR3_Primary_S_S34_L002_R2_001.fastq.gz

### Get fastq files from NCBI SRA for seminal roots
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_seminal_G1_br1 | | BR1_Seminal_G1_S28_L002_R1_001.fastq.gz | BR1_Seminal_G1_S28_L002_R2_001.fastq.gz
B73_seminal_S_br1 | | BR1_Seminal_S_S29_L002_R1_001.fastq.gz | BR1_Seminal_S_S29_L002_R2_001.fastq.gz
B73_seminal_G1_br2 | | BR2_Seminal_G1_S32_L002_R1_001.fastq.gz | BR2_Seminal_G1_S32_L002_R2_001.fastq.gz
B73_seminal_S_br2 | | BR2_Seminal_S_S25_L002_R1_001.fastq.gz | BR2_Seminal_S_S25_L002_R2_001.fastq.gz
B73_seminal_G1_br3 | | BR3_Seminal_G1_S35_L002_R1_001.fastq.gz | BR3_Seminal_G1_S35_L002_R2_001.fastq.gz
B73_seminal_S_br3 | | BR3_Seminal_S_S36_L002_R1_001.fastq.gz | BR3_Seminal_S_S36_L002_R2_001.fastq.gz


### Re-name files
```bash
mv BR1_Primary_G1_S26_L002_R1_001.fastq.gz B73_Primary_G1_br1_R1.fastq.gz
mv BR1_Primary_G1_S26_L002_R2_001.fastq.gz B73_Primary_G1_br1_R2.fastq.gz
mv BR1_Primary_S_S27_L002_R1_001.fastq.gz B73_Primary_S_br1_R1.fastq.gz
mv BR1_Primary_S_S27_L002_R2_001.fastq.gz B73_Primary_S_br1_R2.fastq.gz
mv BR1_Seminal_G1_S28_L002_R1_001.fastq.gz B73_Seminal_G1_br1_R1.fastq.gz

mv BR1_Seminal_G1_S28_L002_R2_001.fastq.gz B73_Seminal_G1_br1_R2.fastq.gz
mv BR1_Seminal_S_S29_L002_R1_001.fastq.gz B73_Seminal_S_br1_R1.fastq.gz
mv BR1_Seminal_S_S29_L002_R2_001.fastq.gz B73_Seminal_S_br1_R2.fastq.gz
mv BR2_Primary_G1_S30_L002_R1_001.fastq.gz B73_Primary_G1_br2_R1.fastq.gz
mv BR2_Primary_G1_S30_L002_R2_001.fastq.gz B73_Primary_G1_br2_R2.fastq.gz

mv BR2_Primary_S_S31_L002_R1_001.fastq.gz B73_Primary_S_br2_R1.fastq.gz
mv BR2_Primary_S_S31_L002_R2_001.fastq.gz B73_Primary_S_br2_R2.fastq.gz
mv BR2_Seminal_G1_S32_L002_R1_001.fastq.gz B73_Seminal_G1_br2_R1.fastq.gz
mv BR2_Seminal_G1_S32_L002_R2_001.fastq.gz B73_Seminal_G1_br2_R2.fastq.gz
mv BR2_Seminal_S_S25_L002_R1_001.fastq.gz B73_Seminal_S_br2_R1.fastq.gz

mv BR2_Seminal_S_S25_L002_R2_001.fastq.gz B73_Seminal_S_br2_R2.fastq.gz
mv BR3_Primary_G1_S33_L002_R1_001.fastq.gz B73_Primary_G1_br3_R1.fastq.gz
mv BR3_Primary_G1_S33_L002_R2_001.fastq.gz B73_Primary_G1_br3_R2.fastq.gz
mv BR3_Primary_S_S34_L002_R1_001.fastq.gz B73_Primary_S_br3_R1.fastq.gz
mv BR3_Primary_S_S34_L002_R2_001.fastq.gz B73_Primary_S_br3_R2.fastq.gz

mv BR3_Seminal_G1_S35_L002_R1_001.fastq.gz B73_Seminal_G1_br3_R1.fastq.gz
mv BR3_Seminal_G1_S35_L002_R2_001.fastq.gz B73_Seminal_G1_br3_R2.fastq.gz
mv BR3_Seminal_S_S36_L002_R1_001.fastq.gz B73_Seminal_S_br3_R1.fastq.gz
mv BR3_Seminal_S_S36_L002_R2_001.fastq.gz B73_Seminal_S_br3_R2.fastq.gz
```

### Merge primary and seminal root samples
```bash

```

## EdU seqeuncing files
### Get fastq files from NCBI SRA for first round of sequencing
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_EdU_G1_br1 | xx | BR1_G1_S1_L001_R1_001.fastq.gz | BR1_G1_S1_L001_R2_001.fastq.gz
B73_EdU_G1_br2 | xx | BR2_G1_S3_L001_R1_001.fastq.gz | BR2_G1_S3_L001_R2_001.fastq.gz
B73_EdU_G1_br3 | xx | BR3_G1_S5_L001_R1_001.fastq.gz | BR3_G1_S5_L001_R2_001.fastq.gz
B73_EdU_S_br1 | xx | BR1_S_S2_L001_R1_001.fastq.gz | BR1_S_S2_L001_R2_001.fastq.gz
B73_EdU_S_br2 | xx | BR2_S_S4_L001_R1_001.fastq.gz | BR2_S_S4_L001_R2_001.fastq.gz
B73_EdU_S_br3 | xx | BR3_S_S6_L001_R1_001.fastq.gz | BR3_S_S6_L001_R2_001.fastq.gz

### Get fastq files from NCBI SRA for second round of sequencing
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_EdU_G1_br1 | xx | B73_EdU_SG1_BR1_G1_S1_L003_R1_001.fastq.gz | B73_EdU_SG1_BR1_G1_S1_L003_R2_001.fastq.gz
B73_EdU_G1_br2 | xx | B73_EdU_SG1_BR2_G1_S3_L003_R1_001.fastq.gz | B73_EdU_SG1_BR2_G1_S3_L003_R2_001.fastq.gz
B73_EdU_G1_br3 | xx | B73_EdU_SG1_BR3_G1_S5_L003_R1_001.fastq.gz | B73_EdU_SG1_BR3_G1_S5_L003_R2_001.fastq.gz
B73_EdU_S_br1 | xx | B73_EdU_SG1_BR1_S_S2_L003_R1_001.fastq.gz | B73_EdU_SG1_BR1_S_S2_L003_R2_001.fastq.gz
B73_EdU_S_br2 | xx | B73_EdU_SG1_BR2_S_S4_L003_R1_001.fastq.gz | B73_EdU_SG1_BR2_S_S4_L003_R2_001.fastq.gz
B73_EdU_S_br3 | xx | B73_EdU_SG1_BR3_S_S6_L003_R1_001.fastq.gz | B73_EdU_SG1_BR3_S_S6_L003_R2_001.fastq.gz

### Merge first and second round of sequencing:
```bash
cat B73_EdU_SG1_BR1_G1_S1_L003_R1_001.fastq.gz BR1_G1_S1_L001_R1_001.fastq.gz > B73_G1_br1_r1.fastq.gz
cat B73_EdU_SG1_BR1_G1_S1_L003_R2_001.fastq.gz BR1_G1_S1_L001_R2_001.fastq.gz > B73_G1_br1_r2.fastq.gz
cat B73_EdU_SG1_BR1_S_S2_L003_R1_001.fastq.gz BR1_S_S2_L001_R1_001.fastq.gz > B73_S_br1_r1.fastq.gz
cat B73_EdU_SG1_BR1_S_S2_L003_R2_001.fastq.gz BR1_S_S2_L001_R2_001.fastq.gz > B73_S_br1_r2.fastq.gz
cat B73_EdU_SG1_BR2_G1_S3_L003_R1_001.fastq.gz BR2_G1_S3_L001_R1_001.fastq.gz > B73_G1_br2_r1.fastq.gz
cat B73_EdU_SG1_BR2_G1_S3_L003_R2_001.fastq.gz BR2_G1_S3_L001_R2_001.fastq.gz > B73_G1_br2_r2.fastq.gz
cat B73_EdU_SG1_BR2_S_S4_L003_R1_001.fastq.gz BR2_S_S4_L001_R1_001.fastq.gz > B73_S_br2_r1.fastq.gz
cat B73_EdU_SG1_BR2_S_S4_L003_R2_001.fastq.gz BR2_S_S4_L001_R2_001.fastq.gz > B73_S_br2_r2.fastq.gz
cat B73_EdU_SG1_BR3_G1_S5_L003_R1_001.fastq.gz BR3_G1_S5_L001_R1_001.fastq.gz > B73_G1_br3_r1.fastq.gz
cat B73_EdU_SG1_BR3_G1_S5_L003_R2_001.fastq.gz BR3_G1_S5_L001_R2_001.fastq.gz > B73_G1_br3_r2.fastq.gz
cat B73_EdU_SG1_BR3_S_S6_L003_R1_001.fastq.gz BR3_S_S6_L001_R1_001.fastq.gz > B73_S_br3_r1.fastq.gz
cat B73_EdU_SG1_BR3_S_S6_L003_R2_001.fastq.gz BR3_S_S6_L001_R2_001.fastq.gz > B73_S_br3_r2.fastq.gz
```

## Repli-seq sequencing files
### Get fastq files from NCBI SRA for first round of sequencing
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
Early_br1 | --- | --- | ---
Early_br2 | --- | --- | ---
Early_br3 | --- | --- | ---
Mid_br1 | --- | --- | ---
Mid_br2 | --- | --- | ---
Mid_br3 | --- | --- | ---
Late_br1 | --- | --- | ---
Late_br2 | --- | --- | ---
Late_br3 | --- | --- | ---



## Trim fastq files
```bash
for x in *_r1.fastq.gz;
do
  trimmomatic PE -threads 48 -phred33 -validatePairs \
  -summary $SCRATCH/$(basename $x _r1.fastq.gz)_summary.txt \
  $SCRATCH/B73_SG1_nonEdU/02_merge_fastq/$(basename $x _r1.fastq.gz)_r1.fastq.gz \
  $SCRATCH/B73_SG1_nonEdU/02_merge_fastq/$(basename $x _r1.fastq.gz)_r2.fastq.gz \
  $SCRATCH/B73_SG1_nonEdU/03_trim_fastq/$(basename $x _r1.fastq.gz)_r1_trimmed.fastq \
  $SCRATCH/B73_SG1_nonEdU/03_trim_fastq/$(basename $x _r1.fastq.gz)_r1_trimmed_unpaired.fastq \
  $SCRATCH/B73_SG1_nonEdU/03_trim_fastq/$(basename $x _r1.fastq.gz)_r2_trimmed.fastq \
  $SCRATCH/B73_SG1_nonEdU/03_trim_fastq/$(basename $x _r1.fastq.gz)_r2_trimmed_unpaired.fastq \
  ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:5 TRAILING:5 MINLEN:20 ;
done
```

## Run Fastqc on trimmed fastq files
```bash
for x in *.fastq
do
  fastqc_latest.sif fastqc -t 48 $x -o . ;
done 
```



# Map and filter trimmed reads 
## Map
Make bowtie 2 index:
```bash
bowtie2-build --threads 48 Zm-B73-REFERENCE-NAM-5.0.fa Zm-B73-REFERENCE-NAM-5.0_bowtie_index
```
Map trimed reads to B73 v5 NAM genome:
```bash
for x in *_r1_trimmed.fastq;
do 
  DO_NOT_USE.bowtie2_v2.5.1.sif bowtie2 -p 48 --very-sensitive \
  -x Zm-B73-REFERENCE-NAM-5.0_bowtie_index \
  -1 $(basename $x _r1_trimmed.fastq)_r1_trimmed.fastq  \
  -2 $(basename $x _r1_trimmed.fastq)_r2_trimmed.fastq | \
  samtools view -@ 48 -bS - > $(basename $x _r1_trimmed.fastq)_withscaffold_prefiltered.bam " ;
done 
```


## Sort bam file
```bash
for x in *_withscaffold_prefiltered.bam;
do
  samtools sort -@ 48 $x > $(basename $x .bam)_sorted.bam ;
done
```

## Index bam files:
Get only reads in chromosome, write commands:
```bash
for x in *_withscaffold_prefiltered_sorted.bam ;
do
  samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done
```

## Remove duplicated reads
```bash
for x in *_sorted.bam
do
  sambamba markdup --io-buffer-size=1920 --overflow-list-size=1000000 -r --nthreads=48 $x $(basename $x _withscaffold_prefiltered_sorted.bam)_rmDup.bam ;
done 
```


## Get only reads that are properly paired and filter out reads that are MAPQ <6 (ie MAPQ<=5), and sort bam file
```bash
for x in *_rmDup.bam ;
do
  time apptainer exec samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -bf 0x2 $x | apptainer exec samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -q 6 -b - | apptainer exec samtools_v1.9-4-deb_cv1.sif samtools sort -@ 48 - > $(basename $x _rmDup.bam)_filtered.bam ;
done
```

## Remove reads from scaffold and low coverage droplist

Create bedfile of KEEP regions, ie the chromosomes minus droplist
```bash
bedtools complement -i ../B73_noCovPercent_gt60_10kb.bed -g ../B73_chrOnly_size.txt > B73_keep_regions.bed
```

Remove scaffolds and mappability droplist
```bash
for x in *_filtered.bam ;
do
  samtools view -bh -@ 48 -L B73_keep_regions.bed $x > $(basename $x _filtered.bam)_droplistRm.bam ;
done
```

## Sort final bam files
```bash
for x in *_droplistRm.bam
do
  samtools sort -@ 48 $x > $(basename $x .bam)_final.bam ;
done
```

## Index final bam files
```bash
for x in *_final.bam
do
  samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done
```

## Calculate number of reads

Number of reads in raw merged fastq files:
```bash
for x in *_r1.fastq.gz ;
do
  zgrep "@A00600" $x | wc -l > reads_merged_$(basename $x .fastq.gz)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' reads_merged_$(basename $x .fastq.gz)_tmp.txt > reads_merged_$(basename $x .fastq.gz)_names.txt ;
  rm *_tmp.txt ;
done

cat reads_merged_*_names.txt > reads_merged_all.txt 
rm reads_merged_*_names.txt
```

Number of reads in trimmed fastq files:
```bash
for x in *_r1_trimmed.fastq ;
do
  zgrep "@A00600" $x | wc -l > reads_trimmed_$(basename $x .fastq.gz)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' reads_trimmed_$(basename $x .fastq.gz)_tmp.txt > reads_trimmed_$(basename $x .fastq.gz)_names.txt ;
done

cat reads_trimmed_*_names.txt > reads_trimmed_all.txt 
rm reads_trimmed_*_names.txt *temp*
```


Number of reads post trimming:
```bash
for x in *_r1_trimmed.fastq ;
do
  zgrep "@A00600" $x | wc -l > reads_trimmed_$(basename $x .fastq.gz)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' reads_trimmed_$(basename $x .fastq.gz)_tmp.txt > reads_trimmed_$(basename $x .fastq.gz)_names.txt ;
done

cat reads_trimmed_*_names.txt > reads_trimmed_all.txt 
rm reads_trimmed_*_names.txt *temp*
```


Number of reads mapped:
```bash
for x in *_withscaffold_prefiltered_sorted.bam
do
  samtools view -@ 48 -c $x > $(basename $x .bam)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x .bam)_tmp.txt > mapped_reads_inital_$(basename $x .bam).txt ;	
done

cat mapped_reads_inital_*txt > mapped_reads_inital_all.txt
```

Number of reads after duplicates removed:
```bash
for x in *_rmDup.bam
do
  samtools view -@ 48 -c $x > $(basename $x _rmDup.bam)_rmDup_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x _rmDup.bam)_rmDup_tmp.txt > mapped_reads_rmDup_$(basename $x .bam).txt ;
  rm *_rmDup_tmp.txt ;
done

cat mapped_reads_rmDup_*txt > mapped_reads_rmDup_all.txt
```

Number of reads after filtering:
```bash
for x in *_filtered.bam
do
  samtools view -@ 48 -c $x > $(basename $x .bam)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x .bam)_tmp.txt > mapped_reads_filtered_$(basename $x .bam).txt ;
  rm *_tmp.txt ;
done

cat mapped_reads_filtered_*txt > mapped_reads_filtered_all.txt
```

Number of reads after removing scaffolds and low mappability regions:
```bash
for x in *_droplistRm_final.bam
do
  samtools view -@ 48 -c $x > $(basename $x .bam)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x .bam)_tmp.txt > mapped_reads_final_$(basename $x .bam).txt ;
  rm *_tmp.txt ;
done

cat mapped_reads_final_*txt > mapped_reads_final_all.txt
```



# High and low coverage droplists
## High coverage droplist 
Various cutoffs were compared, but the highest coverage 0.25% bins from each biorep will be added to droplist. 
```bash
# Count number of 10kb bins and bp in B73 genome:
wc -l ./B73_G1_br1_coverage_10kb.bedgraph
213190  # 10 kb bins in genome

# Calculate the number of bins in 0.25%:
213190 * 0.0025 = 532  # 532 bins are 0.25% of the bins in the genome

# Collect the coordinate information for the highest 0.25% bins according to coverage from each G1 biorep:
for x in *G1*_coverage_10kb.bedgraph
do
  sort -n -k 4 -r $x | head -n 532 | cut -f 1-3 | sort -k 1,1 -k2,2n > $(basename $x _coverage_10kb.bedgraph)_individual_high_droplist.bed ;
done

# Combine all G1 biorep bins to create pooled droplist:
cat *_individual_high_droplist.bed | sort -k 1,1 -k2,2n | bedtools_2.31.0.sif bedtools merge -i - > B73_pooled_high_droplist.bed
```

## Low coverage droplist
See genome mappability repository ([link](https://github.com/ewheeler7/genome_mappability)).

The B73 NAM v5 genome was used to create artifical reads, which were then mapped back to the refrence genome and filtered. Any 10 kb bin with ≥ 60% of base pairs not covered by any artifical reads was added to the low coverage droplist.

This step is important because in the S/G1 approach low read coverage indicates later replication timing. And we want to call a region late because it has low reads, not because it has low inherent mappability. 

## Merge high and low droplists
```bash
cat B73_pooled_high_droplist.bed B73_noCovPercent_gt60_10kb.bed | sort -k 1V,1 -k 2n,2 | bedtools_2.31.0.sif bedtools merge -i - > B73_final_droplist.bed
```


--> This is the point where repliseq data was run through repliscan



# Normalize reads
Using bam files with high and low coverage droplists removed, 1X normalize reads using deeptools RPCG mode.

Effective genome size is calculated by the number of bp in the chromosomal genome - number of bp in final droplist.
```bash
for x in *final.bam;
do
  bamCoverage -p 48 -b $x -o $(basename $x _droplistRm_final.bam)_10kb_1X_norm.bedgraph -of bedgraph -bs 10000 --effectiveGenomeSize 2071945143 --normalizeUsing RPGC --blackListFileName B73_final_droplist.bed --ignoreDuplicates --extendReads ;
done
```
****should I have removed schaffolds before normalizing?????? *****

Need to re-map normalised data, because deeptools merges nextdoor bins with identical values, but because next step is to make the S/G1 ratio, I need all data in 10kb bins. 

Remove placeholder scaffold bins first, then map, re-sort, create bw, then remove intermediate files
(The scaffold bins and plastid coordinates are in this file, but there are no reads in them because they were removed from bam files. It is just for convience that we remove them at this point)
```bash
for x in *_10kb_1X_norm.bedgraph ;
do
	grep -v "scaf" $x | grep -v "chrM" - | grep -v "chrC" - > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_rmScaf.bedgraph;
	bedtools map -c 4 -o sum -null 0 -a B73_10kb.bed -b $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_rmScaf.bedgraph  > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bedgraph ;
	sort -k1,1 -k2,2n $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bedgraph > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_sorted.bedgraph ;
	bedGraphToBigWig $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_sorted.bedgraph B73_chrOnly_size.txt $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bw ;
	rm *_rmScaf.bedgraph *_sorted.bedgraph ;
done
```

# S/G1 ratio
## Pipeline:
- Make S/G1 ratio for each 10kb bin, using individual biorep G1s
- Make S/G1 ratio for each 10kb bin, using an average G1
- Find all non-droplist zeros in all bioreps, create common "zero coverage droplist", remove this droplist from all bioreps
- Remove "standalone" bins of data, ie if one bin of data is flanked by droplists, drop that bin as well (make sure identical between bioreps)


## Make S/G1 ratio for each 10kb bin, using individual biorep G1s:
```bash
paste B73_S_br1_1X_10kb.bedgraph ../11_1x_normalize/B73_G1_br1_1X_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t""0" ; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_br1_10kb.bedgraph
paste ../11_1x_normalize/B73_S_br2_1X_10kb.bedgraph ../11_1x_normalize/B73_G1_br2_1X_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t""0" ; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_br2_10kb.bedgraph
paste ../11_1x_normalize/B73_S_br3_1X_10kb.bedgraph ../11_1x_normalize/B73_G1_br3_1X_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t""0" ; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_br3_10kb.bedgraph





add Non-Edu





```


## Make S/G1 ratio for each 10kb bin

Using all biorep G1 samples, make an average G1
```bash
paste ../11_1x_normalize/B73_G1_br1_1X_10kb.bedgraph ../11_1x_normalize/B73_G1_br2_1X_10kb.bedgraph ../11_1x_normalize/B73_G1_br3_1X_10kb.bedgraph | \
awk '{ print $1"\t"$2"\t"$3"\t"($4+$8+$12)/3 }' - > B73_avgG1_10kb.bedgraph
```

Calculate the S/G1 ratio, using the average G1 made above
```bash
paste ../11_1x_normalize/B73_S_br1_1X_10kb.bedgraph B73_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_avgG1_br1_10kb.bedgraph
paste ../11_1x_normalize/B73_S_br2_1X_10kb.bedgraph B73_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_avgG1_br2_10kb.bedgraph
paste ../11_1x_normalize/B73_S_br3_1X_10kb.bedgraph B73_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_avgG1_br3_10kb.bedgraph



Add non-EdU


```

Check: How many 10kb bins have S/G1=0, that are NOT from the droplist?
```bash
for x in ./*10kb.bedgraph ;
do
	ls $x ;
	bedtools subtract -a $x -b B73_final_droplist.bed | awk '{ if ($4==0) print $0 }' | wc -l ;
done

./B73_avgG1_10kb.bedgraph			0
./B73_SG1_ratio_avgG1_br1_10kb.bedgraph 	6
./B73_SG1_ratio_avgG1_br2_10kb.bedgraph 	4
./B73_SG1_ratio_avgG1_br3_10kb.bedgraph 	4
./B73_SG1_ratio_br1_10kb.bedgraph 		7
./B73_SG1_ratio_br2_10kb.bedgraph 		9
./B73_SG1_ratio_br3_10kb.bedgraph 		6
```
This means there are "true" 0 places?




### Find all zeros in all bioreps, this includes the droplist and NON-DROPLIST ZEROS (see above) and create common "zero coverage droplist", remove this droplist from all bioreps

- Use all bioreps to add to droplist = ../10_droplist/B73_final_droplist.bed + non-droplist zeros = B73_droplist_and_zeroCov.bed
- The ../10_droplist/B73_final_droplist.bed regions are already 0's in the following bedgraphs
```bash
cat B73_SG1_ratio_avgG1_br*_10kb.bedgraph | \
awk '{ if ($4==0) print $0}' - | sort -k 1V,1 -k 2n,2 | \
../bedtools_2.31.0.sif bedtools merge -i - > B73_droplist_and_zeroCov.bed

# Remove the non-droplist-zeros from the bioreps:
for x in ./B73_SG1_ratio_avgG1_br*_10kb.bedgraph ; do
	../bedtools_2.31.0.sif bedtools subtract -a $x -b B73_droplist_and_zeroCov.bed > $(basename $x .bedgraph)_dropNonDropZero.bedgraph ;
done
```

### Remove "standalone" bins of data 
If a single 10kb of data is flanked on both sides by a droplist, drop that bin in that sample and in all bioreps:
```bash
# Create total droplist regions bedfile, includes the "standalone" regions:
bedtools complement -i B73_SG1_ratio_avgG1_br1_10kb_dropNonDropZero.bedgraph -g ../B73_chrOnly_size.txt | \
bedtools merge -d 15000 -i - > B73_dropStandalone.bed

# Remove the standalone regions, use a file that contains ../10_droplist/B73_final_droplist.bed + non-droplist zeros + standalone bins
for x in ./*_dropNonDropZero.bedgraph ;
do
	bedtools subtract -a $x -b B73_dropStandalone.bed > $(basename $x _dropNonDropZero.bedgraph)_final.bedgraph ;
done
```


# Haar Wavelet smoothing

Haar wavelet level 2 smooth each biorep.
The input data for wavelets has droplist removed, non-droplist zeros removed, and "standalone" bins removed, because if not those regions will be seen as a valid 0 and bring down the real values around it. 
```bash
for x in ../12_SG1_ratio/*br*_final.bedgraph ;
do
	# Get only bed location info from final SG1 ratio bedgraph:
	cut -f 1-3 $x > $(basename $x .bedgraph).bedgraph.loc ;
	# Get only value info from final SG1 ratio bedgraph:
	cut -f 4 $x > $(basename $x .bedgraph).bedgraph.val ;
	# Run wavelet smoothing at level two on droplist removed values:
	$WORK/wavelets/bin/wavelets --level 2 --to-stdout --boundary reflected --filter Haar $(basename $x .bedgraph).bedgraph.val > $(basename $x .bedgraph).bedgraph.smooth;
	# Combine the location info from the final SG1 ratio data with the wavelet smoothed output:
	paste $(basename $x .bedgraph).bedgraph.loc $(basename $x .bedgraph).bedgraph.smooth > $(basename $x .bedgraph)_H2.bedgraph ;
	# Sort for kentUtils:
	sort -k1,1 -k2,2n $(basename $x .bedgraph)_H2.bedgraph > $(basename $x .bedgraph)_H2_sorted.bedgraph ;
	# Create bigwig file:
	bedGraphToBigWig $(basename $x .bedgraph)_H2_sorted.bedgraph ../B73_chrOnly_size.txt $(basename $x .bedgraph)_H2.bw ;
	# Remove temporary files:
	rm *val *loc *smooth *_sorted.bedgraph ;
done
```

Make an average of the bioreps profile:
```bash
bedtools map -c 4 -o mean -null 0 -a ../B73_10kb.bed -b B73_SG1_ratio_avgG1_br1_10kb_final_H2.bedgraph | \
bedtools map -c 4 -o mean -null 0 -a - -b B73_SG1_ratio_avgG1_br2_10kb_final_H2.bedgraph | \
bedtools map -c 4 -o mean -null 0 -a - -b B73_SG1_ratio_avgG1_br3_10kb_final_H2.bedgraph | \
awk '{ if ($4>0) print $1"\t"$2"\t"$3"\t"($4+$5+$6)/3}' - > B73_SG1_ratio_avgG1_avg_10kb_final_H2.bedgraph


Add non-edu



```








