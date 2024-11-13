# Pipeline
- Read pre-processing
- Map and filter trimmed reads
- High and low coverage droplists
- Normalize reads
- S/G1 ratio
- Haar Wavelet smoothing

# Software required
Software | File name | Link 
--- | --- | --- 
trimomatic | trimmomatic_latest.sif | link 
fastqc | fastqc_latest.sif | 
samtools | samtools_v1.9-4-deb_cv1.sif | 
bowtie2 | DO_NOT_USE.bowtie2_v2.5.1.sif |
bedtools | bedtools_2.31.0.sif |
kentUtils | kentutils_1.04.00.sif |
deeptools | deeptools_3.5.4--pyhdfd78af_1.sif |




# Read pre-processing

## Non-EdU
### Get fastq files from primary roots:
Sample  | Accesion | Fw | Rv
--- | --- | --- | --- 
B73_EdU_G1_br1 | xx | BR1_G1_S1_L001_R1_001.fastq.gz | BR1_G1_S1_L001_R2_001.fastq.gz

### Get fastq files from seminal roots:
Sample  | Accesion | Fw | Rv
--- | --- | --- | --- 
B73_EdU_G1_br1 | xx | BR1_G1_S1_L001_R1_001.fastq.gz | BR1_G1_S1_L001_R2_001.fastq.gz


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

### Merge primary and seminal root samples:
```bash

```

## EdU
### Get fastq files from first round of sequencing
Sample  | Accesion | Fw | Rv
--- | --- | --- | --- 
B73_EdU_G1_br1 | xx | BR1_G1_S1_L001_R1_001.fastq.gz | BR1_G1_S1_L001_R2_001.fastq.gz
B73_EdU_G1_br2 | xx | BR2_G1_S3_L001_R1_001.fastq.gz | BR2_G1_S3_L001_R2_001.fastq.gz
B73_EdU_G1_br3 | xx | BR3_G1_S5_L001_R1_001.fastq.gz | BR3_G1_S5_L001_R2_001.fastq.gz
B73_EdU_S_br1 | xx | BR1_S_S2_L001_R1_001.fastq.gz | BR1_S_S2_L001_R2_001.fastq.gz
B73_EdU_S_br2 | xx | BR2_S_S4_L001_R1_001.fastq.gz | BR2_S_S4_L001_R2_001.fastq.gz
B73_EdU_S_br3 | xx | BR3_S_S6_L001_R1_001.fastq.gz | BR3_S_S6_L001_R2_001.fastq.gz

### Get fastq files from second round of sequencing
Sample  | Accesion | Fw | Rv
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

## Trim fastq files
```bash
for x in *_r1.fastq.gz;
do
  ../trimmomatic_latest.sif trimmomatic PE -threads 48 -phred33 -validatePairs \
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
DO_NOT_USE.bowtie2_v2.5.1.sif bowtie2-build --threads 48 Zm-B73-REFERENCE-NAM-5.0.fa Zm-B73-REFERENCE-NAM-5.0_bowtie_index
```
Map trimed reads to B73 v5 NAM genome:
```bash
for x in $SCRATCH/B73_SG1_nonEdU/03_trim_fastq/*_r1_trimmed.fastq;
do 
  DO_NOT_USE.bowtie2_v2.5.1.sif bowtie2 -p 48 --very-sensitive \
  -x Zm-B73-REFERENCE-NAM-5.0_bowtie_index \
  -1 $(basename $x _r1_trimmed.fastq)_r1_trimmed.fastq  \
  -2 $(basename $x _r1_trimmed.fastq)_r2_trimmed.fastq | \
  time apptainer exec samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -bS - > $(basename $x _r1_trimmed.fastq)_withscaffold_prefiltered.bam " ;
done 
```


## Sort bam file
```bash
for x in *_withscaffold_prefiltered.bam;
do
  time apptainer exec samtools_v1.9-4-deb_cv1.sif samtools sort -@ 48 $x > $(basename $x .bam)_sorted.bam ;
done
```

## Index bam files:
Get only reads in chromosome, write commands:
```bash
for x in *_withscaffold_prefiltered_sorted.bam ;
do
  time apptainer exec samtools_v1.9-4-deb_cv1.sif samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done
```

## Remove duplicated reads:
```bash
for x in *_sorted.bam
do
  time apptainer exec sambamba_1.0.0.sif sambamba markdup --io-buffer-size=1920 --overflow-list-size=1000000 -r --nthreads=48 $x $(basename $x _withscaffold_prefiltered_sorted.bam)_rmDup.bam ;
done 
```


## Get only reads that are properly paired and filter out reads that are MAPQ <6 (ie MAPQ<=5), and sort bam file
```bash
for x in $SCRATCH/B73_SG1_nonEdU/05_map/*_rmDup.bam ;
do
  time apptainer exec samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -bf 0x2 $x | apptainer exec samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -q 6 -b - | apptainer exec samtools_v1.9-4-deb_cv1.sif samtools sort -@ 48 - > $(basename $x _rmDup.bam)_filtered.bam ;
done
```

## Remove reads from scaffold and low coverage droplist

### Create bedfile of KEEP regions, ie the chromosomes minus droplist
```bash
bedtools_2.31.0.sif bedtools complement -i ../B73_noCovPercent_gt60_10kb.bed -g ../B73_chrOnly_size.txt > B73_keep_regions.bed
```

### Remove scaffolds and mappability droplist
for x in *_filtered.bam ;
do
  samtools_v1.9-4-deb_cv1.sif samtools view -bh -@ 48 -L B73_keep_regions.bed $x > $(basename $x _filtered.bam)_droplistRm.bam ;
done

## Sort final bam files
```bash
for x in *_droplistRm.bam
do
  samtools_v1.9-4-deb_cv1.sif samtools sort -@ 48 $x > $(basename $x .bam)_final.bam ;
done
```

## Index final bam files
```bash
for x in *_final.bam
do
  samtools_v1.9-4-deb_cv1.sif samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done
```

## Calculate number of reads

Number of reads in raw merged fastq files:
```bash
for x in ../02_merge_fastq/*_r1.fastq.gz ;
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
for x in ../03_trim_fastq/*_r1_trimmed.fastq ;
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
  samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -c $x > $(basename $x .bam)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x .bam)_tmp.txt > mapped_reads_inital_$(basename $x .bam).txt ;	
done

cat mapped_reads_inital_*txt > mapped_reads_inital_all.txt
```

Number of reads after duplicates removed:
```bash
for x in *_rmDup.bam
do
  samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -c $x > $(basename $x _rmDup.bam)_rmDup_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x _rmDup.bam)_rmDup_tmp.txt > mapped_reads_rmDup_$(basename $x .bam).txt ;
  rm *_rmDup_tmp.txt ;
done

cat mapped_reads_rmDup_*txt > mapped_reads_rmDup_all.txt
```

Number of reads after filtering:
```bash
for x in *_filtered.bam
do
  samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -c $x > $(basename $x .bam)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x .bam)_tmp.txt > mapped_reads_filtered_$(basename $x .bam).txt ;
  rm *_tmp.txt ;
done

cat mapped_reads_filtered_*txt > mapped_reads_filtered_all.txt
```

Number of reads after removing scaffolds and low mappability regions:
```bash
for x in *_droplistRm_final.bam
do
  samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -c $x > $(basename $x .bam)_tmp.txt ;
  awk '{print FILENAME"\t"$1}' $(basename $x .bam)_tmp.txt > mapped_reads_final_$(basename $x .bam).txt ;
  rm *_tmp.txt ;
done

cat mapped_reads_final_*txt > mapped_reads_final_all.txt
```

Create table of mapped reads:
*this is only for non_EdU!!
```bash
paste reads_primary_all.txt reads_seminal_all.txt reads_merged_all.txt mapped_reads_inital_all.txt mapped_reads_rmDup_all.txt mapped_reads_filtered_all.txt mapped_reads_final_all.txt | cut -f 2,4,6,8,10,12,14 > calc_reads_table_tmp.txt
echo -e "primary\tseminal\tmerged\tinital\trmDup\tfiltered\tfinal" | cat - calc_reads_table_tmp.txt > B73_calc_reads_table_tmp.txt
echo -e "Sample\nB73_G1_br1\nB73_G1_br2\nB73_G1_br3\nB73_S_br1\nB73_S_br2\nB73_S_br3" | paste - B73_calc_reads_table_tmp.txt > B73_calc_reads.txt
rm *_tmp.txt
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
See genome mappability repository (link)
The filtered read coverage was calculated in 10 kb bins and any bin with ≥ 60% ZERO coverage was added to the low coverage droplist
This step is important because in the S/G1 approach low read coverage indicates later replication timing. And we want to call a region late because it has low reads, not because it has low inherent mappability. 

## Merge high and low droplists
```bash
cat B73_pooled_high_droplist.bed B73_noCovPercent_gt60_10kb.bed | sort -k 1V,1 -k 2n,2 | bedtools_2.31.0.sif bedtools merge -i - > B73_final_droplist.bed
```

## Remove standalone bins
*** didn't do until after making ratio, but probably should have done this step here??****


# Normalize reads
Using bam files with high and low coverage droplists removed, 1X normalize reads:
effective genome size is calculated by the number of bp in the chromosomal genome - number of bp in final droplist:
```bash
for x in *final.bam;
do
  time apptainer exec deeptools_3.5.4--pyhdfd78af_1.sif bamCoverage -p 48 -b $x -o $(basename $x _droplistRm_final.bam)_10kb_1X_norm.bedgraph -of bedgraph -bs 10000 --effectiveGenomeSize 2071945143 --normalizeUsing RPGC --blackListFileName B73_final_droplist.bed --ignoreDuplicates --extendReads ;
done
```
****should I have removed schaffolds before normalizing?????? *****

Need to re-map normalised data 
Remove placeholder scaffold bins first, then map, re-sort, create bw, then remove intermediate files
```bash
for x in *_10kb_1X_norm.bedgraph ;
do
	grep -v "scaf" $x | grep -v "chrM" - | grep -v "chrC" - > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_rmScaf.bedgraph;
	time apptainer exec bedtools_2.31.0.sif bedtools map -c 4 -o sum -null 0 -a B73_10kb.bed -b $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_rmScaf.bedgraph  > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bedgraph ;
	sort -k1,1 -k2,2n $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bedgraph > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_sorted.bedgraph ;
	time apptainer exec kentutils_1.04.00.sif bedGraphToBigWig $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_sorted.bedgraph B73_chrOnly_size.txt $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bw ;
	rm *_rmScaf.bedgraph *_sorted.bedgraph ;
done
```


# S/G1 ratio


# Haar Wavelet smoothing
## create an average biorep file
