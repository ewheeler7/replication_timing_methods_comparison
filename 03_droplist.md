# Define high coverage and low coverage droplists
Defining regions of high coverage is important because localized spikes in read coverage are often due to 
collapsed terminal repeats and therefore not biological relevant but will skew the read normalization steps. 
So removing these spikes is necessary. 

Defining regions of low coverage is important becase in the S/G1 ratio method for measuring DNA replication timing a depletion of reads indicates later 
replication time. However, we do not want to call a region late replicating if the real reason for low reads is due to poor mappability of the area.
The pipeline for our approach to this mabbability issue is under the repository "genome_mappability" ([link](https://github.com/ewheeler7/genome_mappability))



## High coverage droplist 
The highest coverage 0.25% bins from each biorep will be added to droplist. Various percent cutoffs were analyized, and 0.25% was the best option.
This step is done for each experiment separately, ie the top 0.25% from Repli-seq G1s are pooled togehter, while the EdU-G1s are pooled separately, and the non-EdU-G1s are also pooled separately.


**work on this want to make experiments easier to read**
```bash
# Count number of 10kb bins in B73 genome
wc -l B73_10kb.bed
213190  # 10 kb bins in genome

# Calculate the number of bins in 0.25%:
213190 * 0.0025 = 532  # 532 bins are 0.25% of the bins in the genome

-> Top 0.25% bins = 213190 * 0.0025 = 532 

# For Replieq experiment: Collect the coordinate information for the highest 0.25% bins according to coverage from each G1 biorep:
for x in *G1*_coverage_10kb.bedgraph
do
  sort -n -k 4 -r $x | head -n 532 | cut -f 1-3 | sort -k 1,1 -k2,2n > $(basename $x _coverage_10kb.bedgraph)_individual_high_droplist.bed ;
done

# For EdU-S/G1 experiment: Collect the coordinate information for the highest 0.25% bins according to coverage from each G1 biorep:
for x in *G1*_coverage_10kb.bedgraph
do
  sort -n -k 4 -r $x | head -n 532 | cut -f 1-3 | sort -k 1,1 -k2,2n > $(basename $x _coverage_10kb.bedgraph)_individual_high_droplist.bed ;
done

# For S/G1 experiment: Collect the coordinate information for the highest 0.25% bins according to coverage from each G1 biorep:
for x in *G1*_coverage_10kb.bedgraph
do
  sort -n -k 4 -r $x | head -n 532 | cut -f 1-3 | sort -k 1,1 -k2,2n > $(basename $x _coverage_10kb.bedgraph)_individual_high_droplist.bed ;
done


# For Repli-seq experiment: Combine all G1 biorep bins to create pooled droplist:
cat *_individual_high_droplist.bed | sort -k 1,1 -k2,2n | bedtools_2.31.0.sif bedtools merge -i - > B73_pooled_high_droplist.bed

# For EdU-S/G1 experiment: Combine all G1 biorep bins to create pooled droplist:
cat *_individual_high_droplist.bed | sort -k 1,1 -k2,2n | bedtools_2.31.0.sif bedtools merge -i - > B73_pooled_high_droplist.bed

# For S/G1 experiment: Combine all G1 biorep bins to create pooled droplist:
cat *_individual_high_droplist.bed | sort -k 1,1 -k2,2n | bedtools_2.31.0.sif bedtools merge -i - > B73_pooled_high_droplist.bed


```


## Low coverage droplist
See genome mappability repository ([link](https://github.com/ewheeler7/genome_mappability)).

Breifly, the B73 NAM v5 genome was used to create artifical reads, which were then mapped back to the refrence genome and filtered. 
Any 10 kb bin with ≥ 60% of base pairs not covered by any artifical reads was added to the low coverage droplist.

The file name with ≥ 60% of base pairs not covered by any artifical reads is called  "B73_noCovPercent_gt60_10kb.bed"

## Merge high and low droplists

This was done separatly for each experiment:

```bash
# For Repli-seq experiment:
cat B73_pooled_high_droplist.bed B73_noCovPercent_gt60_10kb.bed | sort -k 1V,1 -k 2n,2 | bedtools_2.31.0.sif bedtools merge -i - > B73_final_droplist.bed

# For EdU-S/G1 experiment:
cat B73_pooled_high_droplist.bed B73_noCovPercent_gt60_10kb.bed | sort -k 1V,1 -k 2n,2 | bedtools_2.31.0.sif bedtools merge -i - > B73_final_droplist.bed

# For S/G1 experiment:
cat B73_pooled_high_droplist.bed B73_noCovPercent_gt60_10kb.bed | sort -k 1V,1 -k 2n,2 | bedtools_2.31.0.sif bedtools merge -i - > B73_final_droplist.bed



```

## Remove droplist regions from mapped bam files:


```bash





```

-> Late has many more reads removed - is this a problem? or does late biologically have more reads in G1 spike regions? 20-30% less????



## Index bam files:
```bash
for x in ./*_E.bam ;
do
	samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done

for x in ./*_M.bam ;
do
	samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done

for x in ./*_L.bam ;
do
	samtools index -@ 48 $x > $(basename $x .bam).bam.bai ;
done














```








