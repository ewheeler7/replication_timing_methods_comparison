# Calculate S/G1 ratio
Only the EdU-S/G1 and S/G1 experiment data is processed from this point. The Repli-seq data is completely processed after runing Repliscan.

Pipeline:
- normalize mapped reads
- make average of G1 biorep data 
- make S/G1 ratio between each S biorep and the average G1 data


## Normalize reads
Using bam files with the high and low coverage droplists removed, 1X normalize reads using deeptools RPCG mode.
Note: 1X normalization requires an effective genome size. 

We use this calculation: Effective genome size = (# of bp in the chromosomal genome assembly - # of bp in final droplist)
```bash
for x in *final.bam;
do
  bamCoverage -p 48 -b $x -o $(basename $x _droplistRm_final.bam)_10kb_1X_norm.bedgraph -of bedgraph -bs 10000 --effectiveGenomeSize 2071945143 --normalizeUsing RPGC --blackListFileName B73_final_droplist.bed --ignoreDuplicates --extendReads ;
done
```

Remove placeholder scaffold bins, re-map data to 10kb bins, sort the bedgraph file, and remove intermediate files. 
Deeptools merges nextdoor bins with identical values, so we need to re-map normalised data into 10 kb bins so the G1 and S bins align.
(The scaffold bins and plastid coordinates are in this file, but there are no reads in them because they were removed from bam files. It is just for convience that we remove them at this point)
```bash
for x in *_10kb_1X_norm.bedgraph ;
do
	grep -v "scaf" $x | grep -v "chrM" - | grep -v "chrC" - > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_rmScaf.bedgraph;
	bedtools map -c 4 -o sum -null 0 -a B73_10kb.bed -b $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb_rmScaf.bedgraph | sort -k1,1 -k2,2n > $(basename $x _10kb_1X_norm.bedgraph)_1X_10kb.bedgraph ;
	rm *_rmScaf.bedgraph;
done
```

## S/G1 ratio
Pipeline:
- Create an average G1 file using all G1 bioreps
- Make S/G1 ratio for each 10kb bin
- Find all non-droplist zeros in all bioreps, create common "zero coverage droplist", remove this droplist from all bioreps
- Remove "standalone" bins of data, ie if one bin of data is flanked by droplists, drop that bin as well (make sure identical between bioreps)


Using all biorep G1 samples, make an average G1:
```bash
# For EdU-S/G1
paste B73_EdU_G1_br1_1X_10kb.bedgraph B73_EdU_G1_br2_1X_10kb.bedgraph B73_EdU_G1_br3_1X_10kb.bedgraph | \
awk '{ print $1"\t"$2"\t"$3"\t"($4+$8+$12)/3 }' - > B73_EdU_avgG1_10kb.bedgraph

# For S/G1
paste B73_G1_br1_1X_10kb.bedgraph B73_G1_br2_1X_10kb.bedgraph B73_G1_br3_1X_10kb.bedgraph | \
awk '{ print $1"\t"$2"\t"$3"\t"($4+$8+$12)/3 }' - > B73_avgG1_10kb.bedgraph
```

Calculate the S/G1 ratio, using the average G1 sample made above. If either S or G1 has a value of zero, return a 0 as the ratio value:
```bash
# For EdU-S/G1
paste B73_EdU_S_br1_1X_10kb.bedgraph B73_EdU_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_EdU_SG1_ratio_avgG1_br1_10kb.bedgraph
paste B73_EdU_S_br2_1X_10kb.bedgraph B73_EdU_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_EdU_SG1_ratio_avgG1_br2_10kb.bedgraph
paste B73_EdU_S_br3_1X_10kb.bedgraph B73_EdU_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_EdU_SG1_ratio_avgG1_br3_10kb.bedgraph

# For S/G1
paste B73_S_br1_1X_10kb.bedgraph B73_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_avgG1_br1_10kb.bedgraph
paste B73_S_br2_1X_10kb.bedgraph B73_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_avgG1_br2_10kb.bedgraph
paste B73_S_br3_1X_10kb.bedgraph B73_avgG1_10kb.bedgraph | awk '{ if ($4==0 || $8==0) print $1"\t"$2"\t"$3"\t"0; else print $1"\t"$2"\t"$3"\t"$4/$8 }' > B73_SG1_ratio_avgG1_br3_10kb.bedgraph
```



### Find all zeros in all bioreps, this includes the droplist and NON-DROPLIST ZEROS (see above) and create common "zero coverage droplist", remove this droplist from all bioreps

- Use all bioreps to add to droplist = B73_final_droplist.bed + non-droplist zeros = B73_droplist_and_zeroCov.bed
- The B73_final_droplist.bed regions are already 0's in the following bedgraphs
```bash
# For EdU-S/G1
cat B73_EdU_SG1_ratio_avgG1_br*_10kb.bedgraph | \
awk '{ if ($4==0) print $0}' - | sort -k 1V,1 -k 2n,2 | \
bedtools merge -i - > B73_EdU_droplist_and_zeroCov.bed

## Remove the non-droplist-zeros from the bioreps:
for x in ./B73_EdU_SG1_ratio_avgG1_br*_10kb.bedgraph ; do
	bedtools subtract -a $x -b B73_EdU_droplist_and_zeroCov.bed > $(basename $x .bedgraph)_dropNonDropZero.bedgraph ;
done



# For S/G1
cat B73_SG1_ratio_avgG1_br*_10kb.bedgraph | \
awk '{ if ($4==0) print $0}' - | sort -k 1V,1 -k 2n,2 | \
bedtools merge -i - > B73_droplist_and_zeroCov.bed

## Remove the non-droplist-zeros from the bioreps:
for x in ./B73_SG1_ratio_avgG1_br*_10kb.bedgraph ; do
	bedtools subtract -a $x -b B73_droplist_and_zeroCov.bed > $(basename $x .bedgraph)_dropNonDropZero.bedgraph ;
done
```

### Remove "standalone" bins of data 
If a single 10 kb bin of data is flanked on both sides by a droplist, drop that bin in that sample and in all bioreps:
```bash
# For EdU-S/G1
## Create total droplist regions bedfile, includes the "standalone" regions:
bedtools complement -i B73_EdU_SG1_ratio_avgG1_br1_10kb_dropNonDropZero.bedgraph -g B73_chrOnly_size.txt | \
bedtools merge -d 15000 -i - > B73_EdU_dropStandalone.bed

## Remove the standalone regions, use a file that contains B73_final_droplist.bed + non-droplist zeros + standalone bins
for x in *_EdU_dropNonDropZero.bedgraph ;
do
	bedtools subtract -a $x -b B73_EdU_dropStandalone.bed > $(basename $x _dropNonDropZero.bedgraph)_final.bedgraph ;
done



# For S/G1
## Create total droplist regions bedfile, includes the "standalone" regions:
bedtools complement -i B73_SG1_ratio_avgG1_br1_10kb_dropNonDropZero.bedgraph -g ../B73_chrOnly_size.txt | \
bedtools merge -d 15000 -i - > B73_dropStandalone.bed

## Remove the standalone regions, use a file that contains B73_final_droplist.bed + non-droplist zeros + standalone bins
for x in *_dropNonDropZero.bedgraph ;
do
	bedtools subtract -a $x -b B73_dropStandalone.bed > $(basename $x _dropNonDropZero.bedgraph)_final.bedgraph ;
done
```
