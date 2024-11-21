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


