# Run Repliscan program on Repli-seq data
The mapped, filtered, and droplist removed mapping files from all the bioreps from the G1, Early, Mid, and Late sorting gates was
processed through Repliscan (Zynda, 2017; https://github.com/zyndagj/repliscan).

Container for Repliscan is repliscan_v1.1_1.1.sif from https://hub.docker.com/u/lconcia.



## Pipeline
- Use bam files from previous step with high and low coverage droplists removed
- Run repliscan 
- Re-remove droplist, because smoothing internal to repliscan treats droplist at 0 and so will smooth over it
- Convert gff3 file to bed (for cleaner viewing on IVG)



## Run repliscan
Make configuration file needed for Repliscan:
```bash
echo "G1	BR2_B73_G1.bam	BR3_B73_G1.bam	BR4_B73_G1.bam
ES	BR2_B73_E.bam	BR3_B73_E.bam	BR4_B73_E.bam
MS	BR2_B73_M.bam	BR3_B73_M.bam	BR4_B73_M.bam
LS	BR2_B73_L.bam	BR3_B73_L.bam	BR4_B73_L.bam " > input.txt

***these names need to be edited!!!!***
```

Repliscan command, using 10kb window size_________:
```bash
repliscan_removing_blacklist_test.py -l 2 -w 10000 -a mean -t value -S genome -v 1 -c proportion --plot -r Zm-B73-REFERENCE-NAM-5.0.fa input.txt >> repliscan_removing_blacklist_test_10kb.log 2>> repliscan_removing_blacklist_test_10kb.err
```


## Re-remove droplist, because smoothing internal to repliscan treats droplist at 0 and so will smooth over it
- Also remove any bin that has an RT value of 0
- Also remove scaffold placeholders

### Find all regions with RT=0:
cat *_2.smooth.bedgraph | grep "chr" - | awk '{ if ($4==0) print $0}' - | \
sort -k 1V,1 -k 2n,2 | bedtools merge -i - > repliscan_zero_to_remove.bed

### Re-remove high coverage droplist, low mappability droplist, and all zeros in the data from bedgraphs:
for x in *_2.smooth.bedgraph; 
do
	grep "chr" $x | \
	bedtools subtract -a - -b B73_pooled_high_droplist_repliscan.bed | \
	bedtools subtract -a - -b B73_noCovPercent_gt60_10kb.bed | \
	bedtools subtract -a - -b repliscan_zero_to_remove.bed > $(basename $x .bedgraph)_final.bedgraph;
done



## Re-remove high coverage droplist, low mappability droplist, and all zeros in the data from gff3 output file:
- When removing regions in a gff3 "gene" it will think there is an intron, so remove "ID=gene" from last column
- Note, there are regions where there is singal but no Segmentation call is made, this is ok
```bash
for x in ./*gff3; 
do
	head -n 2 $x > header.txt; 
	grep -v "scaf" $x | \
	bedtools subtract -a - -b B73_pooled_high_droplist_repliscan.bed | \
	bedtools subtract -a - -b B73_noCovPercent_gt60_10kb.bed | \
	bedtools subtract -a - -b repliscan_zero_to_remove.bed | \
	sed "s/ID=gene\([\w]*\)/\1/g" - | sed -e 's/S//g' > $(basename $x .gff3)_tmp.gff3 ; #need to remove ID= or gff3 thinks its a gene, could be formatted better
	cat header.txt $(basename $x .gff3)_tmp.gff3 > $(basename $x .gff3)_final.gff3 ; 
	rm header.txt $(basename $x .gff3)_tmp.gff3 ;
done
```

Get rid of repliscan segmentation name, for cleaner IGV, and make the EML and EL segments grey: 
```bash
grep "=ES;" ratio_segmentation.gff3 | cut -f 1-8 - | awk '{ print $0"\t""color=#2250F1" }' - > repliseq_E_tmp.gff3
grep "=ESMS;" ratio_segmentation.gff3 | cut -f 1-8 - |  awk '{ print $0"\t""color=#28C5CC" }' - > repliseq_EM_tmp.gff3
grep "=MS;" ratio_segmentation.gff3 | cut -f 1-8 - |  awk '{ print $0"\t""color=#1A8A12" }' - > repliseq_M_tmp.gff3
grep "=MSLS;" ratio_segmentation.gff3 | cut -f 1-8 - |  awk '{ print $0"\t""color=#FFFD33" }' - > repliseq_ML_tmp.gff3
grep "=LS;" ratio_segmentation.gff3 | cut -f 1-8 - |  awk '{ print $0"\t""color=#FB0018" }' - > repliseq_L_tmp.gff3
#grep "=ESMSLS;" ratio_segmentation.gff3 | cut -f 1-8 - |  awk '{ print $0"\t""color=#808080" }' - > repliseq_EML_tmp.gff3
#grep "=ESLS;" ratio_segmentation.gff3 | cut -f 1-8 - |  awk '{ print $0"\t""color=#808080" }' - > repliseq_EL_tmp.gff3
cat repliseq_E_tmp.gff3 repliseq_EM_tmp.gff3 repliseq_M_tmp.gff3 repliseq_ML_tmp.gff3 repliseq_L_tmp.gff3 | \
sort -k 1V,1 -k 2n,2 - > new_repliseq_ratio_seg.gff3
rm *tmp.gff3
```

