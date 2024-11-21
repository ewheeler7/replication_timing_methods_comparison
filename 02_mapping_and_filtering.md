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
  bowtie2 -p 48 --very-sensitive \
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
  samtools view -@ 48 -bf 0x2 $x | apptainer exec samtools_v1.9-4-deb_cv1.sif samtools view -@ 48 -q 6 -b - | samtools sort -@ 48 - > $(basename $x _rmDup.bam)_filtered.bam ;
done
```

## Remove reads from scaffold and low coverage droplist

Create bedfile of KEEP regions, ie the chromosomes minus droplist
**get low coverage droplist from mappability**
```bash
bedtools complement -i B73_noCovPercent_gt60_10kb.bed -g B73_chrOnly_size.txt > B73_keep_regions.bed
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

## Calculate number of reads at each step

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

