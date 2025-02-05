# Read pre-processing

## Download non-EdU sequencing files, BioProject PRJNA1137362
### Get fastq files from NCBI SRA for primary roots
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_primary_G1_br1 | SRR29881197 | BR1_Primary_G1_S26_L002_R1_001.fastq.gz | BR1_Primary_G1_S26_L002_R2_001.fastq.gz
B73_primary_S_br1 | SRR29881202 | BR1_Primary_S_S27_L002_R1_001.fastq.gz | BR1_Primary_S_S27_L002_R2_001.fastq.gz
B73_primary_G1_br2 | SRR29881196 | BR2_Primary_G1_S30_L002_R1_001.fastq.gz | BR2_Primary_G1_S30_L002_R2_001.fastq.gz
B73_primary_S_br2 | SRR29881201 | BR2_Primary_S_S31_L002_R1_001.fastq.gz | BR2_Primary_S_S31_L002_R2_001.fastq.gz
B73_primary_G1_br3 | SRR29881195 | BR3_Primary_G1_S33_L002_R1_001.fastq.gz | BR3_Primary_G1_S33_L002_R2_001.fastq.gz
B73_primary_S_br23 | SRR29881198 | BR3_Primary_S_S34_L002_R1_001.fastq.gz | BR3_Primary_S_S34_L002_R2_001.fastq.gz

### Get fastq files from NCBI SRA for seminal roots
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_seminal_G1_br1 | SRR29881191 | BR1_Seminal_G1_S28_L002_R1_001.fastq.gz | BR1_Seminal_G1_S28_L002_R2_001.fastq.gz
B73_seminal_S_br1 | SRR29881194 | BR1_Seminal_S_S29_L002_R1_001.fastq.gz | BR1_Seminal_S_S29_L002_R2_001.fastq.gz
B73_seminal_G1_br2 | SRR29881200 | BR2_Seminal_G1_S32_L002_R1_001.fastq.gz | BR2_Seminal_G1_S32_L002_R2_001.fastq.gz
B73_seminal_S_br2 | SRR29881193 | BR2_Seminal_S_S25_L002_R1_001.fastq.gz | BR2_Seminal_S_S25_L002_R2_001.fastq.gz
B73_seminal_G1_br3 | SRR29881191 | BR3_Seminal_G1_S35_L002_R1_001.fastq.gz | BR3_Seminal_G1_S35_L002_R2_001.fastq.gz
B73_seminal_S_br3 | SRR29881192| BR3_Seminal_S_S36_L002_R1_001.fastq.gz | BR3_Seminal_S_S36_L002_R2_001.fastq.gz


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
cat B73_Primary_G1_br1_R1.fastq.gz B73_Seminal_G1_br1_R1.fastq.gz > B73_G1_br1_r1.fastq.gz
cat B73_Primary_G1_br1_R2.fastq.gz B73_Seminal_G1_br1_R2.fastq.gz > B73_G1_br1_r2.fastq.gz
cat B73_Primary_S_br1_R1.fastq.gz B73_Seminal_S_br1_R1.fastq.gz > B73_S_br1_r1.fastq.gz
cat B73_Primary_S_br1_R2.fastq.gz B73_Seminal_S_br1_R2.fastq.gz > B73_S_br1_r2.fastq.gz

cat B73_Primary_G1_br2_R1.fastq.gz B73_Seminal_G1_br2_R1.fastq.gz > B73_G1_br2_r1.fastq.gz
cat B73_Primary_G1_br2_R2.fastq.gz B73_Seminal_G1_br2_R2.fastq.gz > B73_G1_br2_r2.fastq.gz
cat B73_Primary_S_br2_R1.fastq.gz B73_Seminal_S_br2_R1.fastq.gz > B73_S_br2_r1.fastq.gz
cat B73_Primary_S_br2_R2.fastq.gz B73_Seminal_S_br2_R2.fastq.gz > B73_S_br2_r2.fastq.gz

cat B73_Primary_G1_br3_R1.fastq.gz B73_Seminal_G1_br3_R1.fastq.gz > B73_G1_br3_r1.fastq.gz
cat B73_Primary_G1_br3_R2.fastq.gz B73_Seminal_G1_br3_R2.fastq.gz > B73_G1_br3_r2.fastq.gz
cat B73_Primary_S_br3_R1.fastq.gz B73_Seminal_S_br3_R1.fastq.gz > B73_S_br3_r1.fastq.gz
cat B73_Primary_S_br3_R2.fastq.gz B73_Seminal_S_br3_R2.fastq.gz > B73_S_br3_r2.fastq.gz
```

## Download EdU seqeuncing files, BioProject PRJNA1136904
### Get fastq files from NCBI SRA for first and second round of sequencing
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
B73_EdU_G1_br1 | SRR29860876 | BR1_G1_S1_L001_R1_001.fastq.gz | BR1_G1_S1_L001_R2_001.fastq.gz
B73_EdU_G1_br2 | SRR29860875 | BR2_G1_S3_L001_R1_001.fastq.gz | BR2_G1_S3_L001_R2_001.fastq.gz
B73_EdU_G1_br3 | SRR29860874 | BR3_G1_S5_L001_R1_001.fastq.gz | BR3_G1_S5_L001_R2_001.fastq.gz
B73_EdU_S_br1 | SRR29860879 | BR1_S_S2_L001_R1_001.fastq.gz | BR1_S_S2_L001_R2_001.fastq.gz
B73_EdU_S_br2 | SRR29860878 | BR2_S_S4_L001_R1_001.fastq.gz | BR2_S_S4_L001_R2_001.fastq.gz
B73_EdU_S_br3 | SRR29860877 | BR3_S_S6_L001_R1_001.fastq.gz | BR3_S_S6_L001_R2_001.fastq.gz



### Merge first and second round of sequencing:
```bash
cat B73_EdU_SG1_BR1_G1_S1_L003_R1_001.fastq.gz BR1_G1_S1_L001_R1_001.fastq.gz > B73_EdU_G1_br1_r1.fastq.gz
cat B73_EdU_SG1_BR1_G1_S1_L003_R2_001.fastq.gz BR1_G1_S1_L001_R2_001.fastq.gz > B73_EdU_G1_br1_r2.fastq.gz
cat B73_EdU_SG1_BR1_S_S2_L003_R1_001.fastq.gz BR1_S_S2_L001_R1_001.fastq.gz > B73_EdU_S_br1_r1.fastq.gz
cat B73_EdU_SG1_BR1_S_S2_L003_R2_001.fastq.gz BR1_S_S2_L001_R2_001.fastq.gz > B73_EdU_S_br1_r2.fastq.gz
cat B73_EdU_SG1_BR2_G1_S3_L003_R1_001.fastq.gz BR2_G1_S3_L001_R1_001.fastq.gz > B73_EdU_G1_br2_r1.fastq.gz
cat B73_EdU_SG1_BR2_G1_S3_L003_R2_001.fastq.gz BR2_G1_S3_L001_R2_001.fastq.gz > B73_EdU_G1_br2_r2.fastq.gz
cat B73_EdU_SG1_BR2_S_S4_L003_R1_001.fastq.gz BR2_S_S4_L001_R1_001.fastq.gz > B73_EdU_S_br2_r1.fastq.gz
cat B73_EdU_SG1_BR2_S_S4_L003_R2_001.fastq.gz BR2_S_S4_L001_R2_001.fastq.gz > B73_EdU_S_br2_r2.fastq.gz
cat B73_EdU_SG1_BR3_G1_S5_L003_R1_001.fastq.gz BR3_G1_S5_L001_R1_001.fastq.gz > B73_EdU_G1_br3_r1.fastq.gz
cat B73_EdU_SG1_BR3_G1_S5_L003_R2_001.fastq.gz BR3_G1_S5_L001_R2_001.fastq.gz > B73_EdU_G1_br3_r2.fastq.gz
cat B73_EdU_SG1_BR3_S_S6_L003_R1_001.fastq.gz BR3_S_S6_L001_R1_001.fastq.gz > B73_EdU_S_br3_r1.fastq.gz
cat B73_EdU_SG1_BR3_S_S6_L003_R2_001.fastq.gz BR3_S_S6_L001_R2_001.fastq.gz > B73_EdU_S_br3_r2.fastq.gz
```

## Download Repli-seq sequencing files, BioProject PRJNA1219493
### Get fastq files from NCBI SRA 
Sample  | Accesion | Forward | Reverse
--- | --- | --- | --- 
G1_br1 | --- | BR2_B73_G1_A7_S1_L003_R1_001.fastq.gz	| BR2_B73_G1_A7_S1_L003_R2_001.fastq.gz
G1_br2 | --- | BR3_B73_G1_A8_S7_L003_R1_001.fastq.gz	| BR3_B73_G1_A8_S7_L003_R2_001.fastq.gz
G1_br3 | --- | BR4_B73_G1_A6_S13_L003_R1_001.fastq.gz	| BR4_B73_G1_A6_S13_L003_R2_001.fastq.gz
Early_br1 | --- | BR2_B73_E_C7_S3_L003_R1_001.fastq.gz	| BR2_B73_E_C7_S3_L003_R2_001.fastq.gz
Early_br2 | --- | BR3_B73_E_C8_S9_L003_R1_001.fastq.gz	| BR3_B73_E_C8_S9_L003_R2_001.fastq.gz
Early_br3 | --- | BR4_B73_E_C6_S15_L003_R1_001.fastq.gz	| BR4_B73_E_C6_S15_L003_R2_001.fastq.gz
Mid_br1 | --- | BR2_B73_M_D7_S4_L003_R1_001.fastq.gz	| BR2_B73_M_D7_S4_L003_R2_001.fastq.gz
Mid_br2 | --- | BR3_B73_M_D8_S10_L003_R1_001.fastq.gz	| BR3_B73_M_D8_S10_L003_R2_001.fastq.gz
Mid_br3 | --- | BR4_B73_M_D6_S16_L003_R1_001.fastq.gz	| BR4_B73_M_D6_S16_L003_R2_001.fastq.gz
Late_br1 | --- | BR2_B73_L_E7_S5_L003_R1_001.fastq.gz	| BR2_B73_L_E7_S5_L003_R2_001.fastq.gz
Late_br2 | --- | BR3_B73_L_E8_S11_L003_R1_001.fastq.gz	| BR3_B73_L_E8_S11_L003_R2_001.fastq.gz
Late_br3 | --- | BR4_B73_L_E6_S17_L003_R1_001.fastq.gz	| BR4_B73_L_E6_S17_L003_R2_001.fastq.gz


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

## Run Fastqc on all trimmed fastq files
```bash
for x in *.fastq.gz
do
  fastqc -t 48 $x -o . ;
done 
```
