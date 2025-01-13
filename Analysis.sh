#!bin/bash

#Download reference genome
mkdir references
cd references/
wget -O mm39.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
#ungunzip
gunzip mm39.fa.gz
#bgzip (samtools needs bgzipped to index)
bgzip mm39.fa.gz
#index the genome:
samtools faidx mm39.fa.gz
#create chrmosome sizes file
cut -f1,2 mm39.fa.gz.fai > mm39.genome 
#Bwa index
bwa index mm39.fa.gz
cd ../
# Step 1 Generate sam file. NOTE: replace R1.fastq R2.fastq with the actual data
bwa mem -5SP -T0 -t8 references/mm39/mm39.fa.gz data/R1.fastq R2.fastq > aligned.sam

## Step 2. Record valid ligation events 
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 5 --nproc-out 5 --chroms-path references/mm39.genome aligned.sam > parsed.pairsam

## Step 3. Sort the parsed.pairsam
pairtools sort --nproc 5 --tmpdir=/tmp parsed.pairsam > sorted.pairsam
