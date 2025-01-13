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
# Step 1 Generate sam file. NOTE: replace R1.fastq R2.fastq with the actual data. Replate -t8 with the appropriate number of threads
bwa mem -5SP -T0 -t8 references/mm39/mm39.fa.gz data/R1.fastq R2.fastq > aligned.sam

## Step 2. Record valid ligation events. NOTE: replace --nproc-in --n-proc-out with t
he appropriate n of cores
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 5 --nproc-out 5 --chroms-path references/mm39.genome aligned.sam > parsed.pairsam

## Step 3. Sort the parsed.pairsam. NOTE: replace --nproc with the appropriate n of cores
pairtools sort --nproc 5 --tmpdir=/tmp parsed.pairsam > sorted.pairsam

## Step 4. Remove PCR duplicates
pairtools dedup --nproc-in 5 --nproc-out 5 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam

# Step 5. Generate .pair and bam files
pairtools split --nproc-in 5 --nproc-out 5 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam

# Step 6
samtools sort -@8 -T temp/temp.bam -o mapped.PT.bam unsorted.bam

samtools index mapped.PT.bam

python3 .get_qc.py -p stats.txt


