#!/bin/bash
# be sure to put the fastq data in the directory data
# Set global variables
GENOME="hg38"
GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
CORES=5
DATA="data"

# Download reference genome
 wget -O "${GENOME}.fa.gz" "${GENOME_URL}"

# Uncompress the genome file
gunzip "${GENOME}.fa.gz"

# Compress with bgzip (samtools needs bgzipped file to index)
bgzip "${GENOME}.fa"

# Index the genome using samtools
samtools faidx "${GENOME}.fa.gz"

# Create chromosome sizes file
cut -f1,2 "${GENOME}.fa.gz.fai" > "${GENOME}.genome"

# Keep only chr from 1 to 22
grep -E '^chr([1-9]|1[0-9]|2[0-2])\t' hg38.genome > hg38.genome

# BWA index
bwa index "${GENOME}.fa.gz"

# Step 1: Generate SAM file. Replace R1.fastq and R2.fastq with actual data.
bwa mem -5SP -T0 -t${CORES} "${GENOME}.fa.gz" "${DATA}/ MicroC_2M_R1.fastq" "${DATA}/ MicroC_2M_R2.fastq" > aligned.sam

# Step 2: Record valid ligation events.
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${CORES} --nproc-out ${CORES} --chroms-path "${GENOME}.genome" aligned.sam > parsed.pairsam

# Step 3: Sort the parsed.pairsam.
pairtools sort --nproc ${CORES} --tmpdir=/tmp parsed.pairsam > sorted.pairsam

# Step 4: Remove PCR duplicates
pairtools dedup --nproc-in ${CORES} --nproc-out ${CORES} --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam

# Step 5: Generate .pair and BAM files
pairtools split --nproc-in ${CORES} --nproc-out ${CORES} --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam

# Step 6: Sort and index BAM file
samtools sort -@${CORES} -T temp/temp.bam -o mapped.PT.bam unsorted.bam
samtools index mapped.PT.bam

# Step 7: Run QC script
mkdir output
python3 get_qc.py -p stats.txt > output/stats.txt
