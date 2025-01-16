#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 3:00:00 # This depends on the size of the fasta
#SBATCH --array=1-2 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name bwa_index # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /cecilia # This directory must exists, this is where will be the error and out files






# Set paths
# 
pathToGenome="$HOME/genomes/hg38.fa.gz"   # put right genome
pathToFastqTable="$HOME/fastq/samples_fastq_table.txt"


# Get the genome name and fasta file from the table
sample=$(cat ${pathToFastqTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFastq1=$(cat ${pathToFastqTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
filePathForFastq2=$(cat ${pathToFastqTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')



mkdir $HOME/output
mkdir $HOME/output/$sample

sample_output_dir="$HOME/output/$sample"

# Step 1: Generate SAM file. Replace R1.fastq and R2.fastq with actual data.
bwa mem -5SP -T0 -t${CORES} $pathToGenome $filePathForFastq1 $filePathForFastq2 > "$sample_output_dir/${sample}.aligned.sam"

# Step 2: Record valid ligation events.
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${CORES} --nproc-out ${CORES} --chroms-path "$pathToGenome.genome" "$sample_output_dir/${sample}.aligned.sam" > "$sample_output_dir/${sample}.parsed.pairsam"

# Step 3: Sort the parsed.pairsam.
pairtools sort --nproc ${CORES} --tmpdir=/tmp "$sample_output_dir/${sample}.parsed.pairsam" > "$sample_output_dir/${sample}.sorted.pairsam"

# Step 4: Remove PCR duplicates
pairtools dedup --nproc-in ${CORES} --nproc-out ${CORES} --mark-dups --output-stats "$sample_output_dir/${sample}stats.txt" --output "$sample_output_dir/${sample}.dedup.pairsam" "$sample_output_dir/${sample}.sorted.pairsam"

# Step 5: Generate .pair and BAM files
pairtools split --nproc-in ${CORES} --nproc-out ${CORES} --output-pairs mapped.pairs --output-sam "$sample_output_dir/${sample}.unsorted.bam" "$sample_output_dir/${sample}.dedup.pairsam"

# Step 6: Sort and index BAM file
samtools sort -@${CORES} -T temp/temp.bam -o "$sample_output_dir/${sample}.mapped.PT.bam" "$sample_output_dir/${sample}.unsorted.bam" 
samtools index "$sample_output_dir/${sample}.mapped.PT.bam"

# Step 7: Run QC script
python3 get_qc.py -p "$sample_output_dir/${sample}stats.txt" > "$sample_output_dir/${sample}.stats.txt
