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


genome_path="$HOME/references/mm39.fa.gz"
images_path="$HOME/images"

# Use singularity to pull from quai: and manage dependencies
singularity pull $images_path/bwa_0.7.18 docker://quay.io/biocontainers/bwa:0.7.18--he4a0461_1 

# Define function to be able to call bwa 
function bwa() {
exec $images_path/bwa_0.7.18 bwa $*
}




# Instead of singularity pull from quay i can use: wget ''http://datacache.galaxyproject.org/singularity/all/bwa:0.7.18--he4a0461_1'' bwa $*




# Index 
samtools faidx mm39.fa.gz 
cut -f1,2 mm39.fa.gz.fai > mm39.genome
bwa index mm39.fa.gz

# keep only chr1-22
grep -E '^chr([1-9]|1[0-9]|2[0-2])\t' mm39.genome > mm39.genome


#_______________________

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')

R1_fastq=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')

R2_fastq=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')

bwa mem -5SP -T0 -t${CORES} "mm39.fa.gz" R1_fastq R2_fastq > aligned.sam