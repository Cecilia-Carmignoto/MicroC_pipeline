#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 5G # The memory needed depends on the size of the genome
#SBATCH --cpus-per-task 1 # This allows to speed the indexing
#SBATCH --time 3:00:00 # This depends on the size of the fasta
#SBATCH --array 1-26 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name runMicroC # Job name that appear in squeue as well as in output and error text files

# dirs
dirPathForFastq="${microcPilot}/fastq/"
filePathForTable="${dirPathForFastq}/samplesFastqTable.txt"
dirPathWithResults="$microcPilot/outputs/"
dirPathDiagnosis="$microcPilot/outputs/diagnose/"

# Get the sample name and fastq file from the table
sample=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')

wc -l relFilePathFastqR1 | awk '{print $1 / 4}' > $dirPathDiagnosis/${sample}_total_reads_R1.txt

wc -l relFilePathFastqR2 | awk '{print $1 / 4}'> $dirPathDiagnosis/${sample}_total_reads_R2.txt


wc -l file.fastq | awk '{print $1 / 4}'