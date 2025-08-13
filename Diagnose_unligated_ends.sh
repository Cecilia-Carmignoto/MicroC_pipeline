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

mkdir -p $dirPathDiagnosis

# Get the sample name and fastq file from the table
sample=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')


zgrep -i -c AGGTTCGTCCATCGATCGATGGACGAACCT $dirPathForFastq/$relFilePathFastqR1 > $dirPathDiagnosis/${sample}_full_seq_R1.txt
zgrep -i -c AGGTTCGTCCATC $dirPathForFastq/$relFilePathFastqR1 > $dirPathDiagnosis/${sample}_start_seq_R1.txt
zgrep -i -c GATGGACGAACCT $dirPathForFastq/$relFilePathFastqR1 > $dirPathDiagnosis/${sample}_end_seq_R1.txt
zgrep -i -c CATCGATCGAT $dirPathForFastq/$relFilePathFastqR2 > $dirPathDiagnosis/${sample}_internal_seq_R2.txt


zgrep -i -c AGGTTCGTCCATCGATCGATGGACGAACCT $dirPathForFastq/$relFilePathFastqR1 > $dirPathDiagnosis/${sample}_full_seq_R2.txt
zgrep -i -c AGGTTCGTCCATC $dirPathForFastq/$relFilePathFastqR2 > $dirPathDiagnosis/${sample}_start_seq_R2.txt
zgrep -i -c GATGGACGAACCT $dirPathForFastq/$relFilePathFastqR2 > $dirPathDiagnosis/${sample}_end_seq_R2.txt
zgrep -i -c CATCGATCGAT $dirPathForFastq/$relFilePathFastqR2 > $dirPathDiagnosis/${sample}_internal_seq_R2.txt

