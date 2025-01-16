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

mkdir -p $HOME/images
images_path="$HOME/images"
pathToGenomesTable="$HOME/genomes_table.txt"

#################
#### SET UP #####
#################

# Set paths
# The table genomes_table.txt has to be already generated. (See README.md)
# first column is the genome name
# second column is the absolute path for fasta
pathToGenomesTable="$HOME/genomes_table.txt"
pathToBwaIndex="$HOME/genomes/__genome__"

# Pull tools/softwares images
# Use singularity to pull from quai: and manage dependencies:
# singularity pull $images_path/bwa_0.7.18 docker://quay.io/biocontainers/bwa:0.7.18--he4a0461_1 
# Or use wget
# Then define function to be able to call the executable 
 
# Note: is it better to use $* or $@

# bwa
wget -O $images_path/bwa_0.7.18.sif 'http://datacache.galaxyproject.org/singularity/all/bwa:0.7.18--he4a0461_1'
function bwa() {
singularity exec $images_path/bwa_0.7.18.sif bwa $*
}

# samtools
wget -O $images_path/samtools.1.11.sif ''http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0''
function bwa() {
singularity exec $images_path/samtools.1.11.sif samtools $*
}

#################
#### SCRIPT #####
#################

# Check set up
bwa --version
if [ $? -ne 0 ]
then
  echo "Bwa is not installed but required. Please install it"
  exit 1
fi

# Get the genome name and fasta file from the table
genome=$(cat ${pathToGenomesTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFasta=$(cat ${pathToGenomesTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

# Adapt pathToBwaIndex to the name of the genome:
pathToBwaIndex=${pathToBwaIndex/__genome__/${genome}}

if [ ! -e ${pathToBwaIndex}.bwt ]; then

    samtools faidx $pathToBwaIndex
    cut -f1,2 "${pathToBwaIndex}.fai" > "${pathToBwaIndex}.genome"
    grep -E '^chr([1-9]|1[0-9]|2[0-2])\t' "${pathToBwaIndex}.genome" > "${pathToBwaIndex}.genome"
    bwa index $pathToBwaIndex
else
    echo "bwa index seems to already exists. If you want to regenerate it. Please remove it before running the job."
fi




#_______________________

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')

R1_fastq=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')

R2_fastq=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')

bwa mem -5SP -T0 -t${CORES} "mm39.fa.gz" R1_fastq R2_fastq > aligned.sam