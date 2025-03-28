#!/bin/bash

#SBATCH --time 2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 12
#SBATCH --mem=20GB # The memory needed depends on the size of the genome
#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id, and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --array=2-2 # Put here the row/rows from the table that need to be processed

#################
#### SET UP #####
#################

echo "=== Starting BWA Indexing Job ==="

# Set paths and dirs
pathToGenomesTable="$SRC/genomes/genomesTable.txt"
mkdir -p $SRC/genomes/bwaIndex/
pathToBwaIndex="$SRC/genomes/bwaIndex/__genome__"
pathToImages="$SRC/images"

echo "Setting up paths and directories..."
echo "Path to genomes table: $pathToGenomesTable"
echo "BWA Index Directory: $pathToBwaIndex"
echo "Singularity Image Directory: $pathToImages"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

#################
#### SCRIPT #####
#################

# Loading Singularity
echo "Loading Singularity module..."
module load singularity
echo "Singularity is loaded."

# Pull images (skip if already there)
echo "Downloading Singularity images if not already available..."
if [ ! -f "$pathToImages/bwa_0.7.18.sif" ]; then
  wget -nc -O $pathToImages/bwa_0.7.18.sif "http://datacache.galaxyproject.org/singularity/all/bwa:0.7.18--he4a0461_1"
else
  echo "BWA image already downloaded."
fi

if [ ! -f "$pathToImages/samtools.1.11.sif" ]; then
  wget -nc -O $pathToImages/samtools.1.11.sif "http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0"
else
  echo "Samtools image already downloaded"
fi

# Check if Singularity has the right permissions
ls -l /shared/projects/microc_pilot/images/
chmod 644 /shared/projects/microc_pilot/images/bwa_0.7.18.sif
chmod 644 /shared/projects/microc_pilot/images/samtools.1.11.sif
ls -l /shared/projects/microc_pilot/images/

# Function to execute singularity command for BWA
function bwa() {
  singularity exec $pathToImages/bwa_0.7.18.sif bwa $*
}

# Function to execute singularity command for samtools
function samtools() {
  singularity exec $pathToImages/samtools.1.11.sif samtools $*
}


# # Check BWA
# echo "Checking BWA..."
# bwa_output=$(singularity exec /shared/projects/microc_pilot/images/bwa_0.7.18.sif bwa 2>&1)
# if [[ $? -eq 0 ]]; then
#   echo "BWA is working fine!"
#   echo "$bwa_output"
# else
#   echo "BWA command failed!"
#   echo "$bwa_output"
#   exit 1
# fi

# # Debugging step to confirm execution
# echo "Passed BWA check, proceeding with samtools check..."

# Check samtools
echo "Checking samtools..."
samtools --version



set -e

# Proceed with genome processing
echo "Fetching genome and fasta..."
genome=$(cat ${pathToGenomesTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
filePathForFasta=$(cat ${pathToGenomesTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
echo "Genome: $genome"
echo "FASTA Path: $filePathForFasta"

# Adapt pathToBwaIndex to the name of the genome:
pathToBwaIndex=${pathToBwaIndex/__genome__/${genome}}

echo "Indexing genome: $genome at $filePathForFasta"

# Index genome
# Check if the BWA index already exists
if [ ! -e ${pathToBwaIndex}.bwt ]; then
    # Index genome if index does not exist
    if [ ! -e ${filePathForFasta}.fai ]; then
        echo "Indexing genome with samtools faidx..."
        samtools faidx "$filePathForFasta"
    fi
    echo "Running bwa index..."
    bwa index -p "$pathToBwaIndex" "$filePathForFasta"
else
    echo "BWA index seems to already exist. If you want to regenerate it, please remove it before running the job."
fi
