#!/bin/bash

#SBATCH --job-name bwa_index 
#SBATCH --time 10:00 
#SBATCH --clusters=mesopsl1
#SBATCH --partition=hi
#SBATCH --qos=mesopsl1_hi_short 
#SBATCH --account=ijerkovic_hi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 #max in mesopsl1 is 16 so let's try
#SBATCH --mem-per-cpu=5G # The memory needed depends on the size of the genome
#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job

SRC=/travail/ijerkovic/NGS
mkdir -p $SRC/images

echo $SRC

whoami

module load singularity

pathToImages="$SRC/images"

cd $pathToImages
singularity exec bwa_0.7.18.sif bwa
