#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G 
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 23:00:00 # This depends on the size of the fasta
#SBATCH --job-name runMicroC # Job name that appear in squeue as well as in output and error text files

pathToImages="$SRC/images"

wget -nc -O $pathToImages/pairtools.0.3.0 "http://datacache.galaxyproject.org/singularity/all/pairtools:0.3.0--py37h4eba2af_0"
function pairtools() {
  singularity exec $pathToImages/pairtools.0.3.0 pairtools $*
}

pairtools --version
if [ $? -ne 0 ]
then
  echo "pairtools is not installed but required. Please install it"
  exit 1
fi

sample=/shared/projects/microc_pilot/microc/pilot2/outputs/mc_120h_R1_unbio/mc_120h_R1_unbio
pairtools sort --nproc ${SLURM_CPUS_PER_TASK} --tmpdir=$(mktemp -d) "${sample}.parsed.pairsam" > "${sample}.sorted.pairsam"
