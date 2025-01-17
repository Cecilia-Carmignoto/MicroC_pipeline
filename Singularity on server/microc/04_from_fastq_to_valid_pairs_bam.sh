#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 3:00:00 # This depends on the size of the fasta
#SBATCH --array=1-10 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name bwa_index # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /cecilia # This directory must exists, this is where will be the error and out files

# Set paths
pathToGenome="$SRC/genomes/hg38.fa.gz"   # put right genome
pathToFastqTable="$SRC/fastq/samplesFastqTable.txt"
pathToImages="$SRC/images"

# IMAGES
wget -nc -O $pathToImages/pairtools.0.3.0 ''http://datacache.galaxyproject.org/singularity/all/pairtools:0.3.0--py37h4eba2af_0''
function samtools() {
singularity exec $pathToImages/pairtools.0.3.0 pairtools $*
}

wget -nc -O $pathToImages/tabulateAndArgparse ''https://depot.galaxyproject.org/singularity/tabulate:0.7.5--py36_0''
function tabulateAndArgparse() {
singularity exec $pathToImages/tabulateAndArgparse $*
}


# Micro-C
RUN git clone https://github.com/dovetail-genomics/Micro-C 


#################
#### SCRIPT #####
#################

# Check installations 
bwa --version
if [ $? -ne 0 ]
then
  echo "Bwa is not installed but required. Please install it"
  exit 1
fi

samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it"
  exit 1
fi

pairtools --version
if [ $? -ne 0 ]
then
  echo "pairtools is not installed but required. Please install it"
  exit 1
fi

tabulate --version
if [ $? -ne 0 ]
then
  echo "tabulate is not installed but required. Please install it"
  exit 1
fi

argparse --version
if [ $? -ne 0 ]
then
  echo "argparse is not installed but required. Please install it"
  exit 1
fi

# Get the genome name and fasta file from the table
sample=$(cat ${pathToFastqTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
pathToFastq1=$(cat ${pathToFastqTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
pathToFastq2=$(cat ${pathToFastqTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')

mkdir $SRC/output/            # can I put these two lines on the top 
mkdir $SRC/output/$sample/    # or I need to define $sample before being able to execute these lines

sample_output_dir="$SRC/output/$sample/"

# Step 1: Generate SAM file. Replace R1.fastq and R2.fastq with actual data.
bwa mem -5SP -T0 -t${CORES} $pathToGenome $pathToFastq1 $pathToFastq2 > "${sample_output_dir}aligned.sam"

# Step 2: Record valid ligation events.
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in ${CORES} --nproc-out ${CORES} --chroms-path "${pathToGenome}.genome" "${sample_output_dir}aligned.sam" > "${sample_output_dir}parsed.pairsam"

# Step 3: Sort the parsed.pairsam.
pairtools sort --nproc ${CORES} --tmpdir=/tmp "${sample_output_dir}parsed.pairsam" > "${sample_output_dir}sorted.pairsam"

# Step 4: Remove PCR duplicates
pairtools dedup --nproc-in ${CORES} --nproc-out ${CORES} --mark-dups --output-stats "${sample_output_dir}stats.txt" --output "${sample_output_dir}dedup.pairsam" "${sample_output_dßr}sorted.pairsam"

# Step 5: Generate .pair and BAM files
pairtools split --nproc-in ${CORES} --nproc-out ${CORES} --output-pairs mapped.pairs --output-sam "${sample_output_dir}unsorted.bam" "${sample_output_dir}dedup.pairsam"

# # Step 6: Sort and index BAM file
# samtools sort -@${CORES} -T temp/temp.bam -o "${sample_output_dir}mapped.PT.bam" "${sample_output_dir}unsorted.bam" 
# samtools index "$sample_output_dir/${sample}.mapped.PT.bam"

# Step 7: Run QC script
tabulateAndArgparse get_qc.py -p "${sample_output_dir}stats.txt" 

# Save the stats in a common file for all samples
mkdir $SRC/output/stats_all/

touch $SRC/output/stats_all/stats_all.txt

echo -e "$sample" >> $SRC/output/stats_all/stats_all.txt
cat "${sample_output_dir}stats.txt" >> "$SRC/output/stats_all/