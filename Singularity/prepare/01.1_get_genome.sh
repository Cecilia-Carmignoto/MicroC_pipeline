#!/bin/bash
set -e # if a line in the middle fails the script fails

# TO RUN ONLY ONCE
# Download of the reference genome
# Choose which reference genome to use by setting the variable genomeLine. (check genomesTable.txt and chose the line)

#################
#### SET UP #####
#################

# Define paths
# CHECK: The table genomesTable.txt has to be already generated. (See README.md)
# first column is the genome name
# second column is the absolute path for fasta
pathToGenomesTable="$SRC/genomes/genomesTable.txt"
genomeLine=2            # Set the line number of genomesTable.txt of the genome to download (line 1 hg38, line 2 mm39)
mkdir -p $SRC/genomes/fasta
mkdir -p $SRC/images
pathToImages="$SRC/images"

#################
#### SCRIPT #####
#################

# Pull Images
# bgzip is in the image of samtools
# samtools
wget -nc -O $pathToImages/samtools.1.11.sif "http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0"
function samtools() {
  singularity exec $pathToImages/samtools.1.11.sif samtools $*
}
# bgzip from samtools
function bgzip() {
  singularity exec $pathToImages/samtools.1.11.sif bgzip $*
}

# Get the genome name and fasta file from the genomesTable.txt
genomeName=$(cat ${pathToGenomesTable} | awk -v i=$genomeLine 'NR==i{print $1}')
filePathForFasta=$(cat ${pathToGenomesTable} | awk -v i=${genomeLine} 'NR==i{print $2}')
genomeURL=$(cat ${pathToGenomesTable} | awk -v i=$genomeLine 'NR==i{print $3}')

# Download genome
cd $SRC/genomes/fasta
wget -nc -O $filePathForFasta $genomeURL
gunzip -c $filePathForFasta > genomeUnzipped.fa
bgzip genomeUnzipped.fa # is it okay to run this in the front end? yes
mv genomeUnzipped.fa.gz $filePathForFasta
