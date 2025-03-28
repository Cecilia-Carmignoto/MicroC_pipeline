#!/bin/bash
#SBATCH -o slurm-%x-%A_%2a.out  # Standard output template
#SBATCH -e slurm-%x-%A_%2a.err  # Standard error template
#SBATCH --nodes 1               # Use 1 node
#SBATCH --ntasks 1              # Sequential tasks
#SBATCH --mem 50G               # Set memory based on genome size
#SBATCH --cpus-per-task 24      # CPU cores to speed up tasks
#SBATCH --time 3:00:00          # Adjust time based on genome size
#SBATCH --array 2-2             # Specify the genome row to process
#SBATCH --job-name microc_pilot    # Job name for SLURM
#SBATCH --chdir /shared/home/ccarmignoto99        # Set working directory for job

# Set environment variable (make sure this path is correct)
export SRC="/shared/projects/microc_pilot"

cat /etc/resolv.conf


# Load Singularity if it's not already loaded (ensure the correct version is loaded)
module load singularity

# Check if Singularity is loaded
if ! module list 2>&1 | grep -q 'singularity'; then
    echo "Singularity is not loaded. Loading Singularity..."
    module load singularity
else
    echo "Singularity is already loaded."
fi

# Set up the paths
pathToGenomesTable="$SRC/genomes/genomesTable.txt"
genomeLine=2  # Specify which genome to download
mkdir -p $SRC/genomes/fasta/
mkdir -p $SRC/images
pathToImages="$SRC/images"

wget -O $SRC/images/sample.jpg https://upload.wikimedia.org/wikipedia/commons/3/3a/Cat03.jpg

# Pull Images (only if not already downloaded)
if [ ! -f $pathToImages/samtools.1.11.sif ]; then
    wget -nc -O $pathToImages/samtools.1.11.sif "https://depot.galaxyproject.org/singularity/samtools:1.11--h6270b1f_0"
fi

# Functions to run samtools and bgzip using Singularity
function samtools() {
    singularity exec $pathToImages/samtools.1.11.sif samtools "$@"
}

function bgzip() {
    singularity exec $pathToImages/samtools.1.11.sif bgzip "$@"
}

samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it"
  exit 1
fi

echo "Testing internet connection..."
ping -c 4 google.com

# If ping is successful, continue with the rest of the script
if [ $? -eq 0 ]; then
    echo "Internet connection successful. Continuing with genome download and processing..."
else
    echo "Internet connection failed. Exiting script."
    exit 1
fi

# Extract genome information from the genomesTable.txt
genomeName=$(awk -v i=$genomeLine 'NR==i{print $1}' ${pathToGenomesTable})
filePathForFasta=$(awk -v i=$genomeLine 'NR==i{print $2}' ${pathToGenomesTable})
genomeURL=$(awk -v i=$genomeLine 'NR==i{print $3}' ${pathToGenomesTable})

# Download the genome file
cd $SRC/genomes/fasta
wget -nc -O $filePathForFasta $genomeURL
gunzip -c $filePathForFasta > genomeUnzipped.fa
bgzip genomeUnzipped.fa
mv genomeUnzipped.fa.gz $filePathForFasta

# Index the genome if necessary
if [ ! -e ${filePathForFasta}.fai ]; then
    samtools faidx "$filePathForFasta"
fi

# Filter to get only the numbered chromosomes
cut -f1,2 "${filePathForFasta}.fai" > "${filePathForFasta}.genome"
grep -P '^chr([1-9]|1[0-9]|2[0-2])\t' "${filePathForFasta}.genome" > "$SRC/genomes/fasta/${genomeName}_chrNumbered.genome"


