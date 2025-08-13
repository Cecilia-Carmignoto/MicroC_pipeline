#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 100G 
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 23:00:00 # This depends on the size of the fasta
#SBATCH --array 1-1 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name runMicroC # Job name that appear in squeue as well as in output and error text files


###################################
#### TO SET FOR EACH ANALYSIS #####
###################################

### Specify the options for your analysis:

# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}
# Which genome to map on
genome=mm39
#bin size, in kb, for the .cool file
binSizeCoolMatrix=1
# Define a test region for a pgt plot (must be inside the captured region if it is a capture)
# chr7:155000000-158000000 SHH hg38
# chr2:174800000-177800000 HOXD hg38
# chr3:65103500-68603411 Shox2 CaptureC mm39
# chr2:73150000-76150000 HoxD mm39
# chr2:73779626-75669724 HoxD mm10
testRegion="chr2:73150000-76150000"
# bins size (in kb) for the plot:
bins="10 20 50"

### Specify the paths to the directories

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$microcPilot2/outputs/"
# Where fastqs are stored:
dirPathForFastq="$microcPilot2/fastq/"
# This script will use get_qc.py
# from Micro-C:
# https://raw.githubusercontent.com/dovetail-genomics/Micro-C/refs/heads/main/get_qc.py
# If it does not exists
# The python script will be put in dirPathForScripts
dirPathForScripts="$SRC/images/"
# All samples are registered into a table where
# first column is the sample name
# second column is the path of the R1 fastq relatively to dirPathForFastq
# third column is same for R2
# Alternatively second column can be SRA number but third column must be filled by anything for example also the SRA number
filePathForTable="${dirPathForFastq}/samplesFastqTable.txt"

pathToBwaIndex=$SRC/genomes/bwaIndex/$genome
filePathForFasta="$SRC/genomes/fasta/${genome}.fa.gz"
# You need to generate a 'genome' file which is made of 2 columns. First column is the chromosom name, second column is the size of the chromosome.
# You can decide to filter chromosomes at this step:
filePathForSizesForBin="$SRC/genomes/fasta/${genome}_chrNumbered.genome"

### Specify the way to deal with dependencies:
# Here we use singularity

# Images
pathToImages="$SRC/images"

wget -nc -O $pathToImages/bwa_0.7.18.sif "http://datacache.galaxyproject.org/singularity/all/bwa:0.7.18--he4a0461_1"
function bwa() {
  singularity exec $pathToImages/bwa_0.7.18.sif bwa $*
}

wget -nc -O $pathToImages/pairtools.0.3.0 "http://datacache.galaxyproject.org/singularity/all/pairtools:0.3.0--py37h4eba2af_0"
function pairtools() {
  singularity exec $pathToImages/pairtools.0.3.0 pairtools $*
}

wget -nc -O $pathToImages/samtools.1.11.sif "http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0"
function samtools() {
  singularity exec $pathToImages/samtools.1.11.sif samtools $*
}
function bgzip() {
  singularity exec $pathToImages/samtools.1.11.sif bgzip $*
}

wget -nc -O $pathToImages/tabulate:0.7.5--py36_0 "https://depot.galaxyproject.org/singularity/tabulate:0.7.5--py36_0"
function python() {
  singularity exec $pathToImages/tabulate:0.7.5--py36_0 python $*
}
function tabulate() {
  singularity exec $pathToImages/tabulate:0.7.5--py36_0 tabulate $*
}

wget -nc -O $pathToImages/cooler.0.10.3 "https://depot.galaxyproject.org/singularity/cooler:0.10.3--pyhdfd78af_0"
function cooler() {
  singularity exec $pathToImages/cooler.0.10.3 cooler $*
}
function pairix() {
  singularity exec $pathToImages/cooler.0.10.3 pairix $*
}

wget -nc -O $pathToImages/pygenometracks.3.9 "https://depot.galaxyproject.org/singularity/pygenometracks:3.9--pyhdfd78af_0"
function pgt() {
  singularity exec $pathToImages/pygenometracks.3.9 pgt $*
}
function pyGenomeTracks() {
  singularity exec $pathToImages/pygenometracks.3.9 pyGenomeTracks $*
}
# Give access to path with index with fastqs etc
export APPTAINER_BIND=$SRC

# python QC script
wget -nc -O $dirPathForScripts/get_qc.py https://raw.githubusercontent.com/dovetail-genomics/Micro-C/refs/heads/main/get_qc.py

# library complexity
wget -nc -O $pathToImages/preseq:3.2.0--hd36ca80_4.sif "https://depot.galaxyproject.org/singularity/preseq:3.2.0--hd36ca80_4"
function preseq() {
  singularity exec $pathToImages/preseq:3.2.0--hd36ca80_4.sif preseq $*
}
 
# Check installations 
v=$(bwa 2>&1)
if [[ "$v" = *"command not found" ]]
then
  echo "Bwa is not installed but required. Please install it"
  exit 1
fi
echo $v

samtools --version
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it"
  exit 1
fi

bgzip --version
if [ $? -ne 0 ]
then
  echo "bgzip is not installed but required. Please install it"
  exit 1
fi

pairtools --version
if [ $? -ne 0 ]
then
  echo "pairtools is not installed but required. Please install it"
  exit 1
fi

tabulate --version
if [ $? -ne 2 ]
then
  echo "tabulate is not installed but required. Please install it"
  exit 1
fi

# This is not working with singularity because of the quotes
# python -c "import argparse;print(argparse.__version__)"
echo "import argparse;print(argparse.__version__)" > test_argparse.py
python test_argparse.py
if [ $? -ne 0 ]
then
  echo "argparse is not installed but required. Please install it"
  exit 1
fi

cooler --version
if [ $? -ne 0 ]
then
  echo "cooler is not installed but required. Please install it"
  exit 1
fi

pairix --help
if [ $? -ne 0 ]
then
  echo "pairix is not installed but required. Please install it"
  exit 1
fi

pgt --version
if [ $? -ne 0 ]
then
  echo "pyGenomeTracks is not installed but required. Please install it"
  exit 1
fi

preseq --version
if [ $? -ne 0 ]
then
  echo "preseq is not installed but required. Please install it"
  exit 1
fi


# Get the sample name and fastq file from the table
sample=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')

# The directory is created (if not existing)
pathResults=${dirPathWithResults}/${sample}/
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo ${sample}

# The analysis part takes part within the pathResults
cd ${pathResults}

# # Generate SAM file
if [ ! -e ${sample}.sam ]; then
  if [ ! -e ${dirPathForFastq}/${relFilePathFastqR1} ]; then
    # If the fastq does not exists we assume it was an SRA ID
    mkdir -p ${dirPathForFastq}
    cd ${dirPathForFastq}
    # Write version to stdout:
    fasterq-dump --version
    if [ $? -ne 0 ]
    then
      echo "fasterq-dump is not installed and fastqFile not found so assumed it was a SRA ID.
Please install it for example in the conda environment (sra-tools>=2.11)."
      exit 1
    fi
    fasterq-dump -o ${sample}.fastq ${relFilePathFastqR1}
    if [ ! -s ${sample}_1.fastq ]; then
        echo "FASTQ R1 IS EMPTY"
        exit 1
    fi
    gzip ${sample}_1.fastq
    gzip ${sample}_2.fastq
    cd $pathResults
    relFilePathFastqR1=${sample}_1.fastq.gz
    relFilePathFastqR2=${sample}_2.fastq.gz
  fi
  if [ ! -s ${dirPathForFastq}/${relFilePathFastqR1} ]; then
    echo "FASTQ R1 IS EMPTY"
    exit 1
  fi
  # -5 is for split alignemnent, takes the alignemnte of the 5' read as primary. -S skips mate rescue, -P skips pairing. -T sets the minimus mapping quality (we want all reads to compute the stats)
  bwa mem -5SP -T0 -t${nbOfThreads} $pathToBwaIndex ${dirPathForFastq}/${relFilePathFastqR1} \
    ${dirPathForFastq}/${relFilePathFastqR2} > ${sample}.sam
else
  echo "${sample}.sam already exists"
fi

# Record valid ligation events.
# --min-mapq is the Mapping quality threshold for defining an alignment as a multi-mapping alignment. 
# --walks-policy is to handle multi mapping alignements.
# 5unique is used to report the 5’-most unique alignment on each side, if present (one or both sides may map to different locations on the genome, producing more than two alignments per DNA molecule)
if [ ! -e ${sample}.parsed.pairsam ]; then
  pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
    --nproc-in ${nbOfThreads} --nproc-out ${nbOfThreads} \
    --chroms-path "$filePathForSizesForBin" "${sample}.sam" > "${sample}.parsed.pairsam"
  
  if [ $? -ne 0 ]; then
    echo "Error: pairtools parse command failed for sample ${sample}"
    exit 1
  else
    echo "Sam size: $(ls ${sample}.sam)"
    echo "Parsed pairsam created, removing ${sample}.sam"
    rm -f "${sample}.sam"
  fi
else 
  echo "${sample}.parsed.pairsam already exists"    
fi

if [ ! -e ${sample}.sorted.pairsam ]; then
  # i specify the tmp dir cause I thinkin the cluster the automatic tmp location might have limits for the tmp file sizes
  tmpdir=$(mktemp -d /shared/projects/microc_pilot/tmp_${sample}_XXXXXX)
  # monitor size of the tmp dir. checks every 300 seconds (5 minutes)
  ( while true; do du -sh "$tmpdir"; sleep 20; done ) &
  monitor_pid=$!

  pairtools sort \
    --nproc ${nbOfThreads} \
    --tmpdir="$tmpdir" \
    "${sample}.parsed.pairsam" > "${sample}.sorted.pairsam"

  status=$?

  kill $monitor_pid

  if [ $status -ne 0 ]; then
    echo "❌ Error: pairtools sort command failed for sample ${sample}"
    rm -rf "$tmpdir"
    exit 1
  else
    echo "✅ Sorted pairsam created, removing parsed"
    rm -f "${sample}.parsed.pairsam"
    rm -rf "$tmpdir"
  fi
else
  echo "✅ ${sample}.sorted.pairsam already exists"
fi


# Remove PCR duplicates
# duplicate pairs are marked as DD in “pair_type” and as a duplicate in the sam entries.
if [ ! -e ${sample}.dedup.pairsam ]; then
  pairtools dedup --nproc-in ${nbOfThreads} --nproc-out ${nbOfThreads} --mark-dups --output-stats "stats.txt" --output "${sample}.dedup.pairsam" "${sample}.sorted.pairsam"
  if [ $? -ne 0 ]; then
    echo "Error: pairtools dedup command failed for sample ${sample}"
    exit 1
  else
    echo "Sorted pairsam size: $(ls ${sample}.sorted.pairsam)"
    echo "Dedup pairsam created, removing ${sample}.sorted.pairsam"
    rm -f "${sample}.sorted.pairsam"
  fi
else
  echo "${sample}.dedup.pairsam already exists"
fi

# Generate .pair and final BAM file
if [ ! -e ${sample}.mapped.pairs ]; then
  pairtools split --nproc-in ${nbOfThreads} --nproc-out ${nbOfThreads} --output-pairs "${sample}.mapped.pairs" --output-sam "${sample}.unsorted.bam" "${sample}.dedup.pairsam"

  if [ $? -ne 0 ]; then
    echo "Error: pairtools split command failed for sample ${sample}"
    exit 1
  else
    echo "dedup pairsam size: $(ls ${sample}.dedup.pairsam)"
    ## Sort and index the final BAM file
    # Can be used to generate a coverage adn library complexity
    tmpdir=$(mktemp -d /shared/projects/microc_pilot/tmp_${sample}_XXXXXX)
    echo "Sorting and then indexing ${sample}.unsorted.bam to calculate library complexity"
    samtools sort -@${nbOfThreads} -T ${tmpdir}/temp.bam -o "${sample}.mapped.PT.bam" "${sample}.unsorted.bam" 
    samtools index "${sample}.mapped.PT.bam"
    echo "Calculating library complexity"
    preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output ${sample}.preseq ${sample}.mapped.PT.bam
    if [ $? -ne 0 ]; then
      echo "Error: preseq lc_extrap command failed for sample ${sample}"
      exit 1
    else
      echo "library complexity sucessuly calculated with preseq. Deleting ${sample}.dedup.pairsam "
      rm -f "${sample}.dedup.pairsam"
    fi
  fi
else
  echo "Skipping ${sample}.mapped.pairs, it already exists"
fi


# Run QC script
python $dirPathForScripts/get_qc.py -p "${sample}.stats.txt"  > ${sample}.pretty.stats.txt

# # Save the stats in a common file for all samples
# mkdir $microc/output/stats_all/
# touch $microc/output/stats_all/stats_all.txt
# echo -e "$sample" >> $microc/output/stats_all/stats_all.txt
# cat "${sample}.stats.txt" >> "$microc/output/stats_all/"

# Generate contact matrix for .pairs with cooler
# Bgzip the pairs:
if [ ! -e ${sample}.pairs.gz ]; then
  bgzip -c "${sample}.mapped.pairs" > "${sample}.pairs.gz"
else
  echo "${sample}.pairs.gz already exists"
fi
# Index them
if [ ! -e ${sample}.pairs.gz.px2 ]; then
  pairix "${sample}.pairs.gz"
else
  echo "${sample}.pairs.gz already indexed"
fi
# Create the cooler at smaller resolution
if [ ! -e ${sample}_raw.${binSizeCoolMatrix}kb.cool ]; then
  cooler cload pairix -p 16 ${filePathForSizesForBin}:${binSizeCoolMatrix}000 "${sample}.pairs.gz" ${sample}_raw.${binSizeCoolMatrix}kb.cool
else
  echo "${sample}_raw.${binSizeCoolMatrix}kb.cool already exists"
fi
# Zoomify
if [ ! -e ${sample}.mcool ]; then
   singularity exec $pathToImages/cooler.0.10.3 cooler zoomify --balance -p ${nbOfThreads} --resolutions ${binSizeCoolMatrix}000N -o ${sample}.mcool --balance-args '--nproc 24 --cis-only' ${sample}_raw.${binSizeCoolMatrix}kb.cool
else
  echo "${sample}.mcool already exists"
fi

# Plot generation
ini_file=${sample}.ini
echo "[x-axis]" > $ini_file
for bin in $bins; do
  echo "[${sample}_${bin}kb]
file = ${sample}.mcool::/resolutions/${bin}000
depth = 2000000
min_value = 0
title = ${sample}_${bin}kb
file_type = hic_matrix
[spacer]
" >> ${ini_file}
done
echo "[x-axis]" >> ${ini_file}
# Generate a basic plot on the testRegion
pgt --tracks ${ini_file} --region ${testRegion} --fontSize 6 -o ${ini_file/.ini/_testRegion.pdf}

# Copy all final files to specific directories
mkdir -p ${dirPathWithResults}/allFinalFiles/cool
cp *.mcool ${dirPathWithResults}/allFinalFiles/cool/
mkdir -p ${dirPathWithResults}/allFinalFiles/reports
cp *.stats.txt ${dirPathWithResults}/allFinalFiles/reports/
mkdir -p ${dirPathWithResults}/allFinalFiles/pairs
cp ${sample}.pairs.gz ${dirPathWithResults}/allFinalFiles/pairs/
mkdir -p ${dirPathWithResults}/allFinalFiles/visualisation
cp *.pdf ${dirPathWithResults}/allFinalFiles/visualisation