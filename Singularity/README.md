
# Introductiom
Here you find the scripts to run a Micro-C analysis. 

To handle dependencies singularity images have been used from https://depot.galaxyproject.org/singularity/.

The pipeline follows https://micro-c.readthedocs.io/en/latest/index.html.

## Server set up
On the mesoPSL server:
Check what /obs and /travail are. 
./obs is to save only scripts, /travail is for everything else.
In mesoPSL srver $HOME is /obs/ijerkovic. Check how much space I have.

Set the dir $SRC as the source dir where everyhing starts from. Then the $microc dir inside where to do the microC analysis.
$SRC/cecilia is for SLURM error messages and outputs if no absolut path is defined inthe scripts.

For github set up: generate a ssh key on the machine with ssh-keygen and put the public key on github.

```bash
cd $HOME 
git clone git@github.com:Cecilia-Carmignoto/Micro-C_pipeline.git
cd Micro-C_pipeline
git checkout Lucille
```

Export variables for dirs in the ./bashrc
```bash
echo '## FOR MICROC' >> $HOME/.bashrc
echo 'export SRC=/scratch/ldelisle/NGS' >> $HOME/.bashrc
echo 'export microc=$SRC/microc' >> $HOME/.bashrc 
echo 'export microcPilot=$SRC/microc/pilot' >> $HOME/.bashrc 
echo 'export microcFullData=$SRC/microc/fullData' >> $HOME/.bashrc 
echo 'export PREP=$HOME/Micro-C_pipeline/Singularity/prepare' >> $HOME/.bashrc
echo 'export RUN=$HOME/Micro-C_pipeline/Singularity/run' >> $HOME/.bashrc
source $HOME/.bashrc 
```

Create all the needed directories.
```bash
mkdir -p $SRC
mkdir -p $microc
mkdir -p $microcPilot
mkdir -p $microcFullData
mkdir -p $SRC/genomes
mkdir -p $SRC/images
```

Make scripts executable
```bash
chmod +x $PREP/01.1_get_genome.sh
chmod +x $PREP/01.2_get_fastq.sh
chmod +x $PREP/03_fastq_table.sh
chmod +x $PREP/02_bwa_index.sh
chmod +x $RUN/04_from_fastq_to_valid_pairs_and_mcool.sh
```

## Create reference genome table
Generate table for the reference genome data.
In the genomes table: first column is the genome name in the format hg38.fa.gz, second column is the http of the fastq to be downloaded.
```bash
echo -e "hg38\t$SRC/genomes/fasta/hg38.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\nmm39\t$SRC/genomes/fasta/mm39.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" > $SRC/genomes/genomesTable.txt
```

## Download data and genome

For the genome. open script and check paths. set the variable genomeLine to 1 if you want hg38, 2 if you want mm39.
```bash
bash $PREP/01.1_get_genome.sh
```

Get the fastq files form s3 with aws.
```bash
bash $PREP/01.2_get_fastq.sh
```

## Index the genome
-chdir $SRC is to set the dir where to write outputs that have no path indication in the scripts and for the .log and error files.
Modify the SBATCH --array=x with the row/rows from the table that need to be processed in the table
```bash
sbatch --chdir $SRC $PREP/02_bwa_index.sh
```

## Create Sequencing data reference table
Generate tables for the sequencing data.
In the samplesFastqTable: first column is the sample name, second column is the fastq1 path, third column is the fastq2 path.
CHECK: fastq names have to end in '1.fq.gz' (for read 1), '2.fq.gz' (for read 2)

```
bash $PREP/03_fastq_table.sh
```

## MicroC analysis

The script needs the indexed reference genome
```bash
sbatch --chdir $SRC $RUN/04_from_fastq_to_valid_pairs_and_mcool.sh
```

## Get public data to compare the quality

First get the fastq:

```bash
pathToImages="$SRC/images"
wget -nc -O "$pathToImages/sra-tools:3.1.1--h4304569_2"  "https://depot.galaxyproject.org/singularity/sra-tools:3.1.1--h4304569_2"
function fasterq-dump() {
  singularity exec "$pathToImages/sra-tools:3.1.1--h4304569_2" fasterq-dump $*
}
export APPTAINER_BIND=$SRC
dirPathForFastq="${microcPilot}/fastq/"
cd $dirPathForFastq
sample=SRR29294642
fasterq-dump -o ${sample}.fastq ${sample}
```

Get only 50 million reads:

```bash
for r in 1 2; do
    head -n 200000000 ${sample}_${r}.fastq | gzip > ${sample}_50M_${r}.fastq.gz
done
```

Trim to 50bp

```bash
toolversion="seqtk:1.4--h577a1d6_3"
wget -nc -O "$pathToImages/${toolversion}"  "https://depot.galaxyproject.org/singularity/${toolversion}"
function seqtk() {
  singularity exec "$pathToImages/seqtk:1.4--h577a1d6_3" seqtk $*
}
for r in 1 2; do
    seqtk trimfq -l 50 ${sample}_50M_${r}.fastq.gz > ${sample}_50M_50bp_${r}.fastq.gz
done
```

Add these 2 to the table with fastqs:

```bash
echo -e "${sample}_50M_full\t${sample}_50M_1.fastq.gz\t${sample}_50M_2.fastq.gz
${sample}_50M_50bp\t${sample}_50M_50bp_1.fastq.gz\t${sample}_50M_50bp_2.fastq.gz" >> samplesFastqTable.txt
```

## Aggregate all reports in one

```bash
pathToImages="$SRC/images"
toolversion="multiqc:1.26--pyhdfd78af_0"
wget -nc -O "$pathToImages/${toolversion}"  "https://depot.galaxyproject.org/singularity/${toolversion}"
function multiqc() {
  singularity exec "$pathToImages/multiqc:1.26--pyhdfd78af_0" multiqc $*
}
export APPTAINER_BIND=$SRC
cd $microcPilot/outputs
multiqc . -m pairtools --force
```

Joins the pretty reports:

```bash
all_pretty=$(ls allFinalFiles/reports/*pretty*.tsv)
file1=""
for file in $all_pretty; do
  if [ -z $file1 ]; then
    sample=$(basename ${file/.pretty.stats.tsv//})
    echo -e "Name\tnb_$sample\tpct_$sample" | cat - $file > join.tmp
    file1="OK"
  else
    sample=$(basename ${file/.pretty.stats.tsv//})
    echo -e "Name\tnb_$sample\tpct_$sample" | cat - $file > join.input2.tmp
    join -t$'\t' --nocheck-order join.tmp join.input2.tmp > join.tmp.1
    mv join.tmp.1 join.tmp
  fi
done
rm join.input2.tmp
mv join.tmp combined.stats.pretty.tsv
```

## Make a general plot

```bash
pathToImages="$SRC/images"
toolversion="pygenometracks:3.9--pyhdfd78af_0"
wget -nc -O "$pathToImages/${toolversion}"  "https://depot.galaxyproject.org/singularity/${toolversion}"
function pgt() {
  singularity exec "$pathToImages/pygenometracks:3.9--pyhdfd78af_0" pgt $*
}
export APPTAINER_BIND=$SRC
cd $microcPilot/outputs
bin=50
testRegion="chr2:73150000-76150000"
ini_file=all_${bin}kb.ini
echo "[x-axis]" > $ini_file
for mcool in allFinalFiles/cool/*; do
    echo "[$mcool]
file = ${mcool}::/resolutions/${bin}000
depth = 2000000
min_value = 0
title = $(basename ${mcool/.mcool/})_${bin}kb
file_type = hic_matrix
[spacer]
" >> ${ini_file}
done
echo "[x-axis]" >> ${ini_file}
# Generate a basic plot on the testRegion
pgt --tracks ${ini_file} --region ${testRegion} --fontSize 6 -o ${ini_file/.ini/_testRegion.pdf}
```