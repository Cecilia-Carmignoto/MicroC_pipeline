
# Introductiom
Here you find the scripts to run a Micro-C analysis. 

To handle dependencies singularity images have been used from http://datacache.galaxyproject.org/singularity/all/.

The pipeline follows https://micro-c.readthedocs.io/en/latest/index.html.

## Server set up
On the mesoPSL server:
Check what /obs and /travail are. 
./obs is to save only scripts, /travail is for everything else.
In mesoPSL srver $HOME is /obs/ijerkovic. Check how much space I have.

Set the dir $SRC as the source dir where everyhing starts from. Then the $microc dir inside where to do the microC analysis.
$SRC/cecilia is for SLURM error messages and outputs if no absolut path is defined inthe scripts.

For github set up: generate tocken form Setting > Developer settings > Personal acess tokens > Tokens (classic) > Generate new tocken. A Personal Access Tocken wil be genrated , that is the password to use when doing git clone.
```
cd $HOME 
git clone https://github.com/Cecilia-Carmignoto/Micro-C_pipeline
checkout Lucille
```

Export variables for dirs in the ./bashrc
```
echo 'export SRC=/scratch/ldelisle/NGS' >> $HOME/.bashrc
echo 'export microc=$SRC/microc' >> $HOME/.bashrc 
echo 'export microcPilot=$SRC/microc/pilot' >> $HOME/.bashrc 
echo 'export microcFullData=$SRC/microc/fullData' >> $HOME/.bashrc 
echo 'export PREP=$HOME/Micro-C_pipeline/Singularity/prepare' >> $HOME/.bashrc
echo 'export RUN=$HOME/Micro-C_pipeline/Singularity/run' >> $HOME/.bashrc
source $HOME/.bashrc 
```

Create all the needed directories.
```
mkdir -p $SRC
mkdir -p $microc
mkdir -p $microcPilot
mkdir -p $microcFullData
mkdir -p $SRC/genomes
mkdir -p $SRC/images
```

Make scripts executable
```
chmod +x $PREP/01.1_get_genome.sh
chmod +x $PREP/01.2_get_fastq.sh
chmod +x $PREP/03_fastq_table.sh
chmod +x $PREP/02_bwa_index.sh
chmod +x $RUN/04_from_fastq_to_valid_pairs_and_mcool.sh
```

## Create reference genome table
Generate table for the reference genome data.
In the genomes table: first column is the genome name in the format hg38.fa.gz, second column is the http of the fastq to be downloaded.
```
echo -e "hg38\t$SRC/genomes/fasta/hg38.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\nmm39\t$SRC/genomes/fasta/mm39.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" > $SRC/genomes/genomesTable.txt
```

## Download data and genome

For the genome. open script and check paths. set the variable genomeLine to 1 if you want hg38, 2 if you want mm39.
```
bash $PREP/01.1_get_genome.sh
```

Get the fastq files form s3 with aws.
```
bash $PREP/01.2_get_fastq.sh
```

## Index the genome
-chdir $SRC is to set the dir where to write outputs that have no path indication in the scripts and for the .log and error files.
Modify the SBATCH --array=x with the row/rows from the table that need to be processed in the table
```
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
```
bash --chdir $SRC $RUN/04_from_fastq_to_valid_pairs_and_mcool.sh
```