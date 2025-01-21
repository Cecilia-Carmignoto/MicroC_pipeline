
# Introductiom
Here you find the scripts to run a Micro-C analysis. 

To handle dependencies singularity images have been used from http://datacache.galaxyproject.org/singularity/all/.

The pipeline follows https://micro-c.readthedocs.io/en/latest/index.html.

## Server set up
On the mesoPSL server:
Check what /obs and /travail are. 
./obs is to save only scripts, /travail is for everything else.
Check $HOME. Check how much space I have.

Set the dir $SRC as the source dir where everyhing starts from. Then the $microc dir inside where to do the microC analysis.
$SRC/cecilia is for SLURM error messages.
```
echo 'export SRC=$HOME/NGS' >> $HOME/.bashrc
echo 'export microc=$HOME/NGS/microc' >> $HOME/.bashrc 
source $HOME/.bashrc 
mkdir -p $SRC
mkdir -p $microc
mkdir $SRC/genomes
$SRC/cecilia 
```

## Create reference genome table
Generate table for the reference genome data.
In the genomes table: first column is the genome name in the format hg38.fa.gz, second column is the http of the fastq to be downloaded.
```
mkdir -p $SRC/genomes
echo -e "hg38.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\nmm39.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" > $SRC/genomes/genomesTable.txt
```

## Download data and genome

For the genome. open script and check paths. set the variable genomeLine to 1 if you want hg38, 2 if you want mm39.
```
./01.1_get_genome.sh
```

get the fastq files form s3 with aws.
```
./01.2_get_fastq.sh
```

## Create Sequencing data reference table
Generate tables for the sequencing data.
In the samplesFastqTable: first column is the sample name, second column is the fastq1 path, third column is the fastq2 path.
CHECK: fastq names have to end in '1.fq.gz' (for read 1), '2.fq.gz' (for read 2)
```
./03_fastq_table.sh
```

## Index the genome


## MicroC analysis

The script needs the indexed reference genome
```
./04_from_fastq_to_valid_pairs_and_mcool.sh
```

