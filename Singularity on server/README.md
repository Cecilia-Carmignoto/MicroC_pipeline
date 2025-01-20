

# Check $HOME on server where I am

#To do first time on the server
```
echo 'export SRC=$HOME/microc' >> $HOME/.bashrc 
source $HOME/.bashrc 
mkdir -p $SRC
mkdir $SRC/genomes
```

## Reference genome table
Generate tables for the reference genome and the sequencing data.
In the genomes table: first column is the genome name in the format hg38.fa.gz, second column is the http of the fastq to be downloaded.
```
mkdir -p $SRC/genomes
echo -e "hg38.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\nmm39.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" > $SRC/genomes/genomesTable.txt
```


## Download data and genome

```
./01.1_get_genome.sh
```


## Sequencing data reference table
In the samplesFastqTable where, first column is the sample name, second column is the fastq1 path, third column is the fastq2 path.
CHECK: fastq names have to end in '1.fq.gz' (for read 1), '2.fq.gz' (for read 2)
```
./03_fastq_table.sh
```


