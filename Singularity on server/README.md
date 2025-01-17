# Check $HOME on server where I am

#To do first time on the server
```
echo 'export SRC=$HOME/microc' >> $HOME/.bashrc 
source $HOME/.bashrc 
mkdir -p $SRC
mkdir $SRC/genomes
```

## Create table for the genomes hg38 and mm39
```
echo -e "hg38.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz\nmm39.fa.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz" > $SRC/genomes/genomes_table.txt
```



