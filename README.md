# Use Docker for Micro-C analysis for MacOS

1. Download and install Docker Desktop from Docker's official website. Ensure Docker is running before proceeding.
Test if everything is ok
```
docker info
```
2. Write the Dockerfile with all the needed infos as in the example here

3. Build the image
```
docker build --platform linux/amd64 -t microc .
```
--platform linux/amd64 tells to build the container with the linux architecture and not with the arm64 of the macOS

-t gives the name of the image

. tells to build the image in the directory where you are

If you want to be able to access som efiles in the external from the image, add --rm -v <dir_you_want_to_acess>:/<where_you_want_it_on_the_image>
Note: it cannot be a prexisting dir (built with the dockerfile), it would overwrite it, deleting everyhting there is inside

4. Run the container
```
docker run --platform linux/amd64 -it microc
``` 
Now you are inside the container and can work in it interactively with the command line

In this specifc Dockerfile there is Micro-C tool (https://micro-c.readthedocs.io/en/latest/before_you_begin.html) with all its dependencies, in the correct versions. 

# Run the analysis step by step

## Before starting

If test data are needed
```
wget https://s3.amazonaws.com/dovetail.pub/HiC/fastqs/MicroC_2M_R1.fastq\nwget https://s3.amazonaws.com/dovetail.pub/HiC/fastqs/MicroC_2M_R2.fastq
```
Reference genome, index file, chromosome file are needed. I downloaded mm39 reference. 

Note: To use samtools faidx the genome has to be bgzipped, if differntly zipped, unzip it and re-zip it with bgxzip
```
samtools faidx mm39.fa.gz \
cut -f1,2 mm39.fa.gz.fai > mm39.genome \
bwa index mm39.fa.gz
```
If the reference genome and the other files are already present somewhere external to the  image you can, while building it do 
```
docker run -it --rm -v /Users/cecilia.carmignoto/references:/references --platform linux/amd64 microc 
```

## Step 1 Generate sam file

Note: the reference genome has to be indexed. 

```
bwa mem -5SP -T0 -t8  ../../references/mm39/mm39.fa.gz ../Micro-C/test_data/MicroC_2M_R1.fastq ../Micro-C/test_data/MicroC_2M_R2.fastq > MicroC_2M_aligned.sam
```

# Step 2. Record valid ligation events 

```
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 5 --nproc-out 5 --chroms-path ../../references/mm39/mm39.genome MicroC_2M_aligned.sam > parsed.pairsam
```
docker run -it --rm -v /path/to/local/dir:/path/in/container image_name




