# Use Docker for Micro-C analysis for MacOS

1. Download and install Docker Desktop from Docker's official website. Ensure Docker is running before proceeding.
Test if everything is ok
```
docker info
```
2. Write the Dockerfile with all the needed infos as in the example here

3. Build the image (add --no-cache if you want to build it from 0, otherwise it will take the last image built and update it if there are new lines in the dockerfile)
```
docker build --platform linux/amd64 -t microc .
```
--platform linux/amd64 tells to build the container with the linux architecture and not with the arm64 of the macOS

-t gives the name of the image

. tells to build the image in the directory where you are


4. Run the container
```
docker run --platform linux/amd64 -it microc
``` 
If you want to be able to access files in the host from the image, you need to mount volumes use the option 
Note: the new dir in the image cannot be named as a prexisting dir (built with the dockerfile), it would overwrite it.
```
docker run -it --rm -v /Users/cecilia.carmignoto/Documents/GitHub/Micro-C_pipeline/data:/data\
-v /Users/cecilia.carmignoto/Documents/GitHub/Micro-C_pipeline/output:/output \
--platform linux/amd64 microc . 
```
Now you are inside the container and can work in it interactively with the command line.

# Run the analysis step by step in the container

## Before starting
If test data are needed
```
wget https://s3.amazonaws.com/dovetail.pub/HiC/fastqs/MicroC_2M_R1.fastq \
wget https://s3.amazonaws.com/dovetail.pub/HiC/fastqs/MicroC_2M_R2.fastq
```
Reference genome, index file, chromosome file are needed. Download mm39 reference from UCSC:
```
wget -O mm39.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
```
Note: To use samtools faidx the genome has to be bgzipped, if differently zipped, unzip it and re-zip it with bgxzip
```
gunzip mm39.fa.gz \
bgzip mm39.fa
```
Index the genome:
```
samtools faidx mm39.fa.gz \
cut -f1,2 mm39.fa.gz.fai > mm39.genome \
bwa index mm39.fa.gz
```

## Step 1 Generate sam file
Note: the reference genome has to be indexed. 
```
bwa mem -5SP -T0 -t8 references/mm39/mm39.fa.gz MicroC_2M_R1.fastq MicroC_2M_R2.fastq > MicroC_2M_aligned.sam
```

## Step 2. Record valid ligation events 
```
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 5 --nproc-out 5 --chroms-path references/mm39.genome MicroC_2M_aligned.sam > parsed.pairsam
```

## Step 3. Sort the parsed.pairsam
```
pairtools sort --nproc 5 --tmpdir=/tmp parsed.pairsamgenompy > sorted.pairsam
```

## Step 4. Remove PCR duplicates
```
pairtools dedup --nproc-in 5 --nproc-out 5 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
```

# Run the analysis directly when running the container

## Analysis.sh
Set the varaibles at the beginning of the script Analysis.sh. 
```
GENOME="hg38"                      # Genome version
GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.fa.gz"
CORES=5                            # Number of processing cores
DATA="data"                        # Directory containing the input FASTQ files
```

Check that the results are saved in the output dir.
Check in the data folder fastq are called R1.fastq and R2.fastq (line 30 Analysis.sh)

## Modify the Dockerfile

Copy the the script with the commands for the analysis (Analysis.sh). 
Check in the Analysis.sh that the results are saved in the output dir.
The file you want to copy in the image, has to be in the directory you are building the container from and it has to be executable in the host. Make it executable
```
chmod +x Analysis.sh
```
Add line to copy in the Dockerfile
```
COPY Analysis.sh /Analysis.sh
```
Set the new entrypoint. ("sh" is make it executable, no need it you did it already in the host, but just to make sure)
```
ENTRYPOINT ["sh","/Analysis.sh"]
```
And comment the line ENTRYPOINT ["/bin/bash"]

## Build, run image and mount the volumes
Build image
```
docker build --platform linux/amd64 -t microc .
```

Explained above. the dir data will have the fastq and the results will be saved in the dir output. 
```
docker run --rm -v /Users/cecilia.carmignoto/Documents/GitHub/Micro-C_pipeline/data:/data \
-v /Users/cecilia.carmignoto/Documents/GitHub/Micro-C_pipeline/output:/output \
--platform linux/amd64 microc . 
```


# For me
- download the data directly on the server? (is it gonna be still there) copy it from our server? --> we have to mount it 
- we have to iterate for all samples. Loop in the Analyis.sh? Parallelize how?
.