# TO RUN ONLY ONCE
# Download of the reference genome
# Download of the samples data
# at line 11 choose the reference genome. (check genomes_table.txt and chose the line)

# Define paths
# CHECK: The table genomes_table.txt has to be already generated. (See README.md)
# first column is the genome name
# second column is the absolute path for fasta
pathToGenomesTable="$SCR/genomes/genomesTable.txt"
genomeLine=1            # Set the line number of genomesTable.txt of the genome to download (line 1 hg38, line 2 mm39)
mkdir -p $SCR/genomes
mkdir -p $SRC/images
pathToImages="$SRC/images"
cd $SCR/genomes

# Pull Images
# bgzip is in the image of samtools
# samtools
wget -nc -O $pathToImages/samtools.1.11.sif ''http://datacache.galaxyproject.org/singularity/all/samtools:1.11--h6270b1f_0''
function samtools() {
singularity exec $pathToImages/samtools.1.11.sif samtools $*
}
# bgzip from samtools
function bgzip() {
singularity exec $pathToImages/samtools.1.11.sif bgzip $*
}

# Get the genome name and fasta file from the genomes_table.txt
genomeName=$(cat ${pathToGenomesTable} | awk -v i=$genomeLine 'NR==i{print $1}') 
genomeURL=$(cat ${pathToGenomesTable} | awk -v i=$genomeLine 'NR==i{print $2}')

# Download genome
wget -nc -O $genomeName $genomeURL
gunzip -k $genomeName > genomeUnzipped.fa
bgzip genomeUnzipped.fa > $genomeName # is it okay to run this in the front end?
rm genomeUnzipped.fa



# Get our data from s3 with AWS

# Save the samples fastq in 
pathToFastq="$SRC/fastq/"

# URL：https://s3.console.aws.amazon.com/s3/buckets/musgzjor-598731762349?region=eu-west-3&tab=objects

# Project：F24A430002451_MUSgzjoR

# Alias ID：598731762349

# S3 Bucket：musgzjor-598731762349

# Account：musgzjor

# Password：nL8]Y|$K|1Q5

# Region：eu-west-3

# Aws_access_key_id：AKIAYWZZRVKWRKSBNUF3

# Aws_secret_access_key：w7c1VsxjyVtg4ryS2FrZHqTFUFd/2E8/q5mbRH/q