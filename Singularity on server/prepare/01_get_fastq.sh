## The table genomes_table.txt has to be already generated. (See README.md)
# first column is the genome name
# second column is the absolute path for fasta

pathToGenomesTable="$HOME/genomes_table.txt"
genome="1"       # Set the line number of genomes_table.txt of the genome to download (line 1 hg38, line 2 mm39)

mkdir -p $HOME/genomes
cd $HOME/genomes

# Get the genome name and fasta file from the table
genome=$(cat ${pathToGenomesTable} | awk -v i=$genome 'NR==i{print $1}')
filePathForFasta=$(cat ${pathToGenomesTable} | awk -v i=$genome 'NR==i{print $2}')

# Download 
wget -O $genome filePathForFasta
gunzip -k $genome > genome_unzipped.fa
bgzip genome_unzipped.fa > $genome # is it okay to run this in the front end?

# Get our data from s3 with AWS
mkdir="$HOME/fastq"
cd $HOME/fastq





URL：https://s3.console.aws.amazon.com/s3/buckets/musgzjor-598731762349?region=eu-west-3&tab=objects

Project：F24A430002451_MUSgzjoR

Alias ID：598731762349

S3 Bucket：musgzjor-598731762349

Account：musgzjor

Password：nL8]Y|$K|1Q5

Region：eu-west-3

Aws_access_key_id：AKIAYWZZRVKWRKSBNUF3

Aws_secret_access_key：w7c1VsxjyVtg4ryS2FrZHqTFUFd/2E8/q5mbRH/q