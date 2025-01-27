# Get our data from s3 with AWS

# Pull aws image

# Save the samples fastq in 
pathToFastq="$microcPilot/fastq/"
pathToImages="$SRC/images"

# Pull Images
# aws
wget -nc -O $pathToImages/awscli.1.8.3.sif "http://datacache.galaxyproject.org/singularity/all/awscli:1.8.3--py35_0"
function aws() {
  singularity exec $pathToImages/awscli.1.8.3.sif aws $*
}

# check it is properly installed
aws --version 

# Configure with infos form the email of BGI
aws configure
# Aws_access_key_id：AKIAYWZZRVKWRKSBNUF3
# Aws_secret_access_key：w7c1VsxjyVtg4ryS2FrZHqTFUFd/2E8/q5mbRH/q
# Region：eu-west-3
# Defaul output format:

# Test if it finds the data and summarize them
aws s3 ls --summarize --human-readable --recursive s3://musgzjor-598731762349/F24A430002451_MUSgzjoR_24DEC2024

# Download
aws s3 sync s3://musgzjor-598731762349/F24A430002451_MUSgzjoR_24DEC2024/ $pathToFastq > permanent_transfert.log 2> permanent_transfert.err  # better sync then cp, sync doesnt start over again if fails
