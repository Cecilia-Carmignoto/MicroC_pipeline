# Get our data from s3 with AWS

# Pull aws image

# Save the samples fastq in 
pathToFastq="$microcPilot/fastq/"
pathToImages="$SRC/images"

# Pull Images
# aws
if [ ! -e $pathToImages/aws-cli_2.23.6.sif ]; then
    singularity pull docker://amazon/aws-cli:2.23.6
    mv aws-cli_2.23.6.sif $pathToImages/
fi
function aws() {
  singularity exec $pathToImages/aws-cli_2.23.6.sif aws $*
}
# Give access to pathToFastq
export APPTAINER_BIND=$SRC

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
