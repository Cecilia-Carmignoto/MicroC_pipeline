mkdir -p references
cd references

# Get mm39
wget -O mm39.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz 
bgzip mm39.fa #is it okay to run this in the front end?

# Get our data from s3 with AWS





URL：https://s3.console.aws.amazon.com/s3/buckets/musgzjor-598731762349?region=eu-west-3&tab=objects

Project：F24A430002451_MUSgzjoR

Alias ID：598731762349

S3 Bucket：musgzjor-598731762349

Account：musgzjor

Password：nL8]Y|$K|1Q5

Region：eu-west-3

Aws_access_key_id：AKIAYWZZRVKWRKSBNUF3

Aws_secret_access_key：w7c1VsxjyVtg4ryS2FrZHqTFUFd/2E8/q5mbRH/q