#!/bin/bash

# Script to generate the samplesFastqTable where:
# firt column is the sample name
# second column is the fastq1 path
# third column is the fasq2 path

# Define the parent directory containing sample folders
mkdir $SRC/fastq/
pathToFastq="$SRC/fastq/"

# Initialize the output file
output_file="$pathToFastq/samplesFastqTable.txt"

# Iterate over each sample directory
for sample_dir in "$pathToFastq"/*/; do
    # Extract the sample name (folder name)
    sample_name=$(basename "$sample_dir")

    # Find the FASTQ files ending with 1.fq.gz and 2.fq.gz
    read1_file=$(find "$sample_dir" -type f -name '*1.fq.gz')
    read2_file=$(find "$sample_dir" -type f -name '*2.fq.gz')

    # Append the sample information to the output file
    echo -e "$sample_name\t$read1_file\t$read2_file" >> "$output_file"
done