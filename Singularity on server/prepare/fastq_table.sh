#!/bin/bash

# Define the parent directory containing sample folders
parent_dir="$HOME/fastq/F24A430002451_MUSgzjoR_24DEC2024"

# Initialize the output file
output_file="$HOME/fastq/samples_fastq_table.txt"

# Iterate over each sample directory
for sample_dir in "$parent_dir"/*/; do
    # Extract the sample name (folder name)
    sample_name=$(basename "$sample_dir")

    # Find the FASTQ files ending with 1.fq.gz and 2.fq.gz
    read1_file=$(find "$sample_dir" -type f -name '*1.fq.gz')
    read2_file=$(find "$sample_dir" -type f -name '*2.fq.gz')

    # Append the sample information to the output file
    echo -e "$sample_name\t$read1_file\t$read2_file" >> "$output_file"
done