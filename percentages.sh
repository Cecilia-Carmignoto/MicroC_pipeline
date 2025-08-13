#!/bin/bash

#SBATCH -o slurm-%x-%A_%a.out  # Correct formatting
#SBATCH -e slurm-%x-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --array=1-26
#SBATCH --job-name=runMicroC

# Directories
dirPathForFastq="${microcPilot}/fastq/"
filePathForTable="${dirPathForFastq}/samplesFastqTable.txt"
dirPathWithResults="$microcPilot/outputs/"
dirPathDiagnosis="$microcPilot/outputs/diagnose/"

# Extract sample name
sample=$(awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}' "$filePathForTable")

# Calculate percentages and save to files
for region in full internal start end; do
    awk 'FNR==1 { if (NR==1) num1=$1; else print (num1 / $1) * 100 }' \
    "$dirPathDiagnosis/${sample}_${region}_R1.txt" "$dirPathDiagnosis/${sample}_total_reads_R1.txt" \
    > "$dirPathDiagnosis/${sample}_percentage_${region}_R1.txt"
done

# Create summary file
summary_file="$dirPathDiagnosis/summary_percentages.txt"
echo -e "Sample\tFull\tInternal\tStart\tEnd" > "$summary_file"

awk 'NR>1 {print $1}' "$filePathForTable" | while read sample; do
    full=$(cat "$dirPathDiagnosis/${sample}_percentage_full_R1.txt")
    internal=$(cat "$dirPathDiagnosis/${sample}_percentage_internal_R1.txt")
    start=$(cat "$dirPathDiagnosis/${sample}_percentage_start_R1.txt")
    end=$(cat "$dirPathDiagnosis/${sample}_percentage_end_R1.txt")

    echo -e "$sample\t$full\t$internal\t$start\t$end" >> "$summary_file"
done
