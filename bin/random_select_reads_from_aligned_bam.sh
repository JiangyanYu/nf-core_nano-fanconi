#!/bin/bash

# Usage: ./process_bam.sh input.bam output.bam

# Input and output files
INPUT_BAM=$1
OUTPUT_BAM=$2

# Temporary files
READ_NAMES_FILE="all_read_names.txt"
SAMPLED_READS_FILE="sampled_read_names.txt"

# Ensure required commands are available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed or not in PATH."
    exit 1
fi

if ! command -v shuf &> /dev/null; then
    echo "Error: shuf is not installed or not in PATH."
    exit 1
fi

# Step 1: Extract all read names from the BAM file
echo "Extracting all read names from the BAM file..."
samtools view "$INPUT_BAM" | awk '{print $1}' > "$READ_NAMES_FILE"

# Step 2: Randomly select 10% of the read names
echo "Sampling 10% of the read names..."
TOTAL_READS=$(wc -l < "$READ_NAMES_FILE")
SAMPLE_SIZE=$((TOTAL_READS / 10))
shuf -n "$SAMPLE_SIZE" "$READ_NAMES_FILE" > "$SAMPLED_READS_FILE"

# Step 3: Extract aligned records for the sampled read names
echo "Extracting aligned records for sampled reads..."
samtools view -h "$INPUT_BAM" | grep -F -f "$SAMPLED_READS_FILE" | samtools view -bo "$OUTPUT_BAM"

# Cleanup
#rm "$READ_NAMES_FILE" "$SAMPLED_READS_FILE"

echo "Process completed. Output BAM file: $OUTPUT_BAM"
