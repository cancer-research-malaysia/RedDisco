#!/usr/bin/env bash

# Check bash version - associative arrays require bash 4.0+
if [ "${BASH_VERSION%%.*}" -lt 4 ]; then
    echo "Error: This script requires bash 4.0 or higher (current: $BASH_VERSION)"
    echo "Please run with: bash $0 $@"
    exit 1
fi

# Script to generate a manifest file from S3 path
# Usage: ./generate_manifest.sh <s3_path> [output_file]

set -e  # Exit on any error

# Function to display usage
usage() {
    echo "Usage: $0 <s3_path> [output_file]"
    echo ""
    echo "Arguments:"
    echo "  s3_path      - S3 bucket path to where the fastq input files are stored (e.g., s3://bucket/path/to/fastq/files/)"
    echo "  output_file  - Output manifest file (default: manifest.tsv)"
    echo ""
    echo "Sample Type Inference:"
    echo "  - Sample names ending with 'T' -> Tumor"
    echo "  - Sample names ending with 'N' -> Normal"
    echo "  - Samples not following this convention will be skipped"
    echo ""
    echo "Examples:"
    echo "  $0 s3://my-bucket/data/RNAseq/ my_manifest.tsv"
    exit 1
}

# Check if AWS CLI is installed
if ! command -v aws &> /dev/null; then
    echo "Error: AWS CLI is not installed or not in PATH"
    echo "Please install AWS CLI: https://aws.amazon.com/cli/"
    exit 1
fi

# Parse arguments
if [ $# -lt 1 ]; then
    usage
fi

S3_PATH="$1"
OUTPUT_FILE="${2:-manifest.tsv}"

# Remove trailing slash from S3 path if present
S3_PATH="${S3_PATH%/}"

# Extract bucket name from S3 path
BUCKET_NAME=$(echo "$S3_PATH" | sed 's|s3://||' | cut -d'/' -f1)

echo "Generating manifest from S3 path: $S3_PATH"
echo "Bucket name: $BUCKET_NAME"
echo "Output file: $OUTPUT_FILE"
echo "Sample types will be inferred from sample name suffix (T=Tumor, N=Normal)"

# Create temporary file for processing
TEMP_FILE=$(mktemp)
trap "rm -f $TEMP_FILE" EXIT

# List all fastq files in the S3 path
echo "Listing files from S3..."
aws s3 ls "$S3_PATH/" --recursive | grep -E "\.(fq|fastq)(\.gz)?$" | awk '{print $4}' > "$TEMP_FILE"

if [ ! -s "$TEMP_FILE" ]; then
    echo "Error: No FASTQ files found in $S3_PATH"
    echo "Make sure the S3 path is correct and contains .fq, .fastq, .fq.gz, or .fastq.gz files"
    exit 1
fi

echo "Found $(wc -l < "$TEMP_FILE") FASTQ files"

# Create the manifest file with header
echo -e "sampleName\tread1Path\tread2Path\tsampleType" > "$OUTPUT_FILE"

# Process files to create sample pairs
# Create temporary files for R1 and R2 lists
R1_FILE=$(mktemp)
R2_FILE=$(mktemp)
trap "rm -f $TEMP_FILE $R1_FILE $R2_FILE" EXIT

while IFS= read -r file_path; do
    # Extract filename from path
    filename=$(basename "$file_path")
    
    # Extract sample name by removing read indicators and file extensions
    # This handles patterns like: sampleName_r1.fq.gz, sampleName_R1.fastq, sampleName_1.fq, etc.
    sample_name=$(echo "$filename" | sed -E 's/[._-][rR]?[12][._-].*$//' | sed -E 's/[._-][12]\..*$//')
    
    # Determine if this is read 1 or read 2
    if echo "$filename" | grep -qE "[._-][rR]?1[._-]"; then
        echo "${sample_name}|s3://$BUCKET_NAME/$file_path" >> "$R1_FILE"
    elif echo "$filename" | grep -qE "[._-][rR]?2[._-]"; then
        echo "${sample_name}|s3://$BUCKET_NAME/$file_path" >> "$R2_FILE"
    else
        echo "Warning: Could not determine read number for file: $filename"
    fi
done < "$TEMP_FILE"

# Generate manifest entries
processed_samples=()
skipped_samples=()
incomplete_pairs=()

# Function to find file path for a sample and read number
find_file() {
    local sample="$1"
    local read_num="$2"
    local file_list="$3"
    
    grep "^${sample}|" "$file_list" | cut -d'|' -f2
}

# Get unique sample names from R1 files
sample_names=$(cut -d'|' -f1 "$R1_FILE" | sort -u)

for sample_name in $sample_names; do
    # Find corresponding R1 and R2 files
    r1_path=$(find_file "$sample_name" "1" "$R1_FILE")
    r2_path=$(find_file "$sample_name" "2" "$R2_FILE")
    
    # Check if we have both R1 and R2 for this sample
    if [[ -z "$r1_path" ]] || [[ -z "$r2_path" ]]; then
        echo "Warning: Skipping sample '$sample_name' - incomplete R1/R2 pair"
        if [[ -z "$r2_path" ]]; then
            echo "  Missing R2 file for sample: $sample_name"
        elif [[ -z "$r1_path" ]]; then
            echo "  Missing R1 file for sample: $sample_name"
        fi
        incomplete_pairs+=("$sample_name")
        continue
    fi
    
    # Infer sample type from sample name suffix
    sample_type=""
    if [[ "$sample_name" =~ T$ ]]; then
        sample_type="Tumor"
    elif [[ "$sample_name" =~ N$ ]]; then
        sample_type="Normal"
    else
        echo "Warning: Skipping sample '$sample_name' - does not end with 'T' or 'N'"
        skipped_samples+=("$sample_name")
        continue
    fi
    
    # Add to manifest
    echo -e "$sample_name\t$r1_path\t$r2_path\t$sample_type" >> "$OUTPUT_FILE"
    processed_samples+=("$sample_name")
done

# Report results
sample_count=${#processed_samples[@]}
skipped_count=${#skipped_samples[@]}
incomplete_count=${#incomplete_pairs[@]}
echo ""
echo "Manifest generation complete!"
echo "Processed $sample_count sample pairs"
if [ $incomplete_count -gt 0 ]; then
    echo "Skipped $incomplete_count samples (incomplete R1/R2 pairs): ${incomplete_pairs[*]}"
fi
if [ $skipped_count -gt 0 ]; then
    echo "Skipped $skipped_count samples (not ending with T or N): ${skipped_samples[*]}"
fi
echo "Output written to: $OUTPUT_FILE"

if [ $sample_count -eq 0 ]; then
    echo ""
    echo "Warning: No valid sample pairs found."
    echo "Please check that your files follow a naming convention like:"
    echo "  - sampleName_r1.fq.gz and sampleName_r2.fq.gz"
    echo "  - sampleName_R1.fastq and sampleName_R2.fastq"
    echo "  - sampleName_1.fq and sampleName_2.fq"
fi

echo ""
echo "Preview of generated manifest:"
head -n 5 "$OUTPUT_FILE"
if [ $(wc -l < "$OUTPUT_FILE") -gt 5 ]; then
    echo "... (showing first 4 samples)"
fi
