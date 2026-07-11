#!/bin/bash

# Check arguments
if [ $# -ne 5 ]; then
    echo "Usage: $0 <source_path> <target_path> <database> <listFile> <num_jobs>"
    exit 1
fi

input_path=$1
output_path=$2
database=$3
listFile=$4
num_jobs=$5

# Function to process a single sample
process_bwa() {
    local F="$1"
    R=${F%_*}_2.fastq.gz
    BASE=${F##*/}
    SAMPLE=${BASE%_*}
    
    echo "Processing: $SAMPLE"

    # Check if sample is in the target list
    if grep -qw "$SAMPLE" "$listFile"; then
        # Resume capability
        if [ -e "$output_path/${SAMPLE}.bam" ]; then
            echo "[SKIP] $SAMPLE BAM already exists."
        else
            # Run bowtie2 and samtools sort
            bowtie2 -p 64 --sensitive -x "$database" -1 "$F" -2 "$R" | samtools sort -@ 14 -o "$output_path/${SAMPLE}.bam"
            echo "[DONE] $SAMPLE alignment finished."
        fi
    else
        echo "[SKIP] $SAMPLE not in list."
    fi
}

processed_files=0

# Run parallel jobs
for F in "$input_path"/*_1.fastq.gz; do
    process_bwa "$F" &
    ((++processed_files))
    [ $((processed_files % num_jobs)) -eq 0 ] && wait
done

# Wait for remaining jobs
wait
echo "All alignment jobs finished!"