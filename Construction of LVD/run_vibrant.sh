#!/bin/bash

# Check args
if [ "$#" -ne 2 ]; then
    echo "Usage: bash $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

# Max parallel jobs
MAX_JOBS=3 

# Create output dir
mkdir -p "$OUTPUT_DIR"

job_count=0

# Iterate over .fa files
for input_file in "$INPUT_DIR"/*.fa; do
    [ -e "$input_file" ] || continue

    # Get sample name
    filename=$(basename -- "$input_file")
    sample="${filename%.fa}"
    sample_outdir="$OUTPUT_DIR/VIBRANT_$sample"

    # Skip if output directory already exists
    if [ -d "$sample_outdir" ]; then
        echo "[Skip] Directory exists for: $sample"
        continue
    fi

    echo "[Start] Submitting: $sample ..."

    # Run VIBRANT in background
    python /public/miniconda3/envs/vibrant/bin/VIBRANT_run.py \
        -i "$input_file" \
        -folder "$sample_outdir" \
        -t 64 \
        -d /public/miniconda3/envs/vibrant/share/vibrant-1.2.1/db/databases \
        -m /public/miniconda3/envs/vibrant/share/vibrant-1.2.1/db/files &

    # Concurrency control
    ((job_count++))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then
        wait -n  
        ((job_count--))
    fi
done

# Wait for all remaining background jobs
wait
echo "All VIBRANT jobs finished!"