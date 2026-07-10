#!/bin/bash

# Max parallel jobs based on original -j 3
MAX_JOBS=3
job_count=0

# Iterate over VIBRANT directories
for dir in 4.vibrant_results/VIBRANT_*; do
    [ -d "$dir" ] || continue

    # Extract sample_id by removing VIBRANT_ prefix
    dirname=$(basename -- "$dir")
    sample_id="${dirname#VIBRANT_}"

    out_dir="4.vibrant_results/VIBRANT_${sample_id}/CheckV_results"
    input_fna="4.vibrant_results/VIBRANT_${sample_id}/VIBRANT_phages_${sample_id}/${sample_id}.phage_combined.fna"

    # Run if output directory does not exist
    if [ ! -d "$out_dir" ]; then
        
        # Run checkV in background
        checkv end_to_end "$input_fna" "$out_dir" -d /public/home/database/checkv-db-v1.5 -t 64 &

        # Concurrency control
        ((job_count++))
        if [ "$job_count" -ge "$MAX_JOBS" ]; then
            wait -n
            ((job_count--))
        fi
    fi
done

# Wait for all remaining background jobs
wait