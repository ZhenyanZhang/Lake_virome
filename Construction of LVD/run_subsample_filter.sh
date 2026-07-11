#!/bin/bash

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <BAM_DIR> <BREADTH_CUTOFF> <OUT_DIR>"
    echo "Example: $0 ./bam_files 0.75 ./Final_Results_bowtie"
    exit 1
fi

INPUT_DIR="$1"
BREADTH_CUTOFF="$2"
OUT_DIR="$3"

# Target depth
TARGET_DEPTH=50000000

# Existing depth log and low-depth sample list
COUNTS_FILE="5.abundance/depth_counts.txt"
SMALL_SAMPLE_FILE="5.abundance/small_sample_75.txt"

mkdir -p "$OUT_DIR"

if [ ! -f "$COUNTS_FILE" ]; then
    echo "Error: $COUNTS_FILE not found!"
    exit 1
fi

echo "Target Depth: $TARGET_DEPTH"
echo "Breadth Cutoff: $BREADTH_CUTOFF"

for bam in "${INPUT_DIR}"/*.bam; do
    [ -e "$bam" ] || continue
    
    SAMPLE=$(basename "$bam" .bam)
    FINAL_STATS="${OUT_DIR}/${SAMPLE}.final_stats.txt"
    
    # Resume capability
    if [ -f "$FINAL_STATS" ]; then
        echo "[Skipping] Output exists for $SAMPLE"
        continue
    fi

    echo "[Processing] $SAMPLE"

    # Get current depth
    CURRENT_DEPTH=$(grep "$bam" "$COUNTS_FILE" | awk '{print $2}')

    if [ -z "$CURRENT_DEPTH" ]; then
        echo "Warning: Depth not found for $bam. Skipping."
        continue
    fi

    # Check depth threshold
    if [ "$CURRENT_DEPTH" -lt "$TARGET_DEPTH" ]; then
        echo "Skipped: Depth $CURRENT_DEPTH <$TARGET_DEPTH"
        # Append skipped sample record
        echo -e "${SAMPLE}\t${CURRENT_DEPTH}" >> "$SMALL_SAMPLE_FILE"
        continue
    fi

    # Calculate subsampling fraction
    FACTOR=$(awk -v t="$TARGET_DEPTH" -v c="$CURRENT_DEPTH" 'BEGIN {printf "%.6f", t/c}')
    SUBSAMPLE_PARAM="100${FACTOR}"
    
    # Subsample and sort
    TEMP_BAM="${OUT_DIR}/${SAMPLE}.temp.bam"
    samtools view -@ 64 -s "$SUBSAMPLE_PARAM" -b "$bam" \vert{} samtools sort -@ 64 - -o "$TEMP_BAM"
    samtools index -@ 64 "$TEMP_BAM"

    # Quantification with CoverM
    COVERM_RAW="${OUT_DIR}/${SAMPLE}.coverm_raw.txt"
    coverm contig \
        --bam-files "$TEMP_BAM" \
        --output-format sparse \
        --methods trimmed_mean covered_fraction count \
        --trim-min 5 --trim-max 95 \
        --min-covered-fraction 0 \
        --output-file "$COVERM_RAW"

    # Filter by breadth cutoff
    awk -v cutoff="$BREADTH_CUTOFF" -F'\t' '
    NR==1 {print $0} 
    NR>1 {
        if ($4 >= cutoff) print $0
    }' "$COVERM_RAW" > "$FINAL_STATS"

    # Clean up temporary files
    rm "$TEMP_BAM" "${TEMP_BAM}.bai" "$COVERM_RAW"
done

echo "Analysis Finished. Skipped samples logged in $SMALL_SAMPLE_FILE"