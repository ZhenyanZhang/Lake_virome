#!/bin/bash

# Check if the user provided exactly two arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_mag_dir> <gtdbtk_summary_dir>"
    echo "Example: $0 bin_drep classify_wf_out"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR=$1
SUMMARY_DIR=$2

# Define output directories
OUT_BAC="3.MAG/drep_bins_new/bac"
OUT_ARC="3.MAG/drep_bins_new/arc"

# Create output directories (creates parent directories if needed, no error if exists)
mkdir -p "$OUT_BAC"
mkdir -p "$OUT_ARC"

echo "========================================"
echo "Starting MAG separation process..."
echo "Input MAGs directory: $INPUT_DIR"
echo "GTDB-Tk summary directory: $SUMMARY_DIR"
echo "========================================"

# ----------------------------------------
# Process Bacteria (bac120)
# ----------------------------------------
BAC_SUMMARY="${SUMMARY_DIR}/gtdbtk.bac120.summary.tsv"

if [ -f "$BAC_SUMMARY" ]; then
    echo "[Info] Processing bacterial MAGs..."
    
    # Skip the header (NR>1) and extract the first column ($1)
    awk -F'\t' 'NR>1 {print $1}' "$BAC_SUMMARY" | while read -r genome_id; do
        # Construct the full path to the source fasta file
        source_file="${INPUT_DIR}/${genome_id}.fa"
        
        # Check if the file exists before copying
        if [ -f "$source_file" ]; then
            cp "$source_file" "$OUT_BAC/"
        else
            echo "[Warning] File not found: $source_file"
        fi
    done
    echo "[Success] Bacterial MAGs copied to $OUT_BAC"
else
    echo "[Notice] $BAC_SUMMARY not found, skipping bacteria."
fi

echo "----------------------------------------"

# ----------------------------------------
# Process Archaea (ar53)
# ----------------------------------------
ARC_SUMMARY="${SUMMARY_DIR}/gtdbtk.ar53.summary.tsv"

if [ -f "$ARC_SUMMARY" ]; then
    echo "[Info] Processing archaeal MAGs..."
    
    awk -F'\t' 'NR>1 {print $1}' "$ARC_SUMMARY" | while read -r genome_id; do
        source_file="${INPUT_DIR}/${genome_id}.fa"
        
        if [ -f "$source_file" ]; then
            cp "$source_file" "$OUT_ARC/"
        else
            echo "[Warning] File not found: $source_file"
        fi
    done
    echo "[Success] Archaeal MAGs copied to $OUT_ARC"
else
    echo "[Notice] $ARC_SUMMARY not found, skipping archaea."
fi

echo "========================================"
echo "All tasks completed!"