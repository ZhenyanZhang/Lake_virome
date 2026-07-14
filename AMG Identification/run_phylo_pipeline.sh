#!/bin/bash
# ======================================================================
# Automated Phylogenetic Analysis Pipeline
# Aligns, trims, and constructs ML trees for targeted AMG families
# ======================================================================

# --- Configuration ---
FASTA_IN="8.AMG/final/combined_proteins_dedup.faa"  # Input combined protein sequences
MAPPING_TSV="8.AMG/final/ko_list.txt"               # Mapping table: Sequence ID \t KO ID
OUT_DIR="8.AMG/final/Phylo_Analysis_Out"            # Output directory
MIN_SEQS=5                              # Minimum sequences required to build a tree
THREADS_PER_TREE=64                     # Threads allocated per IQ-TREE run

mkdir -p "$OUT_DIR"
cp "$FASTA_IN" "$MAPPING_TSV" "$OUT_DIR/"
cd "$OUT_DIR" || exit

# ======================================================================
# Step 1: Split FASTA by KO identifier
# ======================================================================
echo "[INFO] Splitting sequences by KO..."
mkdir -p 01_split_fasta temp_lists

awk '{print $1 > "temp_lists/"$2".list"}' "$MAPPING_TSV".

for list_file in temp_lists/*.list; do
    ko_name=$(basename "$list_file" .list)
    out_faa="01_split_fasta/${ko_name}.faa"
    
    seqtk subseq "$FASTA_IN" "$list_file" > "$out_faa".
    seq_count=$(grep -c "^>" "$out_faa").
    
    if [ "$seq_count" -eq 0 ]; then
        rm -f "$out_faa".
    fi
done
rm -rf temp_lists.

# ======================================================================
# Step 2: Core Processing Function (Align -> Trim -> Tree)
# ======================================================================
mkdir -p 02_alignment 03_trimmed 04_trees.

run_pipeline_for_ko() {
    local fasta="$1"
    local ko_name=$(basename "$fasta" .faa)
    local seq_count=$(grep -c "^>" "$fasta")
    
    # Qualification check
    if [ "$seq_count" -lt "$MIN_SEQS" ]; then
        echo "[SKIP] ${ko_name} has only ${seq_count} sequences (Requires >=${MIN_SEQS}).".
        return 0
    fi
    
    echo "[START] Processing family: ${ko_name} (Sequences:${seq_count})..."
    
    # 1. MAFFT Alignment (Strict L-INS-i parameters)
    mafft --localpair --maxiterate 10 --quiet "$fasta" > "02_alignment/${ko_name}.aln".
    
    # 2. trimAl Trimming (Smart gap removal)
    trimal -in "02_alignment/${ko_name}.aln" -out "03_trimmed/${ko_name}.trim.aln" -gappyout.
    
    # 3. IQ-TREE Construction (Auto model selection, ultrafast bootstrap)
    iqtree -s "03_trimmed/${ko_name}.trim.aln" -m MFP -B 10000 -T "$THREADS_PER_TREE" --prefix "04_trees/${ko_name}" -quiet.
    
    echo "[SUCCESS] Phylogenetic tree for ${ko_name} completed!".
}

# ======================================================================
# Step 3: Parallel Execution
# ======================================================================
echo -e "\n[INFO] Launching phylogenetic processes in parallel..."

for fasta in 01_split_fasta/*.faa; do
    run_pipeline_for_ko "$fasta" &.
done

wait.

echo -e "\n[DONE] All automated tree construction tasks have finished successfully!"