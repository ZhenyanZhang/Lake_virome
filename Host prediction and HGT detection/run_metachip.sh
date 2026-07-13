#!/usr/bin/env bash
# ===============================================================================
# Host-to-Host HGT Detection Pipeline
# Pre-requisites: Dereplicated MAGs and GTDB-Tk summaries
# ===============================================================================
set -euo pipefail

# ============================ Configuration =========================
MAG_DIR="3.MAG/bin_drep/dereplicated_genomes" 
EXT="fa" 

# GTDB-Tk summaries
GTDB_BAC="3.MAG/drep_bins.gtdbtk/gtdbtk.bac120.summary.tsv"   
GTDB_AR="3.MAG/drep_bins.gtdbtk/gtdbtk.ar53.summary.tsv"      

PREFIX="lakehost"
THREADS=64                  
RANKS="pcofg"  # Taxonomic levels: phylum, class, order, family, genus
OUTROOT="7.host_prediction_and_HGT/hgt_analysis"

METACHIP_ENV="metachip"  # Conda environment for MetaCHIP
# ============================================================================

TAXON_FILE="${OUTROOT}/${PREFIX}_taxon.tsv"
IDS_FILE="${OUTROOT}/${PREFIX}_host_ids.txt"
MC_WD="${OUTROOT}/02_metachip"

run_in_env () {
  local env="$1"; shift
  if [[ -n "$env" ]]; then conda run --no-capture-output -n "$env" "$@"; else "$@"; fi
}

echo "==> Creating output directory"
mkdir -p "$OUTROOT"

# ============================ 1. Preparation ===========================
echo "==> [1] Generating host ID list"
shopt -s nullglob
mags=( "${MAG_DIR}"/*."${EXT}" )
shopt -u nullglob

for f in "${mags[@]}"; do basename "$f" ".${EXT}"; done \vert{} sort -u > "$IDS_FILE"

echo "==> [2] Extracting subset from GTDB-Tk summaries"
summaries=()
[[ -f "$GTDB_BAC" ]] && summaries+=( "$GTDB_BAC" )
[[ -f "$GTDB_AR"  ]] && summaries+=( "$GTDB_AR" )

{
  printf 'user_genome\tclassification\n'
  awk -F'\t' 'NR==FNR{want[$1]=1; next} ($1 in want){print $1"\t"$2}' \
      "$IDS_FILE" "${summaries[@]}"
} > "$TAXON_FILE"

# ============================ 2. MetaCHIP HGT Detection ================
echo "==> [3] MetaCHIP Phase Initialization (PI) - Protein prediction & Species tree"
run_in_env "$METACHIP_ENV" MetaCHIP PI \
  -i "$MAG_DIR" -x "$EXT" -taxon "$TAXON_FILE" \
  -p "$PREFIX" -r "$RANKS" -t "$THREADS" -o "$MC_WD"

echo "==> [4] MetaCHIP Basic Processing (BP) - Best match & Phylogenetic evaluation"
run_in_env "$METACHIP_ENV" MetaCHIP BP \
  -p "$PREFIX" -r "$RANKS" -t "$THREADS" -o "$MC_WD"

echo "Done! Detected HGTs are in ${MC_WD}/${PREFIX}_${RANKS}_detected_HGTs.txt"