#!/usr/bin/env bash
# ===============================================================================
# Identification of Virus-Carried Host Genes
# ===============================================================================
set -euo pipefail

# ============================ Configuration ============================
VOTU_FASTA_ALL="4.vibrant_results/vOTUs.fna"
IPHOP_RESULT="7.host_prediction_and_HGT/iPHoP/Host_prediction_to_genome_m75.csv"
CHECKV_SUMMARY="4.vibrant_results/all.quality_summary.tsv"
HOST_FAA="7.host_prediction_and_HGT/hgt_analysis/02_metachip/lakehost_pcofg_prodigal_output"
OUTDIR="7.host_prediction_and_HGT/virus_carried_genes"

THREADS=64
EGGNOG_DB="db/eggnog"   

ENV_TOOLS="kraken2"   
EGGNOG_ENV="emapper"    

# DIAMOND thresholds
EVALUE="1e-5"
MIN_ID=50
QCOV=70
SCOV=70
# =========================================================================

run_in_env(){ local e="$1"; shift; if [[ -n "$e" ]]; then conda run --no-capture-output -n "$e" "$@"; else "$@"; fi; }
mkdir -p "$OUTDIR"

# ---------- 1. Extract Assigned vOTUs ----------
echo "==> [1] Extracting vOTUs with identified hosts from iPHoP results"
# Assuming standard iPHoP CSV format (Virus column is first)
awk -F',' 'NR>1 {print $1}' "$IPHOP_RESULT" \vert{} sort -u > "$OUTDIR/host_assigned_vOTUs.list"

run_in_env "$ENV_TOOLS" seqtk subseq "$VOTU_FASTA_ALL" "$OUTDIR/host_assigned_vOTUs.list" > "$OUTDIR/vOTU_by_host.fna"
echo "    Extracted $(wc -l < "$OUTDIR/host_assigned_vOTUs.list") vOTUs with assigned hosts."

# ---------- 2. Predict vOTU proteins ----------
echo "==> [2] Prodigal protein prediction for targeted vOTUs"
run_in_env "$ENV_TOOLS" prodigal -i "$OUTDIR/vOTU_by_host.fna" -p meta \
  -a "$OUTDIR/votu_proteins.faa" -d "$OUTDIR/votu_genes.fna" \
  -o "$OUTDIR/votu_prodigal.gff" -f gff -q

# ---------- 3. CD-HIT clustering of host proteins ----------
echo "==> [3] Clustering MetaCHIP host proteins (95% AAI)"
if [[ -d "$HOST_FAA" ]]; then
  cat "$HOST_FAA"/*.faa > "$OUTDIR/host_all.faa"
else
  cp "$HOST_FAA" "$OUTDIR/host_all.faa"
fi
run_in_env "$ENV_TOOLS" cd-hit -i "$OUTDIR/host_all.faa" -o "$OUTDIR/host_catalog.faa" \
  -c 0.95 -n 5 -d 0 -M 0 -T "$THREADS" >/dev/null

# ---------- 4. DIAMOND alignment ----------
echo "==> [4] DIAMOND alignment (vOTU vs. Host Catalog)"
run_in_env "$ENV_TOOLS" diamond makedb --in "$OUTDIR/host_catalog.faa" -d "$OUTDIR/host_catalog" --quiet
run_in_env "$ENV_TOOLS" diamond blastp \
  -q "$OUTDIR/votu_proteins.faa" -d "$OUTDIR/host_catalog" \
  -o "$OUTDIR/votu_vs_host.tsv" \
  --evalue "$EVALUE" --id "$MIN_ID" --query-cover "$QCOV" --subject-cover "$SCOV" \
  --max-target-seqs 1 -p "$THREADS" \
  -f 6 qseqid sseqid pident length evalue qcovhsp scovhsp
  
cut -f1 "$OUTDIR/votu_vs_host.tsv" \vert{} sort -u > "$OUTDIR/carried_genes.txt"

# Extract carried protein sequences
awk -v LIST="$OUTDIR/carried_genes.txt" \
  'BEGIN{while((getline l<LIST)>0) keep[l]=1}
   /^>/{id=substr($1,2); pr=(id in keep)} pr' \
  "$OUTDIR/votu_proteins.faa" > "$OUTDIR/carried.faa"

# ---------- 5. eggNOG-mapper annotation ----------
echo "==> [5] Annotating with eggNOG-mapper"
run_in_env "$EGGNOG_ENV" emapper.py -i "$OUTDIR/carried.faa" -o votu \
  --output_dir "$OUTDIR" -m diamond --data_dir "$EGGNOG_DB" --cpu "$THREADS" --override

# ---------- 6. Filter false positives ----------
echo "==> [6] Filtering false positives (COG-J & CheckV boundaries)"

# A. Filter core informational genes (COG category J)
awk -F'\t' '
  /^##/      {next}
  /^#query/  {for(i=1;i<=NF;i++) if($i=="COG_category") c=i; next}
  c && index($c,"J")>0 {print $1}
' "$OUTDIR/votu.emapper.annotations" \vert{} sort -u > "$OUTDIR/core_J_genes.txt"

# B. Filter genes overlapping CheckV host boundaries
python3 -c "
import sys
genes_to_remove = set()

# Parse CheckV summary for provirus coordinates
coords = {}
with open('$CHECKV_SUMMARY') as f_chk:
    header = f_chk.readline().strip().split('\t')
    try:
        id_idx, prov_idx, coord_idx = header.index('contig_id'), header.index('provirus'), header.index('provirus_coords')
        for line in f_chk:
            parts = line.strip().split('\t')
            if len(parts) > coord_idx and parts[prov_idx] == 'Yes':
                try:
                    start, end = map(int, parts[coord_idx].split('-'))
                    coords[parts[id_idx]] = (start, end)
                except:
                    pass
    except ValueError:
        pass

# Compare with Prodigal GFF
with open('$OUTDIR/votu_prodigal.gff') as f_gff:
    for line in f_gff:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        contig = parts[0]
        gene_id = contig + '_' + parts[8].split(';')[0].split('=')[1].split('_')[-1]
        
        if contig in coords:
            g_start, g_end = int(parts[3]), int(parts[4])
            v_start, v_end = coords[contig]
            # Flag gene if it falls completely or partially outside the viral boundary
            if g_start < v_start or g_end > v_end:
                genes_to_remove.add(gene_id)

with open('$OUTDIR/checkv_host_genes.txt', 'w') as out:
    for g in genes_to_remove:
        out.write(g + '\n')
"

cat "$OUTDIR/core_J_genes.txt" "$OUTDIR/checkv_host_genes.txt" \vert{} sort -u > "$OUTDIR/genes_to_remove.txt"
comm -23 "$OUTDIR/carried_genes.txt" "$OUTDIR/genes_to_remove.txt" > "$OUTDIR/virus_carried_genes_final.txt"

# Generate final TSV
awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' \
  "$OUTDIR/virus_carried_genes_final.txt" "$OUTDIR/votu_vs_host.tsv" \
  > "$OUTDIR/virus_carried_genes_final.tsv"

echo "Done! Final carried genes saved to ${OUTDIR}/virus_carried_genes_final.tsv"