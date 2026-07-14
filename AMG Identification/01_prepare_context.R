#!/usr/bin/env Rscript
# =====================================================================
# Script 1: Data Preparation for Genomic Context Visualization
# ---------------------------------------------------------------------
# Input: annotations.tsv (from DRAM-v) & votu_amg_list.txt (List of identified AMGs)
# Output: context_plot_input.tsv (for downstream plotting)
# =====================================================================

# Set working directory (Change this path if needed)
wd <- "/onedrive/zzy/isme_revised/2.AMG"
setwd(wd)

ANN_PATH  <- "all_annotations.tsv"
LIST_PATH <- "votu_amg_list.txt"        
OUT_PATH  <- "context_plot_input.tsv"

# Read inputs
ann <- read.delim(ANN_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                  check.names=FALSE, quote="", comment.char="", row.names=NULL)
lst <- read.delim(LIST_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                  check.names=FALSE)

# Ensure correct column naming for the gene ID
if (names(ann)[1] %in% c("", "V1", "X", "row.names")) {
  names(ann)[1] <- "gene"
}
gene_col <- "gene"
cat(sprintf("[INFO] Loaded annotations: %d rows\n", nrow(ann)))

# Standardize list columns (votu_id, amg_gene_id, ko_id)
names(lst)[1:3] <- c("votu_id", "amg_gene_id", "ko_id")
cat(sprintf("[INFO] Loaded AMG list: %d genes across %d vOTUs\n",
            nrow(lst), length(unique(lst$votu_id))))

# Build AMG lookup table
amg_ko <- setNames(as.character(lst$ko_id), lst$amg_gene_id)

# Filter annotations to keep only vOTUs present in the AMG list
uniq_votu <- unique(lst$votu_id)
rows_list <- list(); matched <- character(0)
for (v in uniq_votu) {
  idx <- startsWith(ann$scaffold, v)
  if (any(idx)) {
    matched <- c(matched, v)
    sub <- ann[idx,,drop=FALSE]; sub$votu_clean <- v
    rows_list[[v]] <- sub
  }
}

if (length(matched) == 0) stop("No vOTUs matched between the list and annotations.")
dat <- do.call(rbind, rows_list)

# Extract coordinates and strand
dat$start <- suppressWarnings(as.numeric(dat$start_position))
dat$end   <- suppressWarnings(as.numeric(dat$end_position))
sv <- as.character(dat$strandedness)
dat$strand <- ifelse(sv %in% c("-1","-"), "-", "+")

# Mark AMG status and associated KO
dat$is_amg <- dat$gene %in% names(amg_ko)
dat$amg_ko <- ifelse(dat$is_amg, amg_ko[dat$gene], NA)

# Categorize genes (AMG, viral, other_func, hypothetical, unannotated)
has_val <- function(x) !is.na(x) & nzchar(trimws(x)) & x != "NA"
hyp_kw  <- "hypothetical|uncharacterized|unknown function|DUF[0-9]"
viral_kw_extra <- paste0("terminase|portal|capsid|tail|baseplate|coat|head|",
  "prohead|holin|lysin|integrase|spike|packaging|tape measure|sheath|virion|phage|viral")

classify <- function(r) {
  if (r$is_amg) return("AMG")
  
  if (has_val(r$vogdb_hits) \vert{}\vert{} has_val(r$vogdb_id) ||
      has_val(r$viral_hit)  \vert{}\vert{} has_val(r$viral_id) ||
      grepl(viral_kw_extra, paste(r$kegg_hit, r$vogdb_categories, r$vogdb_hits), ignore.case=TRUE)) {
    return("viral")
  }
  
  if (has_val(r$ko_id) || has_val(r$kegg_hit) \vert{}\vert{} has_val(r$pfam_hits) ||
      has_val(r$cazy_hits) \vert{}\vert{} has_val(r$peptidase_hit)) {
    return("other_func")
  }
  
  if (grepl(hyp_kw, paste(r$kegg_hit, r$pfam_hits), ignore.case=TRUE)) {
    return("hypothetical")
  }
  
  return("unannotated")
}

dat$category <- vapply(seq_len(nrow(dat)), function(i) classify(dat[i,,drop=FALSE]), character(1))
dat$label    <- ifelse(dat$is_amg, dat$amg_ko, NA_character_)

# Calculate distance to contig ends (for optional quality checks)
scaf_min <- tapply(dat$start, dat$scaffold, min, na.rm=TRUE)
scaf_max <- tapply(dat$end,   dat$scaffold, max, na.rm=TRUE)
dat$dist_to_end <- NA_real_
for (i in seq_len(nrow(dat))) {
  if (dat$is_amg[i]) {
    s <- dat$scaffold[i]
    dat$dist_to_end[i] <- min(dat$start[i] - scaf_min[[s]], scaf_max[[s]] - dat$end[i], na.rm=TRUE)
  }
}

# Output formatted data
out <- data.frame(
  votu = dat$votu_clean, scaffold = dat$scaffold, gene = dat$gene,
  start = dat$start, end = dat$end, strand = dat$strand, 
  category = dat$category, amg_ko = dat$amg_ko, label = dat$label,
  dist_to_end = dat$dist_to_end, stringsAsFactors = FALSE)

out <- out[order(out$votu, out$start), ]
write.table(out, OUT_PATH, sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("\n[SUCCESS] %s generated (%d genes, %d vOTUs)\n", OUT_PATH, nrow(out), length(unique(out$votu))))
