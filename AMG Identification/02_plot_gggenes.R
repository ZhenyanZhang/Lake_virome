#!/usr/bin/env Rscript
# =====================================================================
# Script 2: Visualizing AMG Genomic Contexts
# ---------------------------------------------------------------------
# Generates genome organization diagrams for each AMG-carrying vOTU.
# Dependencies: ggplot2, gggenes 
# =====================================================================
suppressMessages({ library(ggplot2); library(gggenes) })

wd <- "/onedrive/zzy/isme_revised/3.AMG"
setwd(wd)

IN_PATH <- "context_plot_input.tsv"
OUT_DIR <- "context_plots"
dir.create(OUT_DIR, showWarnings=FALSE)

df <- read.delim(IN_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
df$forward <- df$strand == "+"

# Color Palette and Labels
cat_colors <- c(
  AMG          = "#E64B35",  
  viral        = "#3C5488",  
  other_func   = "#00A087",  
  hypothetical = "#B0B0B0",  
  unannotated  = "#E0E0E0"   
)
cat_labels <- c(
  AMG="AMG", viral="Viral-related",
  other_func="Other functional", hypothetical="Hypothetical protein",
  unannotated="Unannotated")

cat_levels <- c("AMG","viral","other_func","hypothetical","unannotated")
df$category <- factor(df$category, levels=cat_levels)

votus <- unique(df$votu)
cat(sprintf("[INFO] Rendering plots for %d vOTUs\n", length(votus)))

for (v in votus) {
  sub <- df[df$votu == v, , drop=FALSE]
  sub$category <- droplevels(sub$category) 
  lab <- sub[!is.na(sub$label) & nzchar(sub$label), , drop=FALSE] 
  
  p <- ggplot(sub, aes(xmin=start, xmax=end, y=votu, fill=category, forward=forward)) +
    geom_gene_arrow(arrowhead_height=unit(4,"mm"), arrow_body_height=unit(3,"mm")) +
    scale_fill_manual(values=cat_colors, labels=cat_labels, drop=TRUE, name=NULL) +
    theme_genes() +
    labs(title=v, x="Position (bp)", y=NULL) +
    theme(legend.position="bottom",
          plot.title=element_text(size=9, face="bold"),
          axis.text.y=element_blank())
  
  # Annotate only AMGs with their KO labels
  if (nrow(lab) > 0) {
    p <- p + geom_text(data=lab, aes(x=(start+end)/2, y=votu, label=label),
                       inherit.aes=FALSE, vjust=-1.6, size=2.6,
                       fontface="bold", color="#B22222")
  }
  
  ko_all <- sub$amg_ko[!is.na(sub$amg_ko) & nzchar(sub$amg_ko)]
  ko_tag <- if (length(ko_all)>0) paste(sort(unique(ko_all)), collapse="_") else "noKO"
  
  safe_votu <- gsub("[^A-Za-z0-9_.-]", "_", v)
  fname <- paste0(safe_votu, "__", ko_tag)
  
  ggsave(file.path(OUT_DIR, paste0(fname,".pdf")), p, width=11, height=2.6, limitsize=FALSE)
}
cat(sprintf("[SUCCESS] All plots saved to %s/\n", OUT_DIR))
