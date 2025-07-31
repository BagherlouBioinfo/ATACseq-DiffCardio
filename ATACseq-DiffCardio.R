# ATAC-seq Analysis: GSE85330 Dataset (Cardiac Differentiation)
# Author: [Insert Your Name]
# Description: This script performs ATAC-seq analysis comparing chromatin accessibility between Day 0 and Day 30 samples.
# Output: All figures are saved to ./figures, all data files to ./output

# -------------------------
# 1. Setup environment
# -------------------------

setwd("C:/Users/OEM/Desktop/ATACseq-DiffCardio/")

if (!dir.exists("output")) dir.create("output", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# Required CRAN and Bioconductor packages
cran_packages <- c("ggplot2", "enrichplot")
bioc_packages <- c(
  "GenomicRanges", "IRanges", "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
  "clusterProfiler", "DESeq2"
)

# Install missing packages
install_if_missing <- function(pkg, is_bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is_bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

lapply(cran_packages, install_if_missing)
lapply(bioc_packages, install_if_missing, is_bioc = TRUE)

# Load libraries
suppressPackageStartupMessages({
  lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)
})

# -------------------------
# 2. Extract raw data
# -------------------------

raw_tar_file <- "./data/GSE85330_RAW.tar"
extract_dir  <- "./data/GSE85330_RAW"

if (!dir.exists(extract_dir)) dir.create(extract_dir, recursive = TRUE)
untar(tarfile = raw_tar_file, exdir = extract_dir)

bed_files <- list.files(extract_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)
if (length(bed_files) == 0) stop("No BED files found.")

# -------------------------
# 3. Convert first BED file to GRanges
# -------------------------

bed_raw <- read.table(bed_files[1], header = FALSE)

gr <- GRanges(
  seqnames = bed_raw$V1,
  ranges   = IRanges(start = bed_raw$V2 + 1, end = bed_raw$V3),
  strand   = "*",
  score    = bed_raw$V5,
  name     = bed_raw$V4
)

saveRDS(gr, file = "output/ATAC_peaks_GRanges_sample1.rds")

# -------------------------
# 4. Annotate peaks
# -------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak_anno <- annotatePeak(gr, TxDb = txdb, annoDb = "org.Hs.eg.db")

write.csv(as.data.frame(peak_anno), "output/Annotated_Peaks_sample1.csv")

# -------------------------
# 5. Save annotation plots
# -------------------------

pdf("figures/PeakAnnotation_PieChart.pdf")
plotAnnoPie(peak_anno)
dev.off()

pdf("figures/PeakAnnotation_Barplot.pdf")
plotAnnoBar(peak_anno)
dev.off()

pdf("figures/PeakAnnotation_DistToTSS.pdf")
plotDistToTSS(peak_anno)
dev.off()

# -------------------------
# 6. GO enrichment
# -------------------------

genes <- unique(na.omit(as.data.frame(peak_anno)$SYMBOL))

ego <- enrichGO(
  gene          = genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

pdf("figures/GO_Enrichment_Barplot.pdf")
barplot(ego, showCategory = 15)
dev.off()

# -------------------------
# 7. Create peak count matrix
# -------------------------

bed_paths <- list(
  D0_1   = file.path(extract_dir, "GSM2264802_C15_0_1.filterBL.bed.gz"),
  D0_2   = file.path(extract_dir, "GSM2264803_C15_0_2.filterBL.bed.gz"),
  D30_1  = file.path(extract_dir, "GSM2264808_C15_30_1.filterBL.bed.gz"),
  D30_2  = file.path(extract_dir, "GSM2264809_C15_30_2.filterBL.bed.gz")
)

read_bed <- function(path) {
  df <- read.table(path, header = FALSE)
  GRanges(seqnames = df$V1, ranges = IRanges(start = df$V2 + 1, end = df$V3), strand = "*")
}

peak_list <- lapply(bed_paths, read_bed)
all_peaks <- reduce(unlist(GRangesList(peak_list)))
consensus_peaks <- resize(all_peaks, width = 250, fix = "center")

count_matrix <- sapply(peak_list, function(peaks) {
  countOverlaps(consensus_peaks, peaks)
})

rownames(count_matrix) <- paste0("Peak_", seq_len(nrow(count_matrix)))
colnames(count_matrix) <- names(peak_list)

write.csv(count_matrix, "output/Consensus_Peak_CountMatrix.csv")

# -------------------------
# 8. DESeq2 differential analysis
# -------------------------

coldata <- data.frame(
  row.names = colnames(count_matrix),
  condition = c("D0", "D0", "D30", "D30")
)

dds <- DESeqDataSetFromMatrix(count_matrix, coldata, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)

res <- results(dds)
res_df <- as.data.frame(res)
res_df$PeakID <- rownames(res_df)

write.csv(res_df, "output/Differential_Accessibility_Results.csv")

# -------------------------
# 9. Volcano plot
# -------------------------

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Differential Accessibility: Day 30 vs Day 0",
       x = "Log2 Fold Change",
       y = "-Log10 P-value")

ggsave("figures/VolcanoPlot_Day30_vs_Day0.pdf", volcano_plot, width = 7, height = 5)
ggsave("figures/VolcanoPlot_Day30_vs_Day0.png", volcano_plot, width = 7, height = 5)

