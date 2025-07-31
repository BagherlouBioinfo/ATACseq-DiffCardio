ATAC-seq Analysis Pipeline for GSE85330 Dataset
================
Dr. Nazanin Bagherlou
2025-07-27

# **Introduction**

In this report, I demonstrate my ability to perform a comprehensive
ATAC-seq data analysis using the publicly available GSE85330 dataset.
This dataset includes ATAC-seq profiles of human induced pluripotent
stem cells (hiPSCs) at two critical timepoints:

Day 0 (undifferentiated state)

Day 30 (differentiated cardiomyocytes)

The main objective of this analysis is to identify and interpret changes
in chromatin accessibility between these two conditions, which reflect
the epigenomic dynamics associated with cardiac lineage differentiation.

To achieve this, I implemented a complete and reproducible pipeline
using R and Bioconductor, including:

- Preprocessing and reading raw BED files from GEO

- Converting peaks to GRanges objects

- Peak annotation with known gene features

- Functional enrichment analysis (GO)

- Consensus peak matrix generation

- Differential accessibility analysis using DESeq2

- Visualization of genomic distributions and statistical results

This project serves as an example of my proficiency in epigenomic data
analysis, reproducible coding, and interpretation of chromatin
accessibility landscapes using real-world biological datasets. All steps
are fully automated and documented to facilitate reproducibility and
reuse.

------------------------------------------------------------------------

# 1. Install and Load Required Packages

``` r
cran_packages <- c("ggplot2", "enrichplot")
bioc_packages <- c(
  "GenomicRanges", "IRanges", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db", "clusterProfiler", "DESeq2"
)

install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

invisible(lapply(cran_packages, install_if_missing_cran))
invisible(lapply(bioc_packages, install_if_missing_bioc))

suppressPackageStartupMessages({
  lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)
})
```

    ## [[1]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[2]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[3]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[4]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[5]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[6]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[7]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[8]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"                             
    ## 
    ## [[9]]
    ##  [1] "DESeq2"                            "SummarizedExperiment"             
    ##  [3] "MatrixGenerics"                    "matrixStats"                      
    ##  [5] "clusterProfiler"                   "org.Hs.eg.db"                     
    ##  [7] "TxDb.Hsapiens.UCSC.hg38.knownGene" "GenomicFeatures"                  
    ##  [9] "AnnotationDbi"                     "Biobase"                          
    ## [11] "ChIPseeker"                        "GenomicRanges"                    
    ## [13] "GenomeInfoDb"                      "IRanges"                          
    ## [15] "S4Vectors"                         "BiocGenerics"                     
    ## [17] "stats4"                            "enrichplot"                       
    ## [19] "ggplot2"                           "stats"                            
    ## [21] "graphics"                          "grDevices"                        
    ## [23] "utils"                             "datasets"                         
    ## [25] "methods"                           "base"

# 2. Extract Raw Data (.tar)

We first extract the compressed .tar file containing the raw ATAC-seq
BED files.

``` r
raw_tar_file <- "C:/users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW.tar"
extract_dir  <- "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW"

if (!dir.exists(extract_dir)) dir.create(extract_dir, recursive = TRUE)
untar(tarfile = raw_tar_file, exdir = extract_dir)

bed_files <- list.files(extract_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)
bed_files
```

    ##  [1] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264802_C15_0_1.filterBL.bed.gz" 
    ##  [2] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264803_C15_0_2.filterBL.bed.gz" 
    ##  [3] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264804_C15_2_1.filterBL.bed.gz" 
    ##  [4] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264805_C15_2_2.filterBL.bed.gz" 
    ##  [5] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264806_C15_4_1.filterBL.bed.gz" 
    ##  [6] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264807_C15_4_2.filterBL.bed.gz" 
    ##  [7] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264808_C15_30_1.filterBL.bed.gz"
    ##  [8] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264809_C15_30_2.filterBL.bed.gz"
    ##  [9] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264810_C20_0_1.filterBL.bed.gz" 
    ## [10] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264811_C20_0_2.filterBL.bed.gz" 
    ## [11] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264812_C20_2_1.filterBL.bed.gz" 
    ## [12] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264813_C20_2_2.filterBL.bed.gz" 
    ## [13] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264814_C20_4_1.filterBL.bed.gz" 
    ## [14] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264815_C20_4_2.filterBL.bed.gz" 
    ## [15] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264816_C20_30_1.filterBL.bed.gz"
    ## [16] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264817_C20_30_2.filterBL.bed.gz"
    ## [17] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264818_H1_0_1.filterBL.bed.gz"  
    ## [18] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264819_H1_0_2.filterBL.bed.gz"  
    ## [19] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264820_H1_2_1.filterBL.bed.gz"  
    ## [20] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264821_H1_2_2.filterBL.bed.gz"  
    ## [21] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264822_H1_4_1.filterBL.bed.gz"  
    ## [22] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264823_H1_4_2.filterBL.bed.gz"  
    ## [23] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264824_H1_30_1.filterBL.bed.gz" 
    ## [24] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264825_H1_30_2.filterBL.bed.gz" 
    ## [25] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264826_H9_0_1.filterBL.bed.gz"  
    ## [26] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264827_H9_0_2.filterBL.bed.gz"  
    ## [27] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264828_H9_2_1.filterBL.bed.gz"  
    ## [28] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264829_H9_2_2.filterBL.bed.gz"  
    ## [29] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264830_H9_4_1.filterBL.bed.gz"  
    ## [30] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264831_H9_4_2.filterBL.bed.gz"  
    ## [31] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264832_H9_30_1.filterBL.bed.gz" 
    ## [32] "C:/Users/OEM/Desktop/ATACseq-DiffCardio/data/GSE85330_RAW/GSM2264833_H9_30_2.filterBL.bed.gz"

# 3. Convert First Sample to GRanges

We now load the first BED file as an example and convert it into a
GRanges object.

``` r
bed_raw <- read.table(bed_files[1], header = FALSE)
gr <- GRanges(
  seqnames = bed_raw$V1,
  ranges   = IRanges(start = bed_raw$V2 + 1, end = bed_raw$V3),
  strand   = "*",
  score    = bed_raw$V5,
  name     = bed_raw$V4
)

output_dir <- "C:/Users/OEM/Desktop/ATACseq-DiffCardio/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
saveRDS(gr, file = file.path(output_dir, "ATAC_peaks_GRanges_sample1.rds"))
```

# 4. Peak Annotation

We annotate peaks using ChIPseeker and the human genome reference
(hg38).

``` r
peak_gr <- readRDS(file.path(output_dir, "ATAC_peaks_GRanges_sample1.rds"))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_anno <- annotatePeak(peak_gr, TxDb = txdb, annoDb = "org.Hs.eg.db")
```

    ## >> preparing features information...      2025-07-31 8:12:38 AM 
    ## >> identifying nearest features...        2025-07-31 8:12:38 AM 
    ## >> calculating distance from peak to TSS...   2025-07-31 8:12:44 AM 
    ## >> assigning genomic annotation...        2025-07-31 8:12:44 AM 
    ## >> adding gene annotation...          2025-07-31 8:12:55 AM

    ## 'select()' returned 1:many mapping between keys and columns

    ## >> assigning chromosome lengths           2025-07-31 8:12:55 AM 
    ## >> done...                    2025-07-31 8:12:55 AM

``` r
head(as.data.frame(peak_anno))
```

# Annotation Plots

``` r
plotAnnoPie(peak_anno)
```

![](README_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
plotAnnoBar(peak_anno)
```

![](README_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
plotDistToTSS(peak_anno)
```

![](README_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->

# 5. GO Enrichment Analysis

We test whether annotated genes are enriched in specific biological
processes.

``` r
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

barplot(ego, showCategory = 15)
```

![](README_files/figure-gfm/unnamed-chunk-29-1.png)<!-- --> \# 6. Build
Consensus Peak Matrix

``` r
bed_paths <- list(
  D0_1   = file.path(extract_dir, "GSM2264802_C15_0_1.filterBL.bed.gz"),
  D0_2   = file.path(extract_dir, "GSM2264803_C15_0_2.filterBL.bed.gz"),
  D30_1  = file.path(extract_dir, "GSM2264808_C15_30_1.filterBL.bed.gz"),
  D30_2  = file.path(extract_dir, "GSM2264809_C15_30_2.filterBL.bed.gz")
)

read_bed <- function(path) {
  df <- read.table(path, header = FALSE)
  GRanges(seqnames = df$V1,
          ranges   = IRanges(start = df$V2 + 1, end = df$V3),
          strand   = "*")
}

peak_list <- lapply(bed_paths, read_bed)
all_peaks <- GenomicRanges::reduce(unlist(GRangesList(peak_list)))
consensus_peaks <- resize(all_peaks, width = 250, fix = "center")

count_matrix <- sapply(peak_list, function(peaks) {
  countOverlaps(consensus_peaks, peaks)
})

rownames(count_matrix) <- paste0("Peak_", seq_len(nrow(count_matrix)))
colnames(count_matrix) <- names(peak_list)
```

# 8. Differential Accessibility Analysis (Day 0 vs Day 30)

We identify peaks that change in accessibility between Day 0 and Day 30.

``` r
coldata <- data.frame(
  row.names = colnames(count_matrix),
  condition = c("D0", "D0", "D30", "D30")
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = ~ condition
)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are
    ## characters, converting to factors

``` r
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- nbinomWaldTest(dds)

res <- results(dds)
head(res)
```

    ## log2 fold change (MLE): condition D30 vs D0 
    ## Wald test p-value: condition D30 vs D0 
    ## DataFrame with 6 rows and 6 columns
    ##         baseMean log2FoldChange     lfcSE         stat    pvalue      padj
    ##        <numeric>      <numeric> <numeric>    <numeric> <numeric> <numeric>
    ## Peak_1      1.00    4.65098e-16   1.44269  3.22382e-16  1.000000  1.000000
    ## Peak_2      0.50   -2.44269e+00   1.76693 -1.38245e+00  0.166834  0.509897
    ## Peak_3      0.25   -1.44269e+00   2.04027 -7.07109e-01  0.479499  0.617500
    ## Peak_4      0.25   -1.44269e+00   2.04027 -7.07109e-01  0.479499  0.617500
    ## Peak_5      0.25   -1.44269e+00   2.04027 -7.07109e-01  0.479499  0.617500
    ## Peak_6      0.25   -1.44269e+00   2.04027 -7.07109e-01  0.479499  0.617500

# 8. Volcano Plot

``` r
res_df <- as.data.frame(res)
res_df$PeakID <- rownames(res_df)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Differential Accessibility: Day 30 vs Day 0",
    x     = "Log2 Fold Change",
    y     = "-Log10 P-value"
  )
```

![](README_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

# Conclusion

This pipeline identifies genomic regions whose accessibility changes
during cardiomyocyte differentiation. Further integration with gene
expression (RNA-seq) data could enhance biological insights.
