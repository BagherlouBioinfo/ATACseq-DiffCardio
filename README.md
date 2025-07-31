# 📊 ATAC-seq Analysis of Cardiac Differentiation (GSE85330)

This project demonstrates a complete and reproducible ATAC-seq data analysis pipeline using the public dataset **[GSE85330](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85330)**. The dataset captures chromatin accessibility changes in human induced pluripotent stem cells (hiPSCs) during cardiac differentiation (Day 0 vs. Day 30).

> The goal is to identify differentially accessible chromatin regions and their functional implications during lineage specification into cardiomyocytes.

---

## Dataset Summary

| Timepoint | Cell Type                        | Description                              |
|-----------|----------------------------------|------------------------------------------|
| Day 0     | hiPSC (undifferentiated)         | Pluripotent state                        |
| Day 30    | Differentiated cardiomyocytes    | Post-differentiation into heart cells    |

The raw data files are downloaded and extracted from GEO as `.bed.gz` files and processed within R.

---

## Analysis Workflow

The pipeline includes the following steps:

1. **Download and extract raw data** (`GSE85330_RAW.tar`)
2. **Convert BED to GRanges** objects
3. **Peak annotation** using ChIPseeker and UCSC hg38 genome
4. **GO enrichment analysis** for genes near peaks
5. **Build consensus peak matrix** across replicates
6. **Differential accessibility analysis** (DESeq2)
7. **Visualization** (peak annotation plots, volcano plot, GO enrichment)

All plots are saved in `./figures/` and results/data in `./output/`.

---

## Required R Packages

This project uses the following R/Bioconductor packages:

- `GenomicRanges`, `IRanges`
- `ChIPseeker`, `TxDb.Hsapiens.UCSC.hg38.knownGene`
- `org.Hs.eg.db`, `clusterProfiler`
- `DESeq2`, `ggplot2`, `enrichplot`

The script automatically installs missing packages.

---

## Directory Structure

```
ATACseq-DiffCardio/
│
├── data/                  # Downloaded raw data from GEO (BED + tar files)
├── figures/               # All PDF plots (annotation, GO, volcano, etc.)
├── output/                # Processed objects, annotation tables, count matrix
├── ATACseq_pipeline.R     # Full R script
└── README.md              # This file
```

---

## Example Plots

All plots are saved as high-quality PDFs:

- `PeakAnnotation_PieChart.pdf`
- `GO_Enrichment_Barplot.pdf`
- `VolcanoPlot_Day30_vs_Day0.pdf`

---

## Key Results

- Chromatin regions with significantly increased or decreased accessibility between Day 0 and Day 30 were identified.
- GO enrichment highlights key biological processes associated with cardiac differentiation.

---

## How to Run

```r
# Run the full pipeline in R
source("ATACseq_pipeline.R")
```

Ensure you are in the project root directory and have internet access to download GEO data.

---

## Citation

**Original dataset**:  
Zhou et al., *Genome-wide maps of chromatin state in pluripotent and differentiated cells*  
GEO Accession: [GSE85330](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85330)

---

## 👤 Author

**Dr. Nazanin Bagherlou**  
PhD in Biomolecular and Health Sciences
n.bagherlou@campus.uniurb.it
