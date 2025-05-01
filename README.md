# golemMB: A Shiny App for 16S rRNA Gene Amplicon Sequence Analysis

**golemMB** is a Shiny application developed to streamline and document 16S rRNA gene amplicon sequencing analysis from the count table onward in a reproducible, modular workflow. It supports data import into `phyloseq`, alpha and beta diversity analysis, differential abundance, functional predictions via PICRUSt2, BugBase phenotype visualization, and correlation analysis between microbial taxa and numeric metadata.

## Features

- Import and manage OTU/ZOTU tables using `phyloseq`
- Alpha diversity calculation with six ecological indices
- Data normalization using multiple transformation methods (e.g., CLR, log1p, RLE)
- Beta diversity analysis with ordinations and distance matrices
- Differential abundance testing using multiple statistical models
- Visualization of functional predictions from **PICRUSt2** and **BugBase**
- Correlation analysis between taxonomic units and numeric metadata
- Session bookmarking and reproducible PDF reports via R Markdown

## Installation

### Option 1: Run with Docker

```bash
git clone https://github.com/yourusername/golemMB.git
docker build -t golemmb .
docker run -p 3838:3838 golemmb
```

### Option 2: Install from Tarball

```r
install.packages("deploy/golemMB_0.1.0.tar.gz", repos = NULL, type = "source")
```

The app is organized into separate tabs corresponding to the previously mentionedanalysis steps. The most important step, the phyloseq import, requires these files in a specific format:

- TU (OTU/ZOTU) count tables (.tab)
- OTU/ZOTU sequences (.fasta)
- Phylogenetic tree (.tree)
- Metadata (.rds)

Non-numeric metadata variables (e.g., treatment groups) must be converted to factors prior to upload using factor() in R.

As analyses go, new input files for subsequent analyses are generated. Instructions for required inputs are shown in each tab. The app includes input validation and will notify users of missing files or parameters before running any analysis.


# PICRUSt2 and BugBase

- BugBase results (predictions.txt and vsearch_stats.txt) must be generated externally (via the BugBase CLI or web tool) and uploaded as input.
- Functional differential abundance requires output from PICRUSt2 (run with the flags below):

```bash
--in_traits EC,KO,COG,PFAM,TIGRFAM
--stratified
```

## Output

- Downloadable PDF reports with rendered R Markdown code and plots
- Result-specific .rds files for chaining between analysis steps
- statistical results and separate plots
