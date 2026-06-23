# 🧬 16S rRNA Microbiome Analysis Pipeline

## Bronchoalveolar Lavage Fluid (BALF) Microbiome in Lung Cancer With and Without COPD

[![QIIME2 Version](https://img.shields.io/badge/QIIME2-2024.10-brightgreen)](https://qiime2.org/)
[![R Version](https://img.shields.io/badge/R-4.3%2B-blue)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive, reproducible bioinformatics pipeline for 16S rRNA gene sequencing analysis of bronchoalveolar lavage fluid (BALF) microbiome, comparing three clinical groups: **benign lung mass controls**, **lung cancer patients with COPD (LC_COPD)**, and **lung cancer patients without COPD (LC_woCOPD)**.

---

## 📋 Table of Contents

1. [Project Overview](#-project-overview)
2. [Study Design](#-study-design)
3. [Requirements](#-requirements)
4. [Installation & Setup](#-installation--setup)
5. [Directory Structure](#-directory-structure)
6. [Pipeline Workflow](#-pipeline-workflow)
7. [Step-by-Step Execution Guide](#-step-by-step-execution-guide)
8. [R Analysis Scripts](#-r-analysis-scripts)
9. [Key Results Summary](#-key-results-summary)
10. [Output Files](#-output-files)
11. [Parameters & Configuration](#-parameters--configuration)
12. [Reproducing the Analysis](#-reproducing-the-analysis)
13. [Troubleshooting](#-troubleshooting)
14. [Citation](#-citation)
15. [Contact](#-contact)

---

## 🔬 Project Overview

This pipeline performs end-to-end 16S rRNA microbiome analysis from raw Illumina paired-end FASTQ files through to publication-quality figures. The analysis investigates microbial community differences in bronchoalveolar lavage fluid (BALF) among lung cancer patients stratified by COPD comorbidity.

### Analysis Modules

| Module | Tool(s) | Description |
|--------|---------|-------------|
| **Quality Control** | FastQC, MultiQC | Pre- and post-trimming quality assessment |
| **Primer Trimming** | Cutadapt | V3-V4 primer removal (341F / 805R+806R) |
| **Denoising** | DADA2 | Amplicon Sequence Variant (ASV) inference |
| **Taxonomy** | SILVA 138 (primary), GTDB r214 (optional) | Taxonomic classification |
| **Phylogeny** | MAFFT + FastTree | Phylogenetic tree construction |
| **Alpha Diversity** | phyloseq, vegan | Within-sample diversity (8 metrics) |
| **Beta Diversity** | phyloseq, vegan | Between-sample diversity (4 distance metrics + PCoA/NMDS) |
| **Composition** | phyloseq, ggplot2 | Phylum/Genus-level barplots, heatmaps, core microbiome, Venn diagrams |
| **Differential Abundance** | ANCOM-BC2 | Log-fold change analysis with bias correction |
| **LEfSe** | microbiomeMarker | Linear discriminant analysis Effect Size for biomarker discovery |
| **Functional Prediction** | PICRUSt2 | MetaCyc pathways, KO, and EC number predictions |
| **Visualization** | ggplot2, ComplexHeatmap | Publication-quality figures (PDF, PNG, TIFF) |

---

## 🧪 Study Design

### Sample Groups

| Group | Samples | Condition | COPD Status | Sample Type |
|-------|---------|-----------|-------------|-------------|
| **Control** | C1–C10 (n=10) | Healthy | No | BALF |
| **LC_COPD** | LC1–LC10 (n=10) | Lung Cancer | Yes | BALF |
| **LC_woCOPD** | LC11–LC20 (n=10) | Lung Cancer | No | BALF |

### Sequencing Details

- **Platform**: Illumina MiSeq
- **Read type**: Paired-end (2 × 300 bp)
- **Target region**: V3-V4 hypervariable region
- **Forward primer**: 341F (`CCTACGGGNGGCWGCAG`, 17 bp)
- **Reverse primers**: Mixed — 805R (`GACTACHVGGGTATCTAATCC`, 21 bp, ~39% detection) and 806R (`GGACTACHVGGGTWTCTAAT`, 20 bp, ~89% detection)
- **Total raw files**: 60 FASTQ files (30 samples × 2 read directions)

---

## 💻 Requirements

### Hardware

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 4 cores | 6+ cores |
| RAM | 8 GB | 11+ GB (with swap enabled) |
| Swap | 8 GB | 16+ GB |
| Storage | 30 GB | 50+ GB free |

> **WSL users:** This pipeline was developed on WSL with **11 GB RAM** and **28 GB swap** allocated. DADA2 denoising and PICRUSt2 are the most memory-intensive steps. If running under WSL, ensure adequate RAM and swap are configured in your `.wslconfig` file:
> ```ini
> # %UserProfile%\.wslconfig
> [wsl2]
> memory=11GB
> swap=28GB
> processors=6
> ```

### Software

- **OS**: Linux or Windows Subsystem for Linux (WSL)
- **Package manager**: Conda or Mamba
- **R**: version 4.3+
- **RStudio** (optional): For interactive R analysis

### Conda Environments

The pipeline uses three separate conda environments to avoid dependency conflicts:

| Environment | Purpose | Key Packages |
|------------|---------|--------------|
| `qiime2-amplicon-2024.10` | QIIME2 processing | QIIME2, DADA2, Cutadapt, FastQC, MultiQC |
| `microbiome-r` | R statistical analysis | R, phyloseq, ANCOMBC, vegan, ggplot2, ComplexHeatmap |
| `picrust2` | Functional prediction | PICRUSt2 |

Environment YAML specs are available in the `envs/` directory.

---

## 🛠 Installation & Setup

### Step 1: Clone the Repository

```bash
git clone https://github.com/<your-username>/LC_COPD_microbiome.git
cd LC_COPD_microbiome
```

### Step 2: Run the Automated Setup Script

```bash
chmod +x scripts/*.sh
bash scripts/00_setup_environment.sh
```

This script will:
- Create the full directory structure (see [Directory Structure](#-directory-structure))
- Set up all three conda environments
- Download the **SILVA 138** classifier (~1 GB)
- Optionally prepare the **GTDB r214** database (requires RESCRIPt training or a pre-trained artifact)

### Step 3: Install Additional R Packages

```bash
conda activate microbiome-r
Rscript scripts/install_r_packages.R
```

This installs Bioconductor and CRAN packages including: `phyloseq`, `ANCOMBC`, `microbiomeMarker`, `ComplexHeatmap`, `tidyverse`, `vegan`, `ape`, `VennDiagram`, and the project's custom plotting theme.

### Step 4: Restore R Package Versions (renv)

This project uses [`renv`](https://rstudio.github.io/renv/) to lock exact R package versions for reproducibility. The `renv.lock` file (included in the repository) records every package and its version. The `renv/` library directory itself is **not included in the repository** because it is generated locally.

To recreate the exact R environment used in this analysis:

```bash
conda activate microbiome-r
cd LC_COPD_microbiome

# 1. Install renv if it is not already available
Rscript --vanilla -e 'install.packages("renv", repos="https://cloud.r-project.org")'

# 2. Restore all package versions recorded in renv.lock
#    This will create a local renv/ folder and install all packages
#    (CRAN + Bioconductor) with the exact versions used in the project
Rscript --vanilla -e 'renv::restore()'
```

> **Note:** The first `renv::restore()` may take **15–30 minutes** because it downloads and installs all dependencies (including Bioconductor packages such as `phyloseq`, `ANCOMBC`, `ComplexHeatmap`, etc.).
>
> All packages are installed into a **project-local library (`renv/library/`)**, leaving your global R environment unchanged.

If `renv::restore()` fails due to missing system libraries, you can fall back to manual installation:

```bash
Rscript scripts/install_r_packages.R
```

### Step 5: Place Raw Data

Copy your raw FASTQ files into the `raw_data/` directory:

```
raw_data/
├── C1_R1.fastq.gz
├── C1_R2.fastq.gz
├── C2_R1.fastq.gz
├── ...
└── LC20_R2.fastq.gz
```

File naming convention: `{SampleID}_R{1|2}.fastq.gz`

---

## 📁 Directory Structure

```
LC_COPD_microbiome/
├── raw_data/                          # Original FASTQ files (read-only, 60 files)
├── scripts/                           # Bash pipeline scripts
│   ├── config.sh                      # Central configuration (paths, parameters)
│   ├── 00_setup_environment.sh        # Environment & directory setup
│   ├── 01_fastqc_pretrim.sh           # Pre-trimming quality check
│   ├── 02_trim_primers.sh             # Cutadapt primer removal
│   ├── 03_fastqc_posttrim.sh          # Post-trimming quality check
│   ├── 04_qiime2_import.sh            # Import FASTQs into QIIME2
│   ├── 05_dada2_denoise.sh            # DADA2 denoising (ASV inference)
│   ├── 06_taxonomy_silva.sh           # SILVA 138 classification
│   ├── 07_taxonomy_gtdb.sh            # GTDB r214 classification (optional)
│   ├── 08_phylogeny.sh                # MAFFT alignment + FastTree
│   ├── 09_diversity.sh                # Alpha & beta diversity (QIIME2)
│   ├── 10_picrust2.sh                 # PICRUSt2 functional prediction
│   ├── install_r_packages.R           # R package installer
│   └── optimize_dada2_params.R        # DADA2 parameter grid-search
├── R_analysis/                        # R analysis scripts
│   ├── 01_import_phyloseq.R           # Import QIIME2 → phyloseq
│   ├── 02_alpha_diversity.R           # Alpha diversity analysis
│   ├── 03_beta_diversity.R            # Beta diversity & ordination
│   ├── 04_composition.R               # Taxonomic composition & core microbiome
│   ├── 05_ancombc2.R                  # ANCOM-BC2 differential abundance
│   ├── 06_picrust_visualization.R     # PICRUSt2 pathway visualization
│   ├── 06b_ko_visualization.R         # KEGG Ortholog (KO) visualization
│   ├── 06c_ec_visualization.R         # Enzyme Commission (EC) visualization
│   ├── 07_publication_plots.R         # Publication-quality composite figures
│   ├── 08_lefse.R                     # LEfSe biomarker analysis
│   ├── functions/
│   │   └── plotting_theme.R           # Shared publication plotting theme
│   └── data/                          # Intermediate R data objects (.rds)
├── results/
│   ├── qc_reports/                    # FastQC & MultiQC outputs
│   ├── trimmed_reads/                 # Primer-trimmed sequences
│   ├── qiime2/                        # QIIME2 artifacts (.qza / .qzv)
│   │   ├── imported/                  # Imported sequences
│   │   ├── dada2_output/              # ASV table & rep seqs
│   │   ├── taxonomy_silva/            # SILVA classification
│   │   ├── taxonomy_gtdb/             # GTDB classification
│   │   ├── phylogeny/                 # Phylogenetic tree
│   │   └── diversity/                 # QIIME2 diversity outputs
│   ├── dada2_optimization/            # Parameter grid-search results
│   ├── picrust2/                      # PICRUSt2 output tables
│   │   ├── picrust2_out/              # Full PICRUSt2 output
│   │   ├── pathways_with_descriptions.tsv
│   │   ├── KO_metagenome_with_descriptions.tsv
│   │   └── EC_metagenome_with_descriptions.tsv
│   └── tables/                        # Exported data tables (CSV/TSV)
│       ├── asv_counts.tsv             # Raw ASV count table
│       ├── taxonomy_silva.tsv         # SILVA taxonomy assignments
│       ├── dada2_stats.tsv            # DADA2 processing statistics
│       ├── alpha_diversity_values.csv
│       ├── alpha_diversity_kruskal_wallis.csv
│       ├── beta_diversity_permanova.csv
│       ├── beta_diversity_pairwise_permanova.csv
│       ├── ancombc2_significant.csv
│       ├── composition_phylum.csv / composition_genus.csv
│       ├── core_microbiome.csv
│       ├── pathway_significant.csv / ko_significant.csv / ec_significant.csv
│       ├── pval_filtered/             # P-value filtered ANCOM-BC2 results
│       └── lefse/                     # LEfSe marker tables
├── figures/
│   ├── alpha_diversity/               # Alpha diversity plots
│   ├── beta_diversity/                # PCoA, NMDS ordinations
│   ├── composition/                   # Barplots, heatmaps, Venn diagrams
│   ├── differential_abundance/        # Volcano plots, ANCOM-BC2
│   ├── functional/                    # PICRUSt2 pathway/KO/EC figures
│   ├── lefse/                         # LEfSe LDA bar plots & cladograms
│   ├── qc/                            # Quality control plots
│   └── publication/                   # Final publication-ready figures
├── metadata/
│   ├── sample_metadata.tsv            # Sample metadata (group, disease, COPD, sample type)
│   └── manifest.tsv                   # QIIME2 import manifest
├── databases/
│   ├── silva/                         # SILVA 138-99 classifier
│   └── gtdb/                          # GTDB r214 classifier
├── envs/                              # Conda environment YAML specs
│   ├── qiime2-env.yml
│   ├── r-microbiome-env.yml
│   └── picrust2-env.yml
├── logs/                              # Pipeline execution logs
├── renv.lock                          # Exact R package versions (tracked)
├── .Rprofile                          # Auto-activates renv on R startup
├── renv/                              # (auto-generated, not in repo)
│   ├── activate.R                     # renv bootstrap script (tracked)
│   ├── settings.json                  # renv config (tracked)
│   └── library/                       # Installed packages (auto-generated, not in repo)
├── LC_COPD_microbiome.Rproj           # RStudio project file
└── README.md                          # This file
```

---

## 🔄 Pipeline Workflow

```mermaid
flowchart TD
    A["📁 Raw FASTQ Files<br/>(60 files, paired-end)"] --> B["01. Pre-trim QC<br/>(FastQC + MultiQC)"]
    B --> C["02. Primer Trimming<br/>(Cutadapt: 341F / 805R+806R)"]
    C --> D["03. Post-trim QC<br/>(FastQC + MultiQC)"]
    D --> E["04. QIIME2 Import"]
    E --> F["05. DADA2 Denoising<br/>(ASV table + Rep seqs)"]
    F --> G["06. SILVA Classification"]
    F --> H["07. GTDB Classification<br/>(optional)"]
    F --> I["08. Phylogeny<br/>(MAFFT + FastTree)"]
    G --> J["09. Diversity Analysis<br/>(Alpha + Beta)"]
    I --> J
    F --> K["10. PICRUSt2<br/>(Functional Prediction)"]

    J --> L["R Analysis Pipeline"]
    G --> L
    K --> L

    L --> L1["01. Import → phyloseq"]
    L1 --> L2["02. Alpha Diversity"]
    L1 --> L3["03. Beta Diversity"]
    L1 --> L4["04. Composition"]
    L1 --> L5["05. ANCOM-BC2"]
    L1 --> L8["08. LEfSe"]
    K --> L6["06. Pathway Visualization"]
    K --> L6b["06b. KO Visualization"]
    K --> L6c["06c. EC Visualization"]
    L2 & L3 & L4 & L5 & L6 & L8 --> L7["07. Publication Figures"]

    style A fill:#2b2d42,color:#edf2f4,stroke:#8d99ae
    style F fill:#2b2d42,color:#edf2f4,stroke:#8d99ae
    style L7 fill:#ef233c,color:#edf2f4,stroke:#d90429
```

---

## 🚀 Step-by-Step Execution Guide

> **Note:** All bash scripts should be run from the project root directory. Each script sources `config.sh` for shared parameters.

### Phase 1: Quality Control & Trimming

```bash
cd /path/to/LC_COPD_microbiome

# Step 0: One-time setup (environments, databases, directories)
bash scripts/00_setup_environment.sh

# Step 1: Pre-trimming quality assessment
bash scripts/01_fastqc_pretrim.sh
# Output: results/qc_reports/pretrim_fastqc/ + MultiQC report

# Step 2: Primer trimming with Cutadapt
# Removes 341F and mixed 805R/806R primers
# Both primers are searched for and trimmed simultaneously
bash scripts/02_trim_primers.sh
# Output: results/trimmed_reads/

# Step 3: Post-trimming quality assessment
bash scripts/03_fastqc_posttrim.sh
# Output: results/qc_reports/posttrim_fastqc/ + MultiQC report
```

### Phase 2: QIIME2 Processing

```bash
# Step 4: Import trimmed reads as QIIME2 artifact
bash scripts/04_qiime2_import.sh
# Output: results/qiime2/imported/

# Step 5: DADA2 denoising (~2-4 hours)
# Uses optimized parameters: truncLen F=275, R=220; maxEE F=2, R=4
bash scripts/05_dada2_denoise.sh
# Output: results/qiime2/dada2_output/, results/tables/dada2_stats.tsv

# Step 6: Taxonomy assignment with SILVA 138
bash scripts/06_taxonomy_silva.sh
# Output: results/qiime2/taxonomy_silva/, results/tables/taxonomy_silva.tsv

# Step 7 (Optional): Taxonomy assignment with GTDB r214
bash scripts/07_taxonomy_gtdb.sh
# Output: results/qiime2/taxonomy_gtdb/

# Step 8: Phylogenetic tree construction
bash scripts/08_phylogeny.sh
# Output: results/qiime2/phylogeny/, results/tables/phylogenetic_tree.nwk
```

### Phase 3: Diversity & Functional Analysis

```bash
# Step 9: Alpha and beta diversity (QIIME2)
bash scripts/09_diversity.sh
# Output: results/qiime2/diversity/, results/tables/

# Step 10: PICRUSt2 functional prediction (~1-2 hours)
bash scripts/10_picrust2.sh
# Output: results/picrust2/
```

### Phase 4: R Statistical Analysis

```bash
conda activate microbiome-r
cd R_analysis

# 01. Import QIIME2 artifacts into R phyloseq object
Rscript 01_import_phyloseq.R

# 02. Alpha diversity (Kruskal-Wallis, Wilcoxon, rarefaction curves)
Rscript 02_alpha_diversity.R

# 03. Beta diversity (PERMANOVA, PERMDISP, PCoA, NMDS)
Rscript 03_beta_diversity.R

# 04. Taxonomic composition (barplots, heatmaps, core microbiome, Venn diagrams)
Rscript 04_composition.R

# 05. Differential abundance (ANCOM-BC2 for each group comparison)
Rscript 05_ancombc2.R

# 06a. Functional pathway visualization (MetaCyc pathways)
Rscript 06_picrust_visualization.R

# 06b. KEGG Ortholog (KO) visualization
Rscript 06b_ko_visualization.R

# 06c. Enzyme Commission (EC) number visualization
Rscript 06c_ec_visualization.R

# 07. Generate all publication-quality composite figures
Rscript 07_publication_plots.R

# 08. LEfSe biomarker analysis (requires microbiomeMarker package)
Rscript 08_lefse.R
```

### One-Command Pipeline (Expert Mode)

```bash
# Run all bash scripts sequentially with logging
cd /path/to/LC_COPD_microbiome
for script in scripts/{01..10}*.sh; do
    echo "Running: $script"
    bash "$script" 2>&1 | tee "logs/$(basename $script .sh).log"
done

# Run all R scripts with logging
conda activate microbiome-r
for script in R_analysis/{01..08}*.R; do
    echo "Running: $script"
    Rscript "$script" 2>&1 | tee "logs/$(basename $script .R).log"
done
```

---

## 📊 R Analysis Scripts

### `01_import_phyloseq.R` — Data Import & Filtering

Imports QIIME2 artifacts (ASV table, taxonomy, phylogenetic tree, metadata) into a **phyloseq** object. Applies filtering to remove low-abundance/prevalence ASVs. Saves both raw and filtered phyloseq objects as `.rds` files.

### `02_alpha_diversity.R` — Alpha Diversity

Computes 8 alpha diversity indices:
- **Observed**, **Chao1**, **ACE** (richness)
- **Shannon**, **Simpson**, **Inverse Simpson**, **Fisher** (evenness/diversity)
- **Pielou's Evenness**

Statistical testing via Kruskal-Wallis (omnibus) and pairwise Wilcoxon rank-sum tests. Generates rarefaction curves to evaluate sequencing depth sufficiency.

### `03_beta_diversity.R` — Beta Diversity

Calculates four distance matrices:
- **Bray-Curtis** (abundance-weighted)
- **Jaccard** (presence/absence)
- **Weighted UniFrac** (phylogenetic + abundance)
- **Unweighted UniFrac** (phylogenetic only)

Statistical testing: PERMANOVA (adonis2, 999 permutations), pairwise PERMANOVA, and PERMDISP (betadisper) for homogeneity of dispersion. Generates PCoA and NMDS ordination plots.

### `04_composition.R` — Taxonomic Composition

Generates phylum- and genus-level stacked barplots (per sample and per group means), a top-30 genus heatmap (ComplexHeatmap), prevalence-abundance scatter, core microbiome analysis (prevalence ≥ 70%, detection threshold ≥ 0.01%), and a Venn diagram of shared/unique core taxa.

### `05_ancombc2.R` — Differential Abundance

Applies ANCOM-BC2 at the genus level with bias correction. Tests LC_COPD vs Control and LC_woCOPD vs Control. Produces volcano plots of log-fold changes with q-value-based significance. Saves both q-value and p-value filtered results.

### `06_picrust_visualization.R` / `06b_ko_visualization.R` / `06c_ec_visualization.R` — Functional Predictions

Visualizes PICRUSt2 predicted functional profiles:
- **MetaCyc pathways** (`06_picrust_visualization.R`)
- **KEGG Orthologs** (`06b_ko_visualization.R`)
- **Enzyme Commission numbers** (`06c_ec_visualization.R`)

Each script generates PCA ordinations, top-30 heatmaps, differential analysis (Kruskal-Wallis), and significant feature barplots. PERMANOVA tests functional community-level differences between groups.

### `07_publication_plots.R` — Publication Figures

Assembles publication-ready composite figures from individual analyses. Outputs in PDF, PNG, and TIFF formats at journal-quality resolution (300 DPI).

### `08_lefse.R` — LEfSe Biomarker Analysis

Runs Linear discriminant analysis Effect Size (LEfSe) using the `microbiomeMarker` R package. Identifies genus-level biomarker taxa distinguishing the three clinical groups. Generates LDA score barplots and cladograms.

---

## 📈 Key Results Summary

### DADA2 Processing Statistics

| Metric | Value |
|--------|-------|
| Total input reads | 1,395,788 |
| Mean input reads/sample | 46,526 |
| Mean non-chimeric reads/sample | 18,504 |
| Mean retention rate | 39.7% |
| DADA2 parameters | truncLen F=275, R=220; maxEE F=2, R=4 |
| Total ASVs identified | 6,409 (after prevalence filtering) |

> **DADA2 parameter optimization**: A grid-search across 36 parameter combinations (truncLen_F: 260/270/280, truncLen_R: 210/220/230/240, maxEE: (2,2)/(2,3)/(3,2)/(3,3)) was performed. The selected parameters (F=275, R=220, maxEE F=2, R=4) balanced read retention (~40%) with sufficient overlap (≥20 bp) for the ~460 bp V3-V4 amplicon, while using slightly relaxed maxEE for reverse reads to account for their typically lower quality.

#### Per-Sample DADA2 Processing

| Sample | Input | Filtered | % Filtered | Merged | % Merged | Non-chimeric | % Retained |
|--------|-------|----------|-----------|--------|----------|-------------|-----------|
| C1 | 53,911 | 40,707 | 75.5% | 17,828 | 33.1% | 14,749 | 27.4% |
| C2 | 48,355 | 35,965 | 74.4% | 26,469 | 54.7% | 18,081 | 37.4% |
| C3 | 45,837 | 30,778 | 67.2% | 23,815 | 52.0% | 16,178 | 35.3% |
| C4 | 39,542 | 22,507 | 56.9% | 17,048 | 43.1% | 9,175 | 23.2% |
| C5 | 47,946 | 35,284 | 73.6% | 25,823 | 53.9% | 23,402 | 48.8% |
| C6 | 46,571 | 35,009 | 75.2% | 23,659 | 50.8% | 19,847 | 42.6% |
| C7 | 49,314 | 38,113 | 77.3% | 27,656 | 56.1% | 24,259 | 49.2% |
| C8 | 46,931 | 36,149 | 77.0% | 24,982 | 53.2% | 22,872 | 48.7% |
| C9 | 48,483 | 36,579 | 75.5% | 25,468 | 52.5% | 24,096 | 49.7% |
| C10 | 47,928 | 36,297 | 75.7% | 24,927 | 52.0% | 23,101 | 48.2% |
| LC1 | 48,161 | 36,333 | 75.4% | 25,107 | 52.1% | 18,904 | 39.3% |
| LC2 | 46,594 | 36,170 | 77.6% | 27,500 | 59.0% | 21,379 | 45.9% |
| LC3 | 46,248 | 35,333 | 76.4% | 24,070 | 52.1% | 14,068 | 30.4% |
| LC4 | 47,680 | 35,305 | 74.1% | 25,542 | 53.6% | 17,969 | 37.7% |
| LC5 | 41,669 | 22,693 | 54.5% | 15,863 | 38.1% | 14,445 | 34.7% |
| LC6 | 46,164 | 33,120 | 71.7% | 24,983 | 54.1% | 18,465 | 40.0% |
| LC7 | 47,027 | 35,715 | 76.0% | 24,525 | 52.2% | 19,211 | 40.9% |
| LC8 | 47,110 | 35,769 | 75.9% | 25,785 | 54.7% | 20,617 | 43.8% |
| LC9 | 48,256 | 37,377 | 77.5% | 27,785 | 57.6% | 23,744 | 49.2% |
| LC10 | 40,644 | 23,360 | 57.5% | 17,065 | 42.0% | 12,853 | 31.6% |
| LC11 | 44,914 | 29,912 | 66.6% | 20,921 | 46.6% | 15,964 | 35.5% |
| LC12 | 46,633 | 35,636 | 76.4% | 25,608 | 54.9% | 16,539 | 35.5% |
| LC13 | 48,419 | 37,299 | 77.0% | 26,695 | 55.1% | 18,874 | 39.0% |
| LC14 | 48,349 | 37,067 | 76.7% | 28,527 | 59.0% | 21,193 | 43.8% |
| LC15 | 47,402 | 35,972 | 75.9% | 26,027 | 54.9% | 20,739 | 43.8% |
| LC16 | 45,767 | 31,246 | 68.3% | 21,746 | 47.5% | 16,789 | 36.7% |
| LC17 | 41,402 | 23,943 | 57.8% | 17,322 | 41.8% | 15,742 | 38.0% |
| LC18 | 46,278 | 34,053 | 73.6% | 23,571 | 50.9% | 19,778 | 42.7% |
| LC19 | 46,254 | 34,604 | 74.8% | 24,913 | 53.9% | 17,449 | 37.7% |
| LC20 | 45,999 | 34,574 | 75.2% | 23,977 | 52.1% | 18,326 | 39.8% |

---

### 16S rRNA Sequencing Summary

The 16S rRNA gene sequencing (V3-V4 region) was used to assess the changes in the BALF microbiota between the Control, LC_COPD, and LC_woCOPD groups. A total of 1,395,788 high-quality raw reads were generated using the Illumina MiSeq platform. We obtained 558,808 clean reads (19,576.0, 18,165.5, and 18,139.3 average reads for Control, LC_COPD, and LC_woCOPD groups, respectively) after filtering and chimera removal using DADA2. A total of 6,409 ASVs were identified after prevalence filtering (>1% of samples). These 6,409 ASVs belonged to 584 different genera in 43 different phyla. Of these ASVs, 227 were shared among the three groups, and 1,796, 1,929, and 1,870 ASVs were specific to the Control, LC_COPD, and LC_woCOPD groups, respectively.

#### Taxonomic Distribution

| Rank | Unique Taxa |
|------|------------|
| Kingdom | 2 |
| Phylum | 43 |
| Class | 107 |
| Order | 227 |
| Family | 352 |
| Genus | 584 |
| Species | 508 |

#### ASV Distribution Across Groups (Venn Diagram)

| Venn Region | ASV Count |
|-------------|-----------|
| Shared (all 3 groups) | 227 |
| Control & LC_COPD only | 298 |
| Control & LC_woCOPD only | 59 |
| LC_COPD & LC_woCOPD only | 230 |
| Unique to Control | 1,796 |
| Unique to LC_COPD | 1,929 |
| Unique to LC_woCOPD | 1,870 |
| **Total** | **6,409** |

> **Figure**: `figures/composition/asv_venn_diagram.png`

### Sample Summary

| Group | N | Mean Reads | SD | Min Reads | Max Reads | Mean ASVs | SD |
|-------|---|-----------|-----|-----------|-----------|-----------|-----|
| Control | 10 | 18,872 | 4,753 | 9,160 | 24,150 | 321 | 167 |
| LC_COPD | 10 | 17,466 | 3,647 | 12,617 | 23,536 | 364 | 143 |
| LC_woCOPD | 10 | 17,175 | 3,123 | 10,217 | 20,891 | 360 | 143 |

---

### Alpha Diversity

All alpha diversity metrics (8 total) showed **no statistically significant differences** (BH-adjusted p > 0.05) across the three groups by Kruskal-Wallis test:

| Metric | Category | H statistic | df | p-value | Adjusted p | Significance |
|--------|----------|-------------|-----|---------|------------|--------------|
| Observed | Richness | 0.903 | 2 | 0.637 | 0.651 | ns |
| Chao1 | Richness | 0.859 | 2 | 0.651 | 0.651 | ns |
| ACE | Richness | 0.947 | 2 | 0.623 | 0.651 | ns |
| Shannon | Diversity | 3.440 | 2 | 0.179 | 0.579 | ns |
| Simpson | Diversity | 2.480 | 2 | 0.289 | 0.579 | ns |
| InvSimpson | Diversity | 2.480 | 2 | 0.289 | 0.579 | ns |
| Fisher | Diversity | 1.025 | 2 | 0.599 | 0.651 | ns |
| Pielou | Evenness | 2.766 | 2 | 0.251 | 0.579 | ns |

> **Interpretation**: Similar species richness and evenness across groups suggests that overall community diversity is not globally altered in lung cancer patients. However, compositional differences may exist at the individual taxon or functional level (see below).
>
> **Figures**: `figures/alpha_diversity/alpha_diversity_main.png`, `figures/alpha_diversity/rarefaction_curves.png`

---

### Beta Diversity

#### PERMANOVA (Overall – `adonis2`, 999 permutations)

Community composition differed significantly between groups for multiple distance metrics:

| Distance Metric | Df | Sum of Squares | R² | F | p-value | BH-adj. p | Significant? |
|-----------------|-----|---------------|-----|---|---------|-----------|-------------|
| **Bray-Curtis** | 2 | 0.742 | 0.132 | 2.050 | 0.021 | **0.042** | ✅ Yes |
| **Jaccard** | 2 | 1.128 | 0.089 | 1.319 | 0.001 | **0.004** | ✅ Yes |
| Weighted UniFrac | 2 | 0.125 | 0.141 | 2.214 | 0.058 | 0.058 | ❌ No (marginal) |
| **Unweighted UniFrac** | 2 | 0.684 | 0.083 | 1.225 | 0.032 | **0.043** | ✅ Yes |

> **Interpretation**: Community composition differs significantly by Bray-Curtis (abundance-weighted), Jaccard (presence/absence), and Unweighted UniFrac (phylogeny-aware presence/absence). Weighted UniFrac is borderline (p = 0.058), suggesting the dominant phylogenetic lineages are relatively conserved across groups. Effect sizes (R² = 8–14%) are modest but typical for human microbiome studies with small sample sizes.

#### PERMDISP (Beta Dispersion – `betadisper` + `permutest`, 999 permutations)

| Distance Metric | F | p-value | Significant? |
|-----------------|---|---------|-------------|
| Bray-Curtis | 0.172 | 0.857 | ❌ No |
| Jaccard | 1.044 | 0.366 | ❌ No |
| Weighted UniFrac | 1.840 | 0.187 | ❌ No |
| Unweighted UniFrac | 0.089 | 0.922 | ❌ No |

> **Interpretation**: All PERMDISP tests are non-significant (all p > 0.18), confirming that within-group dispersions are homogeneous across all three groups. This validates that the significant PERMANOVA results reflect genuine differences in community centroids (compositional location), not artifacts of unequal group heterogeneity.

#### Pairwise PERMANOVA (All Distance Metrics)

| Comparison | Distance | R² | F | Raw p | BH-adj. p | Sig? |
|------------|----------|-----|---|-------|-----------|------|
| Control vs LC_COPD | Bray-Curtis | 0.081 | 1.580 | 0.122 | 0.122 | ❌ |
| **Control vs LC_woCOPD** | **Bray-Curtis** | **0.126** | **2.588** | **0.021** | **0.045** | ✅ |
| **LC_COPD vs LC_woCOPD** | **Bray-Curtis** | **0.100** | **2.008** | **0.030** | **0.045** | ✅ |
| Control vs LC_COPD | Jaccard | 0.058 | 1.112 | 0.135 | 0.135 | ❌ |
| **Control vs LC_woCOPD** | **Jaccard** | **0.077** | **1.507** | **0.002** | **0.005** | ✅ |
| **LC_COPD vs LC_woCOPD** | **Jaccard** | **0.069** | **1.344** | **0.003** | **0.005** | ✅ |
| Control vs LC_COPD | Weighted UniFrac | 0.100 | 1.999 | 0.097 | 0.146 | ❌ |
| Control vs LC_woCOPD | Weighted UniFrac | 0.133 | 2.750 | 0.054 | 0.146 | ❌ |
| LC_COPD vs LC_woCOPD | Weighted UniFrac | 0.089 | 1.761 | 0.151 | 0.151 | ❌ |
| Control vs LC_COPD | Unweighted UniFrac | 0.060 | 1.154 | 0.166 | 0.249 | ❌ |
| **Control vs LC_woCOPD** | **Unweighted UniFrac** | **0.074** | **1.440** | **0.015** | **0.045** | ✅ |
| LC_COPD vs LC_woCOPD | Unweighted UniFrac | 0.057 | 1.083 | 0.252 | 0.252 | ❌ |

> **Key finding**: The **LC_woCOPD (lung cancer without COPD)** group is consistently the most compositionally distinct, separating significantly from both Control and LC_COPD by Bray-Curtis and Jaccard, and from Control by Unweighted UniFrac. The Control vs LC_COPD comparison is never significant after correction, meaning the COPD-comorbid lung cancer group resembles healthy controls more than the LC-only group. Weighted UniFrac shows no significant pairwise differences, reinforcing that dominant phylogenetic lineages are conserved.

#### NMDS Ordination

- **Method**: Non-Metric Multidimensional Scaling on Bray-Curtis distances
- **Stress value**: 0.154 (adequate, < 0.2)
- **Convergence**: Best solution found and repeated after 39 iterations (trymax = 100)

> **Figures**: `figures/beta_diversity/beta_diversity_pcoa_combined.png`, `figures/beta_diversity/nmds_bray_curtis.png`, `figures/beta_diversity/distance_comparison_boxplot.pdf`

---

### Differential Abundance (ANCOM-BC2)

ANCOM-BC2 was applied at the genus level (55 genera tested after prevalence/abundance filtering from 177 genera input) with bias correction. Reference group: Control.

#### Significant Taxa (q-value < 0.05)

| Genus | Phylum | Family | LC_COPD vs Control (LFC ± SE) | q-value | LC_woCOPD vs Control (LFC ± SE) | q-value |
|-------|--------|--------|-------------------------------|---------|--------------------------------|---------|
| ***Neisseria*** | Proteobacteria | Neisseriaceae | **↑ 2.76 ± 0.59** | **0.007** | ↑ 1.28 ± 0.65 | 0.236 |
| ***Cetobacterium*** | Fusobacteriota | Fusobacteriaceae | ↑ 1.71 ± 0.46 | 0.055 | **↑ 1.92 ± 0.43** | **0.027** |
| ***Gemella*** | Firmicutes | Gemellaceae | ↑ 1.08 ± 0.44 | 0.235 | **↑ 1.82 ± 0.48** | **0.035** |

#### Additional Differentially Abundant Taxa (p-value < 0.05)

With relaxed filtering using raw p-values, **7 taxa** were differentially abundant in LC_COPD vs Control and **9 taxa** in LC_woCOPD vs Control:

| Genus | Phylum | LC_COPD vs Control | p-value | LC_woCOPD vs Control | p-value |
|-------|--------|-------------------|---------|---------------------|---------|
| *Neisseria* | Proteobacteria | ↑ LFC = 2.37 | **0.002** | — | ns |
| *Cetobacterium* | Fusobacteriota | ↑ LFC = 1.71 | **0.002** | ↑ LFC = 1.92 | **< 0.001** |
| *Unclassified Prevotellaceae* | Bacteroidota | ↓ LFC = −1.16 | **0.004** | — | ns |
| *Oribacterium* | Firmicutes | ↓ LFC = −0.97 | **0.020** | — | ns |
| *Gemella* | Firmicutes | ↑ LFC = 1.06 | **0.023** | ↑ LFC = 1.82 | **< 0.001** |
| *Streptococcus* | Firmicutes | ↑ LFC = 1.61 | **0.027** | ↑ LFC = 1.70 | 0.014 |
| *Unclassified Comamonadaceae* | Proteobacteria | ↓ LFC = −0.87 | **0.027** | — | ns |
| *Corynebacterium* | Actinobacteriota | — | ns | ↑ LFC = 1.23 | **0.036** |
| *Leptotrichia* | Fusobacteriota | — | ns | ↑ LFC = 1.47 | **0.009** |
| *Prevotellaceae UCG-001* | Bacteroidota | — | ns | ↑ LFC = 0.99 | **0.018** |
| *Rikenellaceae RC9 gut group* | Bacteroidota | — | ns | ↓ LFC = −0.96 | **0.021** |
| *Unclassified Micrococcaceae* | Actinobacteriota | — | ns | ↓ LFC = −1.49 | **0.006** |
| *Unclassified Enterobacteriaceae* | Proteobacteria | — | ns | ↓ LFC = −0.87 | **0.060** |

> **Figures**: `figures/differential_abundance/pval_filtered/volcano_plots_combined.png`, `figures/differential_abundance/pval_filtered/lfc_barplots_combined.png`, `figures/differential_abundance/pval_filtered/heatmap_significant_taxa.pdf`, `figures/differential_abundance/qval_filtered/volcano_plots_combined.png`, `figures/differential_abundance/qval_filtered/lfc_barplots_combined.png`

---

### LEfSe Biomarker Analysis

Linear discriminant analysis Effect Size (LEfSe) identified genus-level biomarkers distinguishing the three clinical groups (LDA score > 2.0, p < 0.05). A total of **71 significant biomarkers** were identified across three pairwise comparisons.

#### LC_COPD vs Control (25 biomarkers)

**Higher in LC_COPD (16 taxa):**

| Genus | Phylum | LDA Score | p-value |
|-------|--------|-----------|---------|
| *Streptococcus* | Firmicutes | 4.97 | 0.019 |
| *Neisseria* | Proteobacteria | 4.64 | 0.010 |
| *Cetobacterium* | Fusobacteriota | 4.06 | 0.001 |
| *Gemella* | Firmicutes | 3.62 | 0.007 |
| *Capnocytophaga* | Bacteroidota | 3.31 | 0.005 |
| *Peptostreptococcus* | Firmicutes | 3.27 | 0.022 |
| *Aggregatibacter* | Proteobacteria | 3.22 | 0.031 |
| *Anaerobacillus* | Firmicutes | 3.08 | 0.013 |
| *Bacteroides* | Bacteroidota | 3.01 | 0.017 |
| *ZOR0006* | Firmicutes | 2.84 | 0.008 |
| *Desulfobulbus* | Desulfobacterota | 2.61 | 0.017 |
| *Dongia* | Proteobacteria | 2.04 | 0.041 |
| *Oleiagrimonas* | Proteobacteria | 2.03 | 0.031 |

**Higher in Control (9 taxa):**

| Genus | Phylum | LDA Score | p-value |
|-------|--------|-----------|---------|
| *Turicibacter* | Firmicutes | 2.58 | 0.026 |
| *Brevundimonas* | Proteobacteria | 2.54 | 0.024 |
| *Candidatus Saccharimonas* | Patescibacteria | 2.47 | 0.031 |
| *[Eubacterium] coprostanoligenes group* | Firmicutes | 2.14 | 0.013 |

#### LC_woCOPD vs Control (26 biomarkers)

**Higher in LC_woCOPD/LC (13 taxa):**

| Genus | Phylum | LDA Score | p-value |
|-------|--------|-----------|---------|
| *Alloprevotella* | Bacteroidota | 4.80 | 0.034 |
| *Streptococcus* | Firmicutes | 4.55 | 0.041 |
| *Corynebacterium* | Actinobacteriota | 4.01 | 0.016 |
| *Cetobacterium* | Fusobacteriota | 3.95 | 0.026 |
| *Gemella* | Firmicutes | 3.82 | < 0.001 |
| *Granulicatella* | Firmicutes | 3.41 | 0.004 |
| *Dolosigranulum* | Firmicutes | 3.29 | 0.005 |
| *Aggregatibacter* | Proteobacteria | 3.28 | 0.013 |
| *Streptobacillus* | Fusobacteriota | 3.19 | 0.005 |

**Higher in Control (13 taxa):**

| Genus | Phylum | LDA Score | p-value |
|-------|--------|-----------|---------|
| *Rikenellaceae RC9 gut group* | Bacteroidota | 2.88 | 0.006 |
| *Proteus* | Proteobacteria | 2.60 | 0.047 |
| *Solobacterium* | Firmicutes | 2.44 | 0.020 |
| *Desulfovibrio* | Desulfobacterota | 2.14 | 0.025 |
| *Mycobacterium* | Actinobacteriota | 2.11 | 0.036 |

#### LC_COPD vs LC_woCOPD (20 biomarkers)

**Higher in LC_COPD (11 taxa):**

| Genus | Phylum | LDA Score | p-value |
|-------|--------|-----------|---------|
| *Clostridium sensu stricto 1* | Firmicutes | 4.13 | 0.031 |
| *Capnocytophaga* | Bacteroidota | 3.23 | 0.018 |
| *Rummeliibacillus* | Firmicutes | 3.21 | 0.013 |
| *Anaerobacillus* | Firmicutes | 3.03 | 0.036 |
| *Desulfovibrio* | Desulfobacterota | 2.54 | 0.025 |
| *Abiotrophia* | Firmicutes | 2.46 | 0.031 |

**Higher in LC_woCOPD (9 taxa):**

| Genus | Phylum | LDA Score | p-value |
|-------|--------|-----------|---------|
| *Corynebacterium* | Actinobacteriota | 4.04 | 0.049 |
| *Actinobacillus* | Proteobacteria | 3.64 | 0.031 |
| *Granulicatella* | Firmicutes | 3.39 | 0.026 |
| *Dolosigranulum* | Firmicutes | 3.32 | 0.013 |
| *Streptobacillus* | Fusobacteriota | 3.23 | 0.005 |

> **Figure**: `figures/lefse/combined/lefse_combined_lda_barplots.png`

---

### Core Microbiome

#### Core Taxa at ≥50% Prevalence (>0.1% abundance)

A total of **31 core ASVs** were identified across all samples at 50% prevalence threshold. The number of core taxa per group:

| Group | Core Taxa Count |
|-------|----------------|
| Control | 21 |
| LC_COPD | 31 |
| LC_woCOPD | 37 |

Key core genera include: *Prevotella* (3 ASVs), *Veillonella* (3 ASVs), *Streptococcus* (3 ASVs), *Bacillus* (2 ASVs), *Dickeya* (6 ASVs), *Collinsella* (2 ASVs), *Alloprevotella*, *Actinomyces*, *Cutibacterium*, and *Aerococcus*.

#### Core Taxa at ≥80% Prevalence (>0.1% abundance)

Only **6 highly conserved core ASVs** were present at 80% prevalence:

| Core Taxon | Phylum | Family | Genus |
|-----------|--------|--------|-------|
| ASV 1 | Firmicutes | Veillonellaceae | *Veillonella* |
| ASV 2 | Firmicutes | Streptococcaceae | *Streptococcus* |
| ASV 3 | Firmicutes | Bacillaceae | *Bacillus* |
| ASV 4 | Firmicutes | Bacillaceae | *Bacillus* |
| ASV 5 | Proteobacteria | Enterobacterales (unclassified) | — |
| ASV 6 | Proteobacteria | Enterobacterales (unclassified) | — |

The core taxa per group at 80%: Control = 3, LC_COPD = 7, LC_woCOPD = 11.

> **Figures**: `figures/composition/core_microbiome_venn_50pct.png`, `figures/composition/core_microbiome_venn_80pct.png`

---

### Functional Predictions (PICRUSt2)

#### PERMANOVA on Functional Profiles

PERMANOVA was performed on Bray-Curtis distances of predicted functional profiles. **No significant overall community-level differences** were observed for any functional category:

| Functional Category | Df | R² | F | p-value | Significant? |
|--------------------|-----|-----|---|---------|-------------|
| MetaCyc Pathways (422 pathways) | 2 | 0.102 | 1.536 | 0.174 | ❌ No |
| KEGG Orthologs (7,487 KOs) | 2 | 0.085 | 1.261 | 0.269 | ❌ No |
| Enzyme Commission (2,253 ECs) | 2 | 0.092 | 1.365 | 0.229 | ❌ No |

#### Differentially Abundant Features (Kruskal-Wallis, p < 0.05)

Individual feature-level testing identified nominally significant differences, though **none survived FDR correction**:

| Category | Features Tested | Significant (p < 0.05) | Significant (q < 0.05) |
|----------|----------------|----------------------|----------------------|
| MetaCyc Pathways | 422 | 48 | 0 |
| KEGG Orthologs | 7,487 | 653 | 0 |
| Enzyme Commission | 2,253 | 178 | 0 |

#### Top Significant MetaCyc Pathways (p < 0.05)

| Pathway | Description | Control (mean) | LC_COPD (mean) | LC_woCOPD (mean) | p-value | Trend |
|---------|-------------|---------------|----------------|-----------------|---------|-------|
| GLYCOCAT-PWY | Glycogen degradation I (bacterial) | 0.336 | 0.399 | 0.441 | 0.003 | ↑ in LC groups |
| P441-PWY | N-acetylneuraminate degradation | 0.273 | 0.388 | 0.357 | 0.004 | ↑ in LC_COPD |
| PPGPPMET-PWY | ppGpp biosynthesis | 0.108 | 0.172 | 0.128 | 0.006 | ↑ in LC_COPD |
| PWY-5913 | TCA cycle VI (obligate autotrophs) | 0.238 | 0.331 | 0.328 | 0.007 | ↑ in LC groups |
| P101-PWY | Ectoine biosynthesis | 0.004 | 0.010 | 0.004 | 0.010 | ↑ in LC_COPD |
| PWY-5910 | Geranylgeranyldiphosphate biosynthesis I | 0.044 | 0.121 | 0.090 | 0.011 | ↑ in LC groups |
| PWY-922 | Mevalonate pathway I | 0.032 | 0.094 | 0.067 | 0.011 | ↑ in LC groups |
| P164-PWY | Purine nucleobases degradation I (anaerobic) | 0.144 | 0.210 | 0.245 | 0.011 | ↑ in LC groups |
| FASYN-ELONG-PWY | Fatty acid elongation (saturated) | 0.627 | 0.599 | 0.575 | 0.012 | ↓ in LC groups |
| PWY-6471 | Peptidoglycan biosynthesis IV (*E. faecium*) | 0.127 | 0.286 | 0.217 | 0.012 | ↑ in LC_COPD |
| PWY-7228 | Guanosine nucleotides de novo biosynthesis I | 0.507 | 0.532 | 0.538 | 0.014 | ↑ in LC groups |
| P23-PWY | Reductive TCA cycle I | 0.293 | 0.369 | 0.360 | 0.015 | ↑ in LC groups |
| LACTOSECAT-PWY | Lactose and galactose degradation I | 0.037 | 0.102 | 0.066 | 0.024 | ↑ in LC_COPD |
| GLCMANNANAUT-PWY | N-acetylglucosamine/-mannosamine/-neuraminate degradation | 0.130 | 0.212 | 0.182 | 0.025 | ↑ in LC groups |
| PWY0-1061 | L-alanine biosynthesis | 0.291 | 0.395 | 0.329 | 0.040 | ↑ in LC_COPD |

> **Figures**: `figures/functional/pathway_pca.png`, `figures/functional/significant_pathways_barplot.png`, `figures/functional/pathway_heatmap_top30.pdf`, `figures/functional/ko_pca.png`, `figures/functional/significant_ko_barplot.png`, `figures/functional/ec_pca.png`, `figures/functional/significant_ec_barplot.png`, `figures/functional/log2fc_heatmap_all3_top10_significance.png`

---

### Summary of All Generated Figures

#### Publication-Quality Figures (`figures/publication/`)

| Figure | Description | Formats |
|--------|-------------|---------|
| `Figure1_alpha_diversity` | Alpha diversity boxplots (Observed, Shannon, Simpson, Pielou) with Kruskal-Wallis statistics | PDF, PNG, TIFF |
| `Figure2_beta_diversity` | Combined PCoA ordination (4 metrics) with PERMANOVA annotations and 95% confidence ellipses | PDF, PNG, TIFF |
| `Figure3_composition` | Phylum-level stacked barplot by group | PDF, PNG, TIFF |
| `Figure4_differential_abundance` | ANCOM-BC2 log-fold change barplots for significant taxa | PDF, PNG, TIFF |
| `Figure5_functional_pathways` | PICRUSt2 pathway PCA + top significant pathway barplots | PDF, PNG, TIFF |
| `FigureS1_rarefaction` | Rarefaction curves for all 30 samples | PDF, PNG |
| `FigureS2_volcano_plots` | Volcano plots for all ANCOM-BC2 pairwise comparisons | PDF, PNG, TIFF |
| `FigureS3a_genus_heatmap` | Top genus-level abundance heatmap (ComplexHeatmap) | PDF, PNG |
| `FigureS3b_core_venn` | Core microbiome Venn diagram (shared/unique taxa) | PNG |
| `FigureS4_ko_ec_barplots` | Significant KO & EC barplots from PICRUSt2 | PDF, PNG, TIFF |
| `Table1_sample_summary` | Sample group summary statistics | CSV |

#### Individual Analysis Figures

| Directory | Contents |
|-----------|----------|
| `figures/alpha_diversity/` | Main boxplots, supplementary metrics, violin plots, rarefaction curves |
| `figures/beta_diversity/` | Individual PCoA plots (4 metrics), combined PCoA panel, NMDS (Bray-Curtis), distance comparison boxplot |
| `figures/composition/` | Phylum barplots (per sample + per group), genus barplots, top-20 and top-30 genus heatmaps, prevalence-abundance scatter, core microbiome Venn diagrams (50% and 80%), ASV Venn diagram |
| `figures/differential_abundance/pval_filtered/` | Volcano plots (5 comparisons + combined), LFC barplots (5 + combined), waterfall plots (5), significance heatmap |
| `figures/differential_abundance/qval_filtered/` | Volcano plots (3 comparisons + combined), LFC barplots (3 + combined), waterfall plots (3), significance heatmap |
| `figures/functional/` | PCA ordinations (pathway, KO, EC), annotated PCAs, top-30 heatmaps, significant feature barplots, log2FC heatmaps (magnitude and significance ranking) |
| `figures/lefse/` | LDA score barplots for each pairwise comparison and combined panel |

---

## 📄 Output Files

### Key Data Tables

| File | Description |
|------|-------------|
| `results/tables/asv_counts.tsv` | Raw ASV count table (all samples × ASVs) |
| `results/tables/asv_counts_filtered_silva.tsv` | Filtered ASV table after taxonomy-based QC |
| `results/tables/taxonomy_silva.tsv` | SILVA taxonomic assignments |
| `results/tables/dada2_stats.tsv` | Per-sample DADA2 processing statistics |
| `results/tables/alpha_diversity_values.csv` | Per-sample alpha diversity indices |
| `results/tables/alpha_diversity_kruskal_wallis.csv` | Kruskal-Wallis test results |
| `results/tables/beta_diversity_permanova.csv` | Omnibus PERMANOVA results |
| `results/tables/beta_diversity_pairwise_permanova.csv` | Pairwise PERMANOVA results |
| `results/tables/beta_diversity_dispersion.csv` | PERMDISP results |
| `results/tables/ancombc2_results.csv` | Full ANCOM-BC2 output |
| `results/tables/ancombc2_significant.csv` | Significant taxa (q < 0.05) |
| `results/tables/composition_phylum.csv` | Phylum-level relative abundance |
| `results/tables/composition_genus.csv` | Genus-level relative abundance |
| `results/tables/core_microbiome.csv` | Core microbiome taxa |
| `results/tables/pathway_significant.csv` | Significant MetaCyc pathways |
| `results/tables/ko_significant.csv` | Significant KEGG Orthologs |
| `results/tables/ec_significant.csv` | Significant Enzyme Commission numbers |
| `results/tables/distance_matrix_*.csv` | Beta diversity distance matrices |
| `results/tables/phylogenetic_tree.nwk` | Newick-format phylogenetic tree |

### Publication Figures

| Figure | Description | Formats |
|--------|-------------|---------|
| `Figure1_alpha_diversity` | Alpha diversity boxplots across groups | PDF, PNG, TIFF |
| `Figure2_beta_diversity` | PCoA ordination (4 metrics combined) | PDF, PNG, TIFF |
| `Figure3_composition` | Taxonomic composition barplots | PDF, PNG, TIFF |
| `Figure4_differential_abundance` | ANCOM-BC2 differential abundance | PDF, PNG, TIFF |
| `Figure5_functional_pathways` | PICRUSt2 pathway analysis | PDF, PNG, TIFF |
| `FigureS1_rarefaction` | Rarefaction curves | PDF, PNG |
| `FigureS2_volcano_plots` | Volcano plots per comparison | PDF, PNG, TIFF |
| `FigureS3a_genus_heatmap` | Genus-level abundance heatmap | PDF, PNG |
| `FigureS3b_core_venn` | Core microbiome Venn diagram | PNG |
| `FigureS4_ko_ec_barplots` | KO & EC significant barplots | PDF, PNG, TIFF |
| `Table1_sample_summary` | Sample group summary statistics | CSV |

All figures are generated at **300 DPI** resolution for journal submission.

---

## ⚙ Parameters & Configuration

All pipeline parameters are centrally defined in `scripts/config.sh`. Key parameters:

### Computational Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `THREADS` | 6 | CPU threads (Ryzen 5 8645HS: 6 cores) |
| `MAX_MEMORY` | 10G | Memory limit (~11 GB available) |

### Primer Trimming

| Parameter | Value | Description |
|-----------|-------|-------------|
| `FORWARD_PRIMER` | `CCTACGGGNGGCWGCAG` | 341F forward primer (17 bp) |
| `REVERSE_PRIMER_805R` | `GACTACHVGGGTATCTAATCC` | 805R reverse primer (21 bp) |
| `REVERSE_PRIMER_806R` | `GGACTACHVGGGTWTCTAAT` | 806R reverse primer (20 bp) |

### DADA2 Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `TRUNC_LEN_F` | 275 | Optimized via grid-search; maintains quality & overlap |
| `TRUNC_LEN_R` | 220 | Ensures ≥20 bp overlap for ~460 bp V3-V4 amplicon |
| `MAX_EE_F` | 2 | Standard stringency for forward reads |
| `MAX_EE_R` | 4 | Slightly relaxed for lower-quality reverse reads |
| `TRUNC_Q` | 2 | Truncate reads at first base with quality ≤ 2 |
| `CHIMERA_METHOD` | consensus | Chimera removal strategy |
| `RANDOM_SEED` | 42 | For reproducibility |

### Diversity Parameters

| Parameter | Value |
|-----------|-------|
| `ALPHA_METRICS` | observed_features, shannon, simpson, faith_pd, pielou_e, chao1 |
| `BETA_METRICS` | bray_curtis, jaccard, weighted_unifrac, unweighted_unifrac |

### Databases

| Database | Version | Purpose |
|----------|---------|---------|
| SILVA | 138-99-515-806 | Primary taxonomic classifier |
| GTDB | r214 | Alternative / validation classifier |

---

## 🔁 Reproducing the Analysis

### Full Reproducibility Checklist

1. **Random seed**: Set to `42` throughout all scripts (bash and R)
2. **R packages**: Exact versions locked via `renv.lock` (R 4.5.2, Bioconductor 3.22) — see [Step 4](#step-4-restore-r-package-versions-renv) for restore instructions
3. **Conda environments**: Exported in `envs/*.yml`
4. **Parameters**: All in `scripts/config.sh` (no hard-coded values in individual scripts)
5. **Logs**: Execution logs in `logs/` directory with timestamps

### Restoring the R Environment from GitHub

When you clone this repository, the `renv/library/` folder will **not** be present (it is gitignored). To recreate the exact R environment:

```bash
conda activate microbiome-r
cd LC_COPD_microbiome

# Install renv (if needed)
Rscript -e 'install.packages("renv")'

# Restore all packages from the lockfile
# This reads renv.lock and installs every package at the exact recorded version
Rscript -e 'renv::restore()'
```

If you encounter issues (e.g., compilation errors on certain platforms), you can install packages manually:

```r
# In R, within the project directory:
install.packages("BiocManager")
BiocManager::install(version = "3.22")

# Install key Bioconductor packages
BiocManager::install(c(
  "phyloseq", "ANCOMBC", "microbiomeMarker",
  "ComplexHeatmap", "Biostrings"
), ask = FALSE)

# Install CRAN packages
install.packages(c(
  "tidyverse", "vegan", "ape", "VennDiagram",
  "ggrepel", "patchwork", "pheatmap", "RColorBrewer"
))
```

### Export Environment for Reproducibility

```bash
# Export QIIME2 environment
conda activate qiime2-amplicon-2024.10
conda env export > envs/qiime2-env-export.yml

# Export R environment
conda activate microbiome-r
conda env export > envs/r-microbiome-env-export.yml

# Export PICRUSt2 environment
conda activate picrust2
conda env export > envs/picrust2-env-export.yml
```

### Software Versions

Record all versions:

```bash
conda activate qiime2-amplicon-2024.10
qiime info
conda list > logs/qiime2_packages.txt

conda activate microbiome-r
R --version
Rscript -e 'sessionInfo()' > logs/r_session_info.txt

conda activate picrust2
picrust2_pipeline.py --version
```

---

## 🔧 Troubleshooting

### Common Issues

**1. QIIME2 environment not found**

```bash
# Install using official method
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py39-linux-conda.yml
conda env create -n qiime2-amplicon-2024.10 \
    --file qiime2-amplicon-2024.10-py39-linux-conda.yml
```

**2. Out of memory during DADA2 denoising**

```bash
# Reduce threads and memory in config.sh
export THREADS=4
export MAX_MEMORY="8G"
# Or run on a machine with more RAM
```

**3. R package installation fails**

```r
# Install Bioconductor packages directly
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC", ask = FALSE)
BiocManager::install("phyloseq", ask = FALSE)
BiocManager::install("microbiomeMarker", ask = FALSE)

# For ComplexHeatmap
BiocManager::install("ComplexHeatmap", ask = FALSE)
```

**4. GTDB classifier not available**

- The GTDB classifier requires training with RESCRIPt or downloading a pre-trained version.
- SILVA classification is sufficient for most analyses and is the primary classifier used here.
- GTDB step (`07_taxonomy_gtdb.sh`) is optional and can be skipped.

**5. PICRUSt2 taking too long**

```bash
# Use stratified mode off (default) for faster execution
# Ensure correct BIOM format:
biom convert -i results/tables/asv_counts.tsv \
    -o results/tables/feature-table.biom \
    --table-type="OTU table" --to-hdf5
```

**6. Cutadapt primer detection issues**

This dataset has mixed reverse primers (805R: ~39%, 806R: ~89%). The pipeline handles this by searching for both primers and trimming reads accordingly. If primer detection rates appear low, verify primer sequences match your experiment.

**7. Low DADA2 retention rates**

A retention rate of ~40% is acceptable for V3-V4 amplicon data, particularly with mixed primer sets and slightly lower reverse-read quality. If rates are consistently below 20%, consider:
- Adjusting `TRUNC_LEN_F` and `TRUNC_LEN_R` based on quality profiles
- Relaxing `MAX_EE` values
- Re-running `scripts/optimize_dada2_params.R` with a broader grid

---

## 📚 Citation

If you use this pipeline or its outputs, please cite the following tools:

| Tool | Reference |
|------|-----------|
| **QIIME2** | Bolyen et al. (2019) *Nature Biotechnology*, 37, 852–857 |
| **DADA2** | Callahan et al. (2016) *Nature Methods*, 13, 581–583 |
| **SILVA** | Quast et al. (2013) *Nucleic Acids Research*, 41, D590–D596 |
| **GTDB** | Parks et al. (2022) *Nucleic Acids Research*, 50, D199–D210 |
| **phyloseq** | McMurdie & Holmes (2013) *PLoS ONE*, 8, e61217 |
| **ANCOM-BC2** | Lin & Peddada (2020) *Nature Communications*, 11, 3514 |
| **PICRUSt2** | Douglas et al. (2020) *Nature Biotechnology*, 38, 685–688 |
| **LEfSe** | Segata et al. (2011) *Genome Biology*, 12, R60 |
| **microbiomeMarker** | Cao et al. (2022) *Bioinformatics*, 38, 4643–4651 |
| **vegan** | Oksanen et al. (2022) R package version 2.6 |
| **ComplexHeatmap** | Gu et al. (2016) *Bioinformatics*, 32, 2847–2849 |

---

## 📬 Contact

For questions, issues, or contributions, please [open an issue](../../issues) on this repository.

---

*Pipeline last updated: 2026-03-12*
*QIIME2 version: 2024.10 | R version: 4.3+ | DADA2 parameters optimized via grid-search*
