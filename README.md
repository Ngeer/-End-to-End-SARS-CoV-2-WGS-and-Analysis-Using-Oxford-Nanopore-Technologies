# End-to-End SARS-CoV-2 Whole Genome Sequencing and Analysis Using Oxford Nanopore Technologies

> **Oxford Nanopore Technologies (ONT) sequencing** | 20 SARS-CoV-2 samples | Pipeline: Docker + R | Data: NCBI SRA

---

## Project overview

This project demonstrates a complete, reproducible end-to-end bioinformatics pipeline for SARS-CoV-2 whole genome sequencing (WGS) analysis using Oxford Nanopore Technologies long-read data. Starting from raw FASTQ reads deposited in NCBI SRA, the pipeline performs quality control, reference-based genome assembly, variant calling, lineage assignment, and phylogenetic reconstruction — all containerized in Docker for full reproducibility and deployed on a Linux environment.

This repository is designed as an industrial-ready portfolio demonstrating the ability to build, document, and interpret a clinical-grade genomic surveillance pipeline for a high-priority respiratory pathogen.

---

## Biological and public health context

SARS-CoV-2 genomic surveillance is a cornerstone of pandemic response — identifying circulating lineages, tracking variant emergence, and informing public health decisions. Oxford Nanopore Technologies sequencing offers real-time, portable, and cost-effective genome sequencing particularly suited to resource-limited settings. This pipeline replicates the core analytical workflow used by genomic surveillance programs such as GISAID contributors and the Africa CDC pathogen genomics initiative.

---

## Dataset

| Parameter | Detail |
|-----------|--------|
| Organism | SARS-CoV-2 |
| Sequencing platform | Oxford Nanopore Technologies (ONT) |
| Number of samples | 20 |
| Data source | NCBI SRA |
| SRR accessions | SRR13020991, SRR13020993, SRR13020996, SRR13020998, SRR13020999, SRR13021001, SRR13021003, SRR13021005, SRR13021008, SRR13021015, SRR13021017, SRR13021020, SRR13021022, SRR13021024, SRR13021030, SRR13021032, SRR13021037, SRR13021044, SRR13021047, SRR13021048, SRR13021050 |
| Reference genome | SARS-CoV-2 Wuhan-Hu-1 (NC_045512.2) |

---

## Methodology

### Pipeline overview

```
Raw FASTQ (SRA)
      │
      ▼
  NanoStat QC                  ← nanostat.sh
      │
      ▼
  Reference Mapping            ← mapped.sh
  (Minimap2 → Samtools)
      │
      ▼
  Consensus Assembly           ← consensus.sh
  (Medaka)
      │
      ▼
  Variant Calling              ← variants.sh
  (Medaka / Clair3)
      │
      ▼
  Lineage Assignment           ← Lineage.sh
  (Pangolin)
      │
      ▼
  Phylogenetic Analysis        ← phylogenetic.sh
  (MAFFT + IQ-TREE2)
      │
      ▼
  Visualization and QC         ← SARVS CODE R.r
  (R + ggplot2)
```

### Tools and environment

| Tool | Purpose | Version |
|------|---------|---------|
| Docker | Full pipeline containerization and reproducibility | Latest |
| NanoStat | Raw read quality statistics for ONT data | — |
| Minimap2 | Long-read reference alignment (ONT mode) | — |
| Samtools | BAM sorting, indexing, coverage calculation | — |
| Medaka | Oxford Nanopore consensus polishing and variant calling | — |
| Pangolin | SARS-CoV-2 lineage assignment | — |
| MAFFT | Multiple sequence alignment of consensus genomes | — |
| IQ-TREE2 | Maximum likelihood phylogenetic tree inference | — |
| R + ggplot2 | Publication-grade QC visualization | — |

### Why Docker?

All tools are containerized using Docker to ensure:
- **Full reproducibility** — any researcher can re-run the exact same analysis on any machine
- **No dependency conflicts** — each tool runs in its own isolated environment
- **Industrial standard** — containerized pipelines are the standard in clinical genomics and public health surveillance workflows

### Shell scripts (CODES/)

| Script | Function |
|--------|----------|
| `nanostat.sh` | Runs NanoStat on raw FASTQ to report read length, quality, and total bases |
| `mapped.sh` | Aligns reads to NC_045512.2 using Minimap2, sorts and indexes with Samtools |
| `consensus.sh` | Generates polished consensus genome per sample using Medaka |
| `variants.sh` | Calls variants relative to reference, outputs VCF per sample |
| `Lineage.sh` | Assigns Pangolin lineage to each consensus sequence |
| `phylogenetic.sh` | Aligns consensus sequences with MAFFT and infers ML tree with IQ-TREE2 |

### R analysis (CODES/SARVS CODE R.r)

Produces three publication-grade QC figures using ggplot2:
- Genome coverage bar chart per sample
- Variant count bar chart per sample
- Sample quality assessment bubble plot (coverage vs variants vs genome length)

---

## Results and figure interpretation

### Figure 1 — SARS-CoV-2 Genome Coverage by Sample

`figures/01_genome_coverage.png`

All 20 samples achieved near-complete genome coverage (>99%) against the Wuhan-Hu-1 reference. Nineteen samples are coloured green (coverage 99.4–99.7%), indicating high-quality consensus assemblies suitable for downstream lineage assignment and phylogenetics. One sample (**SRR13021050_1**) is coloured orange (~99.1% coverage), flagging it as a lower-quality assembly relative to the cohort. This sample warrants closer inspection before clinical or surveillance reporting.

---

### Figure 2 — Variants per Sample

`figures/02_variant_counts.png`

Variant counts range from approximately 2,600 to 6,300 per sample. The bar chart is sorted in descending order and colour-coded by count (dark blue = high, light blue = low). **SRR13021050_1** has the highest variant count (~6,300), nearly double the median of the cohort. This is consistent with its lower coverage and confirms it as a statistical outlier — high variant counts in low-coverage assemblies often reflect sequencing errors rather than true biological variation, and this sample would be flagged for manual review in a surveillance context.

---

### Figure 3 — Sample Quality Assessment

`figures/03_quality_assessment.png`

This bubble plot integrates three quality dimensions simultaneously:
- **X-axis:** genome coverage (%)
- **Y-axis:** variant count
- **Bubble size:** consensus genome length (29,600–29,800 bp)
- **Colour:** coverage % (red = low, green = high)

The majority of samples cluster in the top-right quadrant (high coverage, moderate-to-high variants, full-length genomes), indicating good sequencing quality. **SRR13021050_1** is isolated in the top-left (red dot) — low coverage, highest variant count, confirming its outlier status across all three quality dimensions. This multi-dimensional QC plot provides a single-glance summary that would be used in a surveillance setting to flag samples requiring repeat sequencing.

---

### Figure 4 — Phylogenetic Tree

`figures/SARVS.png`

Maximum likelihood phylogenetic tree inferred by IQ-TREE2 from MAFFT-aligned consensus sequences. Bootstrap values (1000 replicates) are shown at each node.

**Key observations:**

- The tree resolves into **two major clades** with strong bootstrap support (100 at the root split), indicating two genetically distinct sub-populations within this sample set
- **Upper clade** (bootstrap 55–99): contains the majority of samples; further subdivides into two sub-clades, suggesting ongoing within-lineage diversification
- **Lower clade** (bootstrap 100): contains SRR13021008, SRR13021032, SRR13021022, SRR1302141, with strong internal support (bootstrap 75–100), indicating a closely related cluster potentially representing a distinct transmission chain
- **SRR13021050 and SRR13021047** form a well-supported sister pair (bootstrap 81) that branches basally — consistent with their divergent variant profiles observed in QC figures
- **SRR13020991** is the most divergent sample, branching as an outgroup to the entire remaining dataset

This tree structure would inform epidemiological investigation — the two major clades could represent distinct introductions or lineages circulating simultaneously in the sampled population.

---

## Repository structure

```
End-to-End-SARS-CoV-2-WGS/
├── README.md
├── CODES/
│   ├── nanostat.sh               QC statistics for raw ONT reads
│   ├── mapped.sh                 Reference alignment (Minimap2 + Samtools)
│   ├── consensus.sh              Consensus genome assembly (Medaka)
│   ├── variants.sh               Variant calling (VCF output)
│   ├── Lineage.sh                Pangolin lineage assignment
│   ├── phylogenetic.sh           MAFFT alignment + IQ-TREE2 ML tree
│   └── SARVS CODE R.r            R visualization pipeline (ggplot2)
├── figures/
│   ├── 01_genome_coverage.png    Coverage bar chart per sample
│   ├── 02_variant_counts.png     Variant count bar chart per sample
│   ├── 03_quality_assessment.png Multi-dimensional QC bubble plot
│   └── SARVS.png                 Maximum likelihood phylogenetic tree
└── analysis/
    ├── 01_genome_coverage.csv    Coverage data table
    ├── 02_variant_counts.csv     Variant count data table
    └── 03_quality_summary.csv    Combined QC summary table
```

---

## How to reproduce

### 1. Pull required Docker containers

```bash
# Minimap2 + Samtools
docker pull quay.io/biocontainers/minimap2:2.24--h7132678_1
docker pull quay.io/biocontainers/samtools:1.17--hd87286a_1

# Medaka
docker pull ontresearch/medaka:latest

# Pangolin
docker pull staphb/pangolin:latest

# MAFFT
docker pull quay.io/biocontainers/mafft:7.520--h031d066_1

# IQ-TREE2
docker pull quay.io/biocontainers/iqtree:2.2.2.6--h21ec9f0_0
```

### 2. Download raw data from SRA

```bash
# Example for one sample — repeat for all SRR accessions
fastq-dump --split-files SRR13021050 -O data/raw/
```

### 3. Run pipeline scripts in order

```bash
bash CODES/nanostat.sh
bash CODES/mapped.sh
bash CODES/consensus.sh
bash CODES/variants.sh
bash CODES/Lineage.sh
bash CODES/phylogenetic.sh
```

### 4. Generate figures in R

```r
source("CODES/SARVS CODE R.r")
```

---

## Key findings summary

| Sample | Coverage | Variants | Status |
|--------|----------|----------|--------|
| SRR13021050_1 | ~99.1% | ~6,300 | Outlier — flag for review |
| SRR13020991 | >99.5% | ~2,800 | Phylogenetically divergent |
| All others | 99.4–99.7% | 2,600–5,400 | Pass QC |

Nineteen of twenty samples passed all QC thresholds and resolved into two well-supported phylogenetic clades, consistent with two co-circulating SARS-CoV-2 lineages. One sample (SRR13021050_1) was flagged as a low-quality outlier across all quality metrics and should be excluded from or annotated in final surveillance reporting.

---

## Author

**[Your Name]**
Bioinformatics Training Program | [Institution]
[GitHub profile] | [Email or LinkedIn]
