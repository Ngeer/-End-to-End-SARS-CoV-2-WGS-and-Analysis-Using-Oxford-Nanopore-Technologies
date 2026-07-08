# End-to-End SARS-CoV-2 Whole Genome Sequencing and Analysis Using Oxford Nanopore Technologies

> **Oxford Nanopore Technologies (ONT) sequencing** | 20 SARS-CoV-2 samples | Pipeline: Docker + R | Reference: MN908947.3

---

## Project overview

This project demonstrates a complete, reproducible end-to-end bioinformatics pipeline for SARS-CoV-2 whole genome sequencing (WGS) analysis using Oxford Nanopore Technologies long-read data. Starting from raw FASTQ reads, the pipeline performs quality control, reference-based read mapping, primer trimming, consensus genome assembly, variant calling, lineage assignment, and phylogenetic reconstruction — every step containerized in Docker for full reproducibility.

This repository is designed as an industrial-ready portfolio demonstrating the ability to build, document, and interpret a clinical-grade genomic surveillance pipeline for a high-priority respiratory pathogen.

---

## Biological and public health context

SARS-CoV-2 genomic surveillance is a cornerstone of pandemic response — identifying circulating lineages, tracking variant emergence, and informing public health decisions. Oxford Nanopore Technologies sequencing offers real-time, portable, and cost-effective genome sequencing particularly suited to resource-limited settings and point-of-care deployment. This pipeline replicates the core analytical workflow used by genomic surveillance programs contributing to GISAID and the Africa CDC pathogen genomics initiative.

---

## Dataset

| Parameter | Detail |
|-----------|--------|
| Organism | SARS-CoV-2 |
| Sequencing platform | Oxford Nanopore Technologies (ONT) |
| Number of samples | 20 |
| Reference genome | SARS-CoV-2 Wuhan-Hu-1 (MN908947.3) |
| Primer scheme | ARTIC nCoV-2019 V3 |
| Raw data format | FASTQ.gz (ONT_downloads/) |

---

## Methodology

### Full pipeline overview

```
Raw ONT FASTQ reads
        │
        ▼
   NanoStat QC                    nanostat.sh
   (per-sample read stats)
        │
        ▼
   Reference Mapping               mapped.sh
   Minimap2 (-ax map-ont)
   → Samtools sort → .sorted.bam
        │
        ▼
   Primer Trimming                 variants.sh
   iVar trim (ARTIC V3 BED)
   → .trimmed.sorted.bam
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
   Consensus Assembly              variants.sh
   consensus.sh                    iVar variants
   samtools mpileup | ivar         (samtools mpileup | ivar variants)
   consensus (-t 0.5 -m 10)        → .variants.tsv per sample
   → .fa per sample
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
   Lineage Assignment              Phylogenetic Analysis
   Lineage.sh                      phylogenetic.sh
   Pangolin → lineage report       MAFFT --auto → aligned.fasta
   Nextclade → clade + mutations   IQ-TREE2 GTR+G -bb 1000
        │                                      │
        └──────────────┬───────────────────────┘
                       ▼
              R Visualization
              SARVS CODE R.r
              ggplot2 QC figures
```

---

### Step-by-step methodology

#### Step 1 — Raw read QC (`nanostat.sh`)

NanoStat is run per sample on raw FASTQ.gz files to report read count, total bases, mean read length, and mean read quality. Results are written to `qc_results/`. This step identifies samples with insufficient read depth or poor quality before committing compute resources to alignment.

```bash
# Docker image used
quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0
```

---

#### Step 2 — Reference mapping (`mapped.sh`)

Each sample's FASTQ.gz is aligned to the SARS-CoV-2 Wuhan-Hu-1 reference genome (MN908947.3) using Minimap2 in ONT mode (`-ax map-ont`). The SAM output is piped directly into Samtools sort to produce a sorted BAM file, avoiding intermediate SAM disk writes.

```bash
# Docker images used
quay.io/biocontainers/minimap2:2.24--h7132678_1
quay.io/biocontainers/samtools:1.17--h00cdaf9_0

# Key parameters
minimap2 -ax map-ont    # ONT long-read preset
samtools sort           # coordinate-sorted BAM output
```

---

#### Step 3 — Primer trimming and variant calling (`variants.sh`)

ARTIC V3 primers are trimmed from aligned reads using iVar trim with a minimum length of 30 bp (`-m 30`). The trimmed BAM is re-sorted, then variants are called using `samtools mpileup` piped into `ivar variants` against MN908947.3. This produces a TSV file of variant positions per sample.

```bash
# Docker images used
quay.io/biocontainers/samtools:1.17--h00cdaf9_0
quay.io/biocontainers/ivar:1.4.4--h077b44d_0

# Key parameters
ivar trim -m 30 -q 0                    # minimum length 30bp, no quality filter
samtools mpileup -A -d 0 -Q 0           # all reads, unlimited depth, no base quality filter
ivar variants -r MN908947.3.fasta       # variant calling against reference
```

---

#### Step 4 — Consensus genome assembly (`consensus.sh`)

Consensus genomes are built using `samtools mpileup` piped into `ivar consensus` with a majority threshold of 0.5 (50% of reads must support a base call) and a minimum depth of 10x (`-m 10`). Positions below 10x depth are masked with N. Output is one FASTA consensus per sample.

```bash
# Docker image used
quay.io/biocontainers/ivar:1.4.4--h077b44d_0

# Key parameters
ivar consensus -t 0.5    # majority base call threshold
ivar consensus -m 10     # minimum depth of 10x; below this → N masking
```

---

#### Step 5 — Lineage assignment (`Lineage.sh`)

Each consensus FASTA is processed by two complementary lineage classifiers:

- **Pangolin** (staphb/pangolin): assigns WHO/Pango lineage designation (e.g. B.1.1.7, BA.2) based on the most current lineage nomenclature
- **Nextclade** (nextstrain/nextclade): assigns Nextstrain clade, identifies amino acid substitutions relative to reference, and flags quality issues per sample

Running both tools provides cross-validated lineage calls and richer mutation-level annotation than either tool alone.

```bash
# Docker images used
staphb/pangolin:latest
nextstrain/nextclade:latest
```

---

#### Step 6 — Phylogenetic analysis (`phylogenetic.sh`)

All consensus FASTA files are concatenated into `all_consensus.fasta`, aligned with MAFFT (`--auto` mode selects the appropriate algorithm for the dataset size), and a maximum likelihood phylogenetic tree is inferred by IQ-TREE2 using the GTR+G substitution model with 1000 ultrafast bootstrap replicates (`-bb 1000`).

```bash
# Docker images used
staphb/mafft:latest
staphb/iqtree2:latest

# Key parameters
mafft --auto                        # automatic algorithm selection
iqtree2 -m GTR+G -bb 1000 -nt AUTO # GTR+Gamma model, 1000 bootstrap, auto threads
```

---

#### Step 7 — Visualization (`SARVS CODE R.r`)

Three publication-grade QC figures are produced in R using ggplot2, integrating coverage, variant, and genome length data across all samples.

---

### Why Docker for every step?

Every tool in this pipeline runs inside a Docker container with a pinned image version. This guarantees:

- **Exact reproducibility** — the same tool version runs on any machine, any operating system, any time
- **No dependency conflicts** — Minimap2, iVar, Pangolin, and IQ-TREE2 have incompatible Python/system dependencies; Docker isolates each one
- **Clinical-grade auditability** — the exact container image tag is logged in each script, creating a complete audit trail of the software environment used
- **Portability** — the pipeline runs identically on a local Linux workstation, an HPC cluster, or a cloud VM

---

## Results and figure interpretation

### Figure 1 — SARS-CoV-2 Genome Coverage by Sample

`figures/01_genome_coverage.png`

All 20 samples achieved near-complete genome coverage against MN908947.3. Nineteen samples show coverage of 99.4–99.7% (green), indicating high-quality consensus assemblies with minimal N-masking. **SRR13021050_1** is the single outlier at approximately 99.1% coverage (orange bar), reflecting a greater proportion of low-depth positions masked as N in the consensus. In a clinical surveillance context, this sample would require either repeat sequencing or explicit annotation as lower confidence before submission to GISAID.

---

### Figure 2 — Variants per Sample

`figures/02_variant_counts.png`

Variant counts across the cohort range from approximately 2,600 to 6,300 per sample. **SRR13021050_1** carries the highest variant burden (~6,300), roughly double the cohort median. This is directly consistent with its lower genome coverage — regions of low or absent coverage produce ambiguous mpileup pileups that iVar may interpret as variant sites, inflating variant counts artificially. This concordance between low coverage and high variant count is the expected signature of a sequencing quality problem, not a biologically distinct viral strain, and further supports flagging this sample for review.

The remaining samples show a gradual decrease from ~5,400 (SRR13021048_1) to ~2,600 (SRR13021024_1), a range consistent with natural between-sample variation in lineage-specific mutation burden.

---

### Figure 3 — Sample Quality Assessment

`figures/03_quality_assessment.png`

This bubble plot integrates three quality dimensions simultaneously:
- **X-axis:** genome coverage (%)
- **Y-axis:** variant count
- **Bubble size:** consensus genome length (29,600–29,800 bp)
- **Colour:** coverage % (red = low, yellow = intermediate, green = high)

The majority of samples cluster in the right-centre region (coverage 99.4–99.7%, variant counts 2,600–5,400, genome lengths near complete), confirming cohort-wide sequencing quality. **SRR13021050_1** is isolated as a red dot in the upper-left corner — the lowest coverage, highest variant count, and visually separated from the rest of the cohort across all three dimensions simultaneously. This multi-axis separation confirms it is a genuine quality outlier and not an artefact of any single metric. This type of integrated QC plot is standard in clinical genomics pipelines before downstream reporting.

---

### Figure 4 — Maximum Likelihood Phylogenetic Tree

`figures/SARVS.png`

Maximum likelihood tree inferred by IQ-TREE2 (GTR+G, 1000 ultrafast bootstrap replicates). Bootstrap support values are shown at each internal node.

**Key observations:**

The tree resolves into **two major clades** separated by a deep internal node with bootstrap support of 100, indicating two genetically distinct sub-populations within this sample set — consistent with either two co-circulating lineages or two separate introduction events.

**Upper clade (bootstrap 55–99):** Contains the majority of samples and further subdivides into sub-clades with strong internal support (bootstrap 88–100). The tight clustering of SRR13021005/SRR13021020 (bootstrap 100) and SRR13021001/SRR13021037 (bootstrap 96) suggests recent common ancestry, potentially indicating direct or indirect transmission links.

**Lower clade (bootstrap 100):** SRR13021008, SRR13021032, SRR13021022, and SRR1302141 form a strongly supported cluster (bootstrap 75–100), representing a closely related group likely derived from a distinct transmission chain or introduction.

**SRR13021050 and SRR13021047** form a well-supported sister pair (bootstrap 81) that branches at the base of the entire tree — phylogenetically the most divergent samples in the dataset. This is consistent with their outlier variant profiles: SRR13021050's high variant count contributes to a longer branch length, placing it basal relative to other samples.

**SRR13020991** branches as the most divergent single sample, separating from all others at the root — a candidate outgroup that may represent an earlier or genetically distinct lineage.

In a surveillance context, this tree would be used to: (1) confirm lineage co-circulation, (2) identify potential transmission clusters for epidemiological follow-up, and (3) flag phylogenetically anomalous samples for deeper investigation.

---

## Repository structure

```
End-to-End-SARS-CoV-2-WGS/
├── README.md
├── CODES/
│   ├── nanostat.sh               NanoStat QC on raw ONT FASTQ
│   ├── mapped.sh                 Minimap2 alignment + Samtools sort
│   ├── variants.sh               iVar primer trim + variant calling
│   ├── consensus.sh              iVar consensus genome assembly
│   ├── Lineage.sh.sh             Pangolin + Nextclade lineage assignment
│   ├── phylogenetic .sh.sh       MAFFT alignment + IQ-TREE2 ML tree
│   └── SARVS CODE R.r            R ggplot2 QC visualization
├── figures/
│   ├── 01_genome_coverage.png    Coverage bar chart
│   ├── 02_variant_counts.png     Variant count bar chart
│   ├── 03_quality_assessment.png Multi-dimensional QC bubble plot
│   └── SARVS.png                 IQ-TREE2 phylogenetic tree
└── analysis/
    ├── 01_genome_coverage.csv    Coverage data table
    ├── 02_variant_counts.csv     Variant count data table
    └── 03_quality_summary.csv    Combined QC summary
```

---

## How to reproduce

### 1. Pull all Docker containers

```bash
docker pull quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0
docker pull quay.io/biocontainers/minimap2:2.24--h7132678_1
docker pull quay.io/biocontainers/samtools:1.17--h00cdaf9_0
docker pull quay.io/biocontainers/ivar:1.4.4--h077b44d_0
docker pull staphb/pangolin:latest
docker pull nextstrain/nextclade:latest
docker pull staphb/mafft:latest
docker pull staphb/iqtree2:latest
```

### 2. Set up directory structure

```bash
mkdir -p /home/stuart/viral/{ONT_downloads,mapped,consensus,variants,qc_results,pangolin_results,nextclade_results,phylo}
```

### 3. Add reference genome and primer scheme

```bash
# Place MN908947.3.fasta in:
/home/stuart/viral/reference/MN908947.3.fasta

# Place ARTIC V3 BED file in:
/home/stuart/viral/primer_schemes/nCoV-2019/V3/nCoV-2019.bed
```

### 4. Run pipeline in order

```bash
bash CODES/nanostat.sh
bash CODES/mapped.sh
bash CODES/variants.sh       # includes primer trimming
bash CODES/consensus.sh
bash CODES/Lineage.sh.sh
bash CODES/phylogenetic.sh.sh
```

### 5. Generate QC figures

```r
source("CODES/SARVS CODE R.r")
```

---

## Key findings summary

| Sample | Coverage | Variant count | Phylogenetic position | Status |
|--------|----------|---------------|-----------------------|--------|
| SRR13021050_1 | ~99.1% | ~6,300 | Basal outgroup pair | Flag for review |
| SRR13020991 | >99.5% | ~2,800 | Most divergent single sample | Monitor |
| SRR13021047 | >99.5% | ~2,700 | Basal pair with SRR13021050 | Pass QC |
| All others | 99.4–99.7% | 2,600–5,400 | Two well-supported clades | Pass QC |

19 of 20 samples passed all QC thresholds. The cohort resolves into two phylogenetically distinct clades consistent with co-circulating lineages. SRR13021050_1 is flagged as a quality outlier across coverage, variant count, and phylogenetic position and should be excluded from or explicitly annotated in final surveillance reporting.

---
