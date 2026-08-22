# Hi-C End-to-End Analysis Pipeline

## From Raw FASTQ Data to 3D Genome Architecture

This repository contains an end-to-end Hi-C data analysis workflow, from raw paired-end FASTQ files through quality control, preprocessing, Hi-C contact map generation, and downstream 3D genome architecture analysis.

The workflow includes analysis of **TADs/domains, A/B compartments, compartment strength, Pearson correlations, chromatin loops, loop-anchor annotations, CTCF association, regulatory burden, statistical testing, and publication-quality visualizations**.

---

## Overview

The complete workflow follows the general structure below:

```text
Hi-C FASTQ
   │
   ▼
1. Raw-read QC
   │
   ▼
2. Adapter / quality trimming
   │
   ▼
3. Reference genome preparation
   │
   ├── hg38 genome
   ├── BWA index
   ├── MboI restriction sites
   └── chromosome sizes
   │
   ▼
4. Juicer
   │
   ├── mapping
   ├── filtering
   ├── duplicate handling
   ├── contact generation
   └── .hic file
   │
   ▼
5. Hi-C QC / validation
   │
   ├───────────────┬───────────────────┐
   ▼               ▼                   ▼
6. TADs       7. A/B compartments   8. Loops
                 │                     │
                 │                     ▼
                 │                 Loop anchors
                 │                     │
                 │                     ▼
                 │              Gene/promoter/
                 │              enhancer/CTCF
                 │                annotation
                 │                     │
                 │                     ▼
                 │                 Loop IDs
                 │                     │
                 └───────────┬─────────┘
                             ▼
                    9. Integrated
                       regulatory analysis
                             │
                             ▼
                        Final QC
```

---

## Scripts

### 1. `script_PART_1.txt`

This is the main end-to-end Hi-C analysis workflow.

It includes:

* Conda environment setup
* Installation of required software and dependencies
* Reference genome preparation
* hg38 genome indexing
* Restriction enzyme site generation
* Chromosome size file generation
* Raw FASTQ data processing
* FastQC quality assessment
* Adapter trimming and quality filtering using Trim Galore
* Hi-C processing using Juicer
* Hi-C contact map generation
* TAD detection using Arrowhead
* Contact matrix extraction
* A/B compartment analysis using FAN-C
* Compartment domain analysis
* Compartment strength analysis
* Pearson correlation analysis
* Chromatin loop detection
* Loop-anchor extraction
* Gene annotation
* Promoter annotation
* Enhancer and cCRE annotation
* CTCF annotation
* Loop-anchor classification
* CTCF association analysis
* Regulatory burden analysis
* Final quality control

### 2. `downstream_R_1.R`

This R script performs downstream analysis of annotated chromatin loops.

It includes:

* Basic data quality control
* Loop-level regulatory burden calculation
* Anchor functional classification
* Loop-pair distribution analysis
* CTCF distribution analysis
* CTCF association across loop classes
* Chi-square statistical testing
* CTCF percentage analysis
* Regulatory burden comparison across loop classes
* Genomic loop-distance analysis
* Statistical summary tables

### 3. `downstream_R_2.R`

This R script performs downstream analysis of large-scale 3D genome organization.

It includes:

* TAD/domain statistics
* TAD size distribution
* TAD size categorization
* A/B compartment analysis
* Eigenvector analysis
* A/B compartment distribution
* Compartment domain analysis
* Compartment strength analysis
* Pearson correlation analysis
* Statistical summaries


---

## Software and Tools

The workflow uses the following major software and tools:

* Conda / Mamba
* Juicer
* Juicer Tools
* FAN-C
* BEDTools
* BWA
* samtools
* FastQC
* Trim Galore
* SRA Toolkit
* GNU Parallel
* R

---

## R Dependencies

The downstream R analysis scripts use:

```r
tidyverse
rstatix
```

## Genome Assembly

**Reference genome:** hg38

The workflow prepares the reference genome for Hi-C analysis and generates the required supporting files, including:

* BWA index
* Restriction enzyme site file
* Chromosome sizes file

---

## Input Data

The pipeline uses:

### Raw sequencing data

* Paired-end Hi-C FASTQ files

### Reference data

* hg38 reference genome
* Restriction enzyme recognition sites
* Chromosome sizes
* GENCODE gene annotation
* ENCODE SCREEN cCRE annotation

### Hi-C analysis outputs

* Hi-C `.hic` contact maps
* TAD/domain files
* A/B compartment files
* Compartment domain files
* Compartment strength results
* Pearson correlation matrices
* Chromatin loop files
* Annotated loop-anchor files

---

## Main Analysis Workflow

### Step 1: Environment Setup

Create and configure the required Conda environments and install the software needed for Hi-C analysis.

### Step 2: Reference Genome Preparation

Download and prepare the hg38 reference genome.

This includes:

* Reference genome indexing with BWA
* Restriction site generation
* Chromosome size generation

### Step 3: Raw FASTQ Quality Control

Raw paired-end Hi-C sequencing reads are assessed using FastQC.

### Step 4: Read Trimming and Filtering

Adapters and low-quality bases are removed using Trim Galore.

### Step 5: Hi-C Processing

The processed reads are analyzed using Juicer to generate Hi-C contact maps.

### Step 6: TAD Analysis

Topologically associating domains are identified and analyzed.

The downstream analysis includes:

* Number of TADs
* TAD size statistics
* TAD size distribution
* TAD size categorization

### Step 7: A/B Compartment Analysis

A/B compartments are calculated using FAN-C.

The analysis includes:

* Compartment assignment
* Eigenvector values
* A/B compartment distribution
* Compartment domains
* Compartment strength
* Compartment enrichment

### Step 8: Pearson Correlation Analysis

Pearson correlation matrices are analyzed to examine chromatin interaction patterns.

### Step 9: Chromatin Loop Analysis

Chromatin loops are detected and loop anchors are extracted for downstream annotation.

### Step 10: Regulatory Annotation

Loop anchors are annotated for:

* Genes
* Promoters
* cCRE promoters
* Enhancers
* CTCF sites

### Step 11: CTCF Analysis

The workflow investigates CTCF association across different loop classes and performs statistical testing.

### Step 12: Regulatory Burden Analysis

Regulatory features are summarized at the loop level and compared across loop classes.

### Step 13: Visualization

The R scripts generate publication-quality figures, including:

* Anchor functional status
* Loop-pair distribution
* CTCF distribution
* CTCF status by loop class
* Regulatory burden across loop classes
* Loop-distance distribution
* Loop-distance comparison by class
* TAD size distribution
* TAD size boxplots
* A/B compartment distribution
* Eigenvector distribution
* Eigenvector values across compartments

---


### Tables

Output tables include summaries for:

* TAD statistics
* TAD size categories
* A/B compartment distribution
* Eigenvector values
* Loop-pair distribution
* Anchor functional status
* CTCF distribution
* CTCF association by loop class
* Regulatory burden
* Loop distance



