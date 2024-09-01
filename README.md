# CNV Analysis Pipeline using GATK and CNVkit

This repository contains a pipeline for performing Copy Number Variation (CNV) analysis on exome sequencing data using GATK and CNVkit. The pipeline is implemented using Snakemake, making it easy to run and manage the analysis.

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Pipeline Overview](#pipeline-overview)
4. [Running the Pipeline](#running-the-pipeline)
5. [Detailed Workflow Explanation](#detailed-workflow-explanation)
6. [References](#references)

## Introduction

Copy Number Variations (CNVs) are a form of structural variation where segments of the genome are duplicated or deleted. CNV analysis is essential for understanding genomic diversity and disease mechanisms, especially in cancer. This pipeline automates the process of detecting CNVs from exome sequencing data using GATK's `gCNV` and CNVkit, two robust tools for this purpose.

## Installation

### Prerequisites

- A Linux-based system with root access.
- Python 3.x installed.

### Installing the Necessary Tools

To install all the necessary tools, run the following script:

```
bash install_tools.sh

```

This script installs the following tools:

    BWA: For read alignment.
    SAMtools: For manipulating alignments in BAM format.
    GATK: For preprocessing and CNV calling.
    CNVkit: For CNV analysis from exome sequencing data.
    FastQC: For quality control of sequencing data.
    MultiQC: For aggregating results from FastQC and other tools.
    Trimmomatic: For trimming low-quality reads and adapters.
    Snakemake: For workflow management.
    ANNOVAR: (Optional, requires manual installation) For annotating CNVs with gene information.

Pipeline Overview

The pipeline includes the following steps:

    Quality Control: Assess the quality of the raw sequencing data using FastQC.
    Read Alignment: Align reads to a reference genome using BWA and sort the resulting BAM files.
    Duplicate Marking and Base Quality Score Recalibration: Mark duplicates and recalibrate base quality scores using GATK.
    CNV Calling using GATK: Use GATK's gCNV module to call CNVs.
    CNV Calling using CNVkit: An alternative approach for CNV analysis using CNVkit.
    Visualization: Generate plots for visualizing CNVs using CNVkit.
    Downstream Analysis: Annotate CNVs with gene information and assess their potential functional impact.

Running the Pipeline
1. Prepare the Reference Data

Download the necessary reference genome files and exome target regions. Place them in the reference/ directory.
2. Organize Your Data

Place your raw sequencing data in the data/ directory. Ensure your data files are named appropriately, for example:

    'sample_R1.fastq.gz'
    'sample_R2.fastq.gz'

3. Run Snakemake

Execute the pipeline with the following command:


```
snakemake --cores 8
``

This command will execute all steps in the workflow, utilizing 8 cores.
4. Output

The results will be saved in the results/ directory, organized by step.
Detailed Workflow Explanation
Step 1: Quality Control

Quality control is essential to ensure the raw sequencing data is of sufficient quality for downstream analysis.

    FastQC generates reports on various quality metrics for each FASTQ file.
    MultiQC aggregates these reports into a single HTML report for easier interpretation.

Step 2: Read Alignment

Reads are aligned to a reference genome using BWA-MEM. The aligned reads are then sorted using SAMtools.

    BWA-MEM aligns the reads to the reference genome.
    SAMtools is used to sort the resulting BAM files and create index files.

Step 3: Duplicate Marking and Base Quality Score Recalibration

These steps are crucial to improve the accuracy of variant calling:

    GATK MarkDuplicates identifies and marks duplicate reads.
    GATK BaseRecalibrator adjusts base quality scores based on known variant sites, improving the accuracy of subsequent analyses.

Step 4: CNV Calling using GATK

GATK's gCNV module is used to detect CNVs:

    PreprocessIntervals creates an interval list from the exome capture kit.
    DenoiseReadCounts normalizes the read counts using a panel of normals (PON).
    CallCopyRatioSegments calls CNVs based on the normalized data.

Step 5: CNV Calling using CNVkit

CNVkit provides an alternative approach to CNV analysis:

    CNVkit coverage calculates the coverage in target regions.
    CNVkit reference creates a reference profile using normal samples.
    CNVkit fix normalizes the sample data against the reference.
    CNVkit segment identifies CNVs by segmenting the data.

Step 6: Visualization

CNVkit can generate visualizations for easier interpretation of CNV data:

    Scatter plot shows the distribution of copy ratios.
    Diagram visualizes CNVs across the genome.

Step 7: Downstream Analysis

The CNV results can be further analyzed to assess their biological significance:

    Gene Annotation: Annotate CNVs with gene information using tools like ANNOVAR or CNVkit's built-in gene annotation functionality.
    Functional Impact: Assess the potential functional impact of the identified CNVs, especially if they overlap with known cancer-related genes.

References

    GATK: GATK Documentation
    CNVkit: CNVkit Documentation
    BWA: BWA Documentation
    SAMtools: SAMtools Documentation
    FastQC: FastQC Documentation
    MultiQC: MultiQC Documentation
    Trimmomatic: Trimmomatic Documentation
