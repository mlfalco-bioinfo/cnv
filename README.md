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

```bash
bash install_tools.sh
