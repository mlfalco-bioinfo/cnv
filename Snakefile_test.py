import os
import glob

# Define directory paths (modify as needed)
DATA_DIR = "data"
RESULTS_DIR = "results"
REFERENCE_DIR = "reference"

# Use `glob_wildcards` for robust sample name retrieval
SAMPLES = [
    sample.split("_R1")[0]
    for sample in sorted(glob.glob(f"{DATA_DIR}/{{sample}}_R1.fastq.gz"))
]

# Step 1: Quality Control (FastQC)
rule fastqc:
    """
    Perform FastQC quality control on raw sequencing data (FASTQ files).
    """
    input:
        R1 = f"{DATA_DIR}/{{sample}}_R1.fastq.gz",
        R2 = f"{DATA_DIR}/{{sample}}_R2.fastq.gz"
    output:
        R1_out = f"{RESULTS_DIR}/fastqc/{{sample}}_R1_fastqc.html",
        R2_out = f"{RESULTS_DIR}/fastqc/{{sample}}_R2_fastqc.html"
    shell:
        "fastqc {input.R1} {input.R2} -o {RESULTS_DIR}/fastqc"
        
# Step 2: Alignment (BWA mem)
rule bwa_mem:
    """
    Align sequencing reads to a reference genome using BWA mem.
    """
    input:
        R1 = f"{DATA_DIR}/{{sample}}_R1.fastq.gz",
        R2 = f"{DATA_DIR}/{{sample}}_R2.fastq.gz",
        ref = f"{REFERENCE_DIR}/reference.fasta"
    output:
        bam = f"{RESULTS_DIR}/aligned/{{sample}}_sorted.bam"
    params:
        rg = "@RG\\tID:{{sample}}\\tSM:{{sample}}\\tPL:ILLUMINA"
    threads: 8
    shell:
        """
        bwa mem -t {threads} -R "{params.rg}" {input.ref} {input.R1} {input.R2} | \
        samtools sort -o {output.bam}
        """

rule index_bam:
    """
    Index the sorted BAM file for efficient access.
    """
    input:
        bam = f"{RESULTS_DIR}/aligned/{{sample}}_sorted.bam"
    output:
        bai = f"{RESULTS_DIR}/aligned/{{sample}}_sorted.bam.bai"
    shell:
        "samtools index {input.bam}"
        
# Step 3: Mark Duplicates and Base Quality Score Recalibration (BQSR)
rule mark_duplicates:
    """
    Identify and remove duplicate reads in the BAM file.
    """
    input:
        bam = f"{RESULTS_DIR}/aligned/{{sample}}_sorted.bam"
    output:
        dedup_bam = f"{RESULTS_DIR}/gatk/{{sample}}_dedup.bam",
        metrics = f"{RESULTS_DIR}/gatk/{{sample}}_metrics.txt"
    shell:
        "gatk MarkDuplicates -I {input.bam} -O {output.dedup_bam} -M {output.metrics}"

rule base_recalibration:
    """
    Recalibrate base quality scores using GATK BaseRecalibrator and ApplyBQSR.
    """
    input:
        bam = f"{RESULTS_DIR}/gatk/{{sample}}_dedup.bam",
        ref = f"{REFERENCE_DIR}/reference.fasta",
        known_sites = f"{REFERENCE_DIR}/known_sites.vcf"
    output:
        recal_data = f"{RESULTS_DIR}/gatk/{{sample}}_recal_data.table",
        bam_out = f"{RESULTS_DIR}/gatk/{{sample}}_recalibrated.bam"
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} -O {output.recal_data}
        gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {output.recal_data} -O {output.bam_out}
        """

# Step 4: CNV Calling with GATK
rule preprocess_intervals:
    """
    Prepare a target interval list for CNV calling.
    """
    input:
        ref = f"{REFERENCE_DIR}/reference.fasta",
        targets = f"{REFERENCE_DIR}/exome_targets.bed"
    output:
        intervals = f"{RESULTS_DIR}/gatk/intervals.interval_list"
    shell:
        "gatk PreprocessIntervals -R {input.ref} -L {input.targets} -imr OVERLAPPING_ONLY -O {output.intervals}"

rule denoise_counts:
    """
    Denoise read counts using GATK DenoiseReadCounts.
    """
    input:
        counts = f"{DATA_DIR}/{{sample}}.counts.hdf5",
        pon = f"{DATA_DIR}/pon.hdf5"
    output:
        standardized = f"{RESULTS_DIR}/gatk/{{sample}}.standardizedCR.tsv",
        denoised = f"{RESULTS_DIR}/gatk/{{sample}}.denoisedCR.tsv"
    shell:
        "gatk DenoiseReadCounts -I {input.counts} --standardized-copy-ratios {output.standardized} --denoised-copy-ratios {output.denoised} --count-panel-of-normals {input.pon}"

rule call_segments:
    """
    Call CNV segments using GATK CallCopyRatioSegments.
    """
    input:
        denoised = f"{RESULTS_DIR}/gatk/{{sample}}.denoisedCR.tsv"
    output:
        segments = f"{RESULTS_DIR}/segments/{{sample}}_segments.seg"
    shell:
        "gatk CallCopyRatioSegments --denoised-copy-ratios {input.denoised} --output {output.segments}"


# Step 5: CNV Calling with CNVkit
rule cnvkit_coverage:
    """
    Calculate coverage for the sample using CNVkit.
    """
    input:
        bam = f"{RESULTS_DIR}/gatk/{{sample}}_recalibrated.bam",
        targets = f"{REFERENCE_DIR}/exome_targets.bed"
    output:
        cnn = f"{RESULTS_DIR}/cnvkit/{{sample}}.cnn"
    shell:
        "cnvkit.py coverage {input.bam} {input.targets} -o {output.cnn}"


rule cnvkit_reference:
    """
    Create a reference copy number profile using CNVkit.
    """
    input:
        normal_cnns = expand(f"{RESULTS_DIR}/cnvkit/{{normal}}.cnn", normal=["normal_1", "normal_2"]),
        ref_fasta = f"{REFERENCE_DIR}/reference.fasta"
    output:
        ref_cnn = f"{RESULTS_DIR}/cnvkit/my_reference.cnn"
    shell:
        "cnvkit.py reference {input.normal_cnns} -f {input.ref_fasta} -o {output.ref_cnn}"


rule cnvkit_fix:
    """
    Correct for GC bias and other artifacts using CNVkit.
    """
    input:
        sample_cnn = f"{RESULTS_DIR}/cnvkit/{{sample}}.cnn",
        ref_cnn = f"{RESULTS_DIR}/cnvkit/my_reference.cnn",
        targets = f"{REFERENCE_DIR}/exome_targets.bed"
    output:
        cnr = f"{RESULTS_DIR}/cnvkit/{{sample}}.cnr"
    shell:
        "cnvkit.py fix {input.sample_cnn} {input.ref_cnn} {input.targets} -o {output.cnr}"


rule cnvkit_segment:
    """
    Segment the data to identify CNV regions using CNVkit.
    """
    input:
        cnr = f"{RESULTS_DIR}/cnvkit/{{sample}}.cnr"
    output:
        cns = f"{RESULTS_DIR}/cnvkit/{{sample}}.cns"
    shell:
        "cnvkit.py segment {input.cnr} -o {output.cns}"

rule cnvkit_visualize:
    """
    Visualize CNV regions using CNVkit.
    """
    input:
        cnr = f"{RESULTS_DIR}/cnvkit/{{sample}}.cnr",
        cns = f"{RESULTS_DIR}/cnvkit/{{sample}}.cns"
    output:
        scatter = f"{RESULTS_DIR}/cnvkit/{{sample}}_scatter.png",
        diagram = f"{RESULTS_DIR}/cnvkit/{{sample}}_diagram.png"
    shell:
        """
        cnvkit.py scatter {input.cnr} -s {input.cns} -o {output.scatter}
        cnvkit.py diagram {input.cns} -o {output.diagram}
        """
