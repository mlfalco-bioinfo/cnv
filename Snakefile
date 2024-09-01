rule all:
    input:
        "results/segments/sample_segments.seg",
        "results/cnvkit/sample.cns",
        "results/cnvkit/scatter.png",
        "results/cnvkit/diagram.png",
        "results/gatk/sample_recalibrated.bam"

# Step 1: Perform Quality Control
rule fastqc:
    input:
        "data/sample_R1.fastq.gz",
        "data/sample_R2.fastq.gz"
    output:
        "results/fastqc/sample_R1_fastqc.html",
        "results/fastqc/sample_R2_fastqc.html"
    shell:
        "fastqc {input} -o results/fastqc"

# Step 2: Alignment
rule bwa_mem:
    input:
        R1="data/sample_R1.fastq.gz",
        R2="data/sample_R2.fastq.gz",
        ref="reference/reference.fasta"
    output:
        "results/aligned/sample_sorted.bam"
    params:
        rg="@RG\\tID:sample\\tSM:sample\\tPL:ILLUMINA"
    threads: 8
    shell:
        """
        bwa mem -t {threads} -R "{params.rg}" {input.ref} {input.R1} {input.R2} | \
        samtools sort -o {output}
        """

rule index_bam:
    input:
        "results/aligned/sample_sorted.bam"
    output:
        "results/aligned/sample_sorted.bam.bai"
    shell:
        "samtools index {input}"

# Step 3: Mark Duplicates and BQSR
rule mark_duplicates:
    input:
        "results/aligned/sample_sorted.bam"
    output:
        "results/gatk/sample_dedup.bam",
        "results/gatk/sample_metrics.txt"
    shell:
        """
        gatk MarkDuplicates -I {input} -O {output[0]} -M {output[1]}
        """

rule base_recalibration:
    input:
        bam="results/gatk/sample_dedup.bam",
        ref="reference/reference.fasta",
        known_sites="reference/known_sites.vcf"
    output:
        recal_data="results/gatk/sample_recal_data.table",
        bam_out="results/gatk/sample_recalibrated.bam"
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} -O {output.recal_data}
        gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {output.recal_data} -O {output.bam_out}
        """

# Step 4: CNV Calling with GATK
rule preprocess_intervals:
    input:
        ref="reference/reference.fasta",
        targets="reference/exome_targets.bed"
    output:
        "results/gatk/intervals.interval_list"
    shell:
        "gatk PreprocessIntervals -R {input.ref} -L {input.targets} -imr OVERLAPPING_ONLY -O {output}"

rule denoise_counts:
    input:
        counts="data/sample.counts.hdf5",
        pon="data/pon.hdf5"
    output:
        standardized="results/gatk/sample.standardizedCR.tsv",
        denoised="results/gatk/sample.denoisedCR.tsv"
    shell:
        "gatk DenoiseReadCounts -I {input.counts} --standardized-copy-ratios {output.standardized} --denoised-copy-ratios {output.denoised} --count-panel-of-normals {input.pon}"

rule call_segments:
    input:
        denoised="results/gatk/sample.denoisedCR.tsv"
    output:
        "results/segments/sample_segments.seg"
    shell:
        "gatk CallCopyRatioSegments --denoised-copy-ratios {input.denoised} --output {output}"

# Step 5: CNV Calling with CNVkit
rule cnvkit_coverage:
    input:
        bam="results/gatk/sample_recalibrated.bam",
        targets="reference/exome_targets.bed"
    output:
        "results/cnvkit/sample.cnn"
    shell:
        "cnvkit.py coverage {input.bam} {input.targets} -o {output}"

rule cnvkit_reference:
    input:
        normal_cnns="results/cnvkit/normal_1.cnn results/cnvkit/normal_2.cnn",
        ref_fasta="reference/reference.fasta"
    output:
        "results/cnvkit/my_reference.cnn"
    shell:
        "cnvkit.py reference {input.normal_cnns} -f {input.ref_fasta} -o {output}"

rule cnvkit_fix:
    input:
        sample_cnn="results/cnvkit/sample.cnn",
        ref_cnn="results/cnvkit/my_reference.cnn",
        targets="reference/exome_targets.bed"
    output:
        "results/cnvkit/sample.cnr"
    shell:
        "cnvkit.py fix {input.sample_cnn} {input.ref_cnn} {input.targets} -o {output}"

rule cnvkit_segment:
    input:
        cnr="results/cnvkit/sample.cnr"
    output:
        "results/cnvkit/sample.cns"
    shell:
        "cnvkit.py segment {input.cnr} -o {output}"

rule cnvkit_visualize:
    input:
        cnr="results/cnvkit/sample.cnr",
        cns="results/cnvkit/sample.cns"
    output:
        scatter="results/cnvkit/scatter.png",
        diagram="results/cnvkit/diagram.png"
    shell:
        """
        cnvkit.py scatter {input.cnr} -s {input.cns} -o {output.scatter}
        cnvkit.py diagram {input.cns} -o {output.diagram}
        """
