## Load config values
configfile: "config/input.yaml"
configfile: "config/parameter.yaml"
configfile: "config/resource.yaml"

# Get ids
input_dir = config["input_dir"]
sample_id = config["sample_id"]
cell_ids = config['cell_ids']

rule all:
    input: expand("output/{sample_id}/{cell_id}/0_qc/multiqc_report.html", cell_id = cell_ids, sample_id = [sample_id])
    localrule: True

###############
# QC Metrics for fastq
###############

# Generate fastq screen config file
rule fastq_screen_config:
    output: "output/fastq.config"
    localrule: True
    script:
        "src/generateFastqScreen.py"

# Run fastq_screen to detect contamination
rule fastq_screen: 
    input: 
        fq1 = input_dir + "/{cell_id}_1.fastq.gz",
        fq2 = input_dir + "/{cell_id}_2.fastq.gz",
        config = "output/fastq.config"
    output:
        tagged_fq1 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_1.tagged.fastq.gz",
        tagged_fq2 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_2.tagged.fastq.gz"
        metric_fq1 = "output/{sample_id}/{cell_id}/0_qc/fastq_screen/{cell_id}_1_screen.txt",
        metric_fq2 = "output/{sample_id}/{cell_id}/0_qc/fastq_screen/{cell_id}_2_screen.txt"
    singularity: "docker://quay.io/biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0"
    threads: config["fastq_screen"]["threads"]
    resources:
        mem_mb = config["fastq_screen"]["memory"],
        cpus = config["fastq_screen"]["threads"]
    log: "output/{sample_id}/{cell_id}/log/fastq_screen.log"
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastq_screen
        
        fastq_screen \
            --aligner bowtie2 \
            --conf {input.config} \
            --outdir output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastp \
            --tag \
            --thread {threads} \
            {input.fq1} {input.fq2} &> {log}
        
        mv output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastp/*tagged.fastq.gz \
            output/{wildcards.sample_id}/{wildcards.cell_id}/1_preprocessing
        """

rule fastqc: 
    input: 
        fq1 = input_dir + "/{cell_id}_1.fastq.gz",
        fq2 = input_dir + "/{cell_id}_2.fastq.gz",
    output:
        zip_1 = "output/{sample_id}/{cell_id}/0_qc/fastqc/{cell_id}_1_fastqc.zip",
        zip_2 = "output/{sample_id}/{cell_id}/0_qc/fastqc/{cell_id}_2_fastqc.zip"
    threads: config["fastqc"]["threads"]
    resources:
        mem_mb = config["fastqc"]["memory"],
        cpus = config["fastqc"]["threads"]
    log: "output/{sample_id}/{cell_id}/log/fastqc.log"
    singularity: "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastqc
        fastqc \
            --threads {threads} \
            --memory {resources.mem_mb} \
            --outdir output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastqc \
            {input.fq1} {input.fq2} &> {log}
        """


###############
# Preprocess fastq
###############

# Trims adapter sequence and poly-G tail
rule fastp:
    input:
        fq1 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_1.tagged.fastq.gz",
        fq2 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_2.tagged.fastq.gz"
    output:
        fq1 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_1_trimmed.fastq.gz",
        fq2 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_2_trimmed.fastq.gz",
        html = "output/{sample_id}/{cell_id}/0_qc/fastp/{cell_id}_fastp.html",
        json = "output/{sample_id}/{cell_id}/0_qc/fastp/{cell_id}_fastp.json"
    params: 
        adapter_sequence_1 = config["adapter_sequence_1"],
        adapter_sequence_2 = config["adapter_sequence_2"]
    singularity: "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    threads: config["fastp"]["threads"]
    resources:
        mem_mb = config["fastp"]["memory"],
        cpus = config["fastp"]["threads"]
    log: "output/{sample_id}/{cell_id}/log/fastp.log"
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/1_preprocessing
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastp

        fastp \
            -i {input.fq1} \
            -I {input.fq2} \
            -o {output.fq1} \
            -O {output.fq2} \
            -h {output.html} \
            -j {output.json} \
            --adapter_sequence {params.adapter_sequence_1} \
            --adapter_sequence_r2 {params.adapter_sequence_2} \
            --cut_right \
            -w {threads} \
            --trim_poly_g &> {log}
        """


###############
# Align reads
###############

# Align to minimap2 or bwa (Recommended minimap2)
if config['aligner'] == "minimap2":
    rule minimap2:
        input:
            fq1 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_1_trimmed.fastq.gz",
            fq2 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_2_trimmed.fastq.gz",
            reference = config["reference_genome"]
        output:
            bam = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.bam"
        singularity: "docker://quay.io/biocontainers/minimap2:2.26--he4a0461_2"
        threads: config["minimap2"]["threads"]
        resources:
            mem_mb=config["minimap2"]["memory"],
            cpus=config["minimap2"]["threads"]
        log: "output/{sample_id}/{cell_id}/log/minimap2.log"
        shell:
            """
            mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/2_alignment
            minimap2 \
                -t {threads} \
                -ax sr \
                {input.reference} \
                {input.fq1} {input.fq2} > {output.bam} 2> {log} 
            """
elif config['aligner'] == "bwa":
    rule bwa:
        input:
            fq1 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_1_trimmed.fastq.gz",
            fq2 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_2_trimmed.fastq.gz",
            reference = config["reference_genome"]
        output:
            bam = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.bam"
        singularity: "docker://quay.io/biocontainers/bwa:0.7.8--he4a0461_9"
        threads: config["bwa"]["threads"]
        log: "output/{sample_id}/{cell_id}/log/bwa.log"
        shell:
            """
            mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/2_alignment
            # TODO 
            """

# Sort and index BAM files
rule samtools_sort:
    input: "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.bam"
    output: 
        sorted_bam = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.sorted.bam",
        index_bam = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.sorted.bam.bai"
    singularity: "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    log: "output/{sample_id}/{cell_id}/log/sort_index_bam.log"
    threads: config["samtools_sort"]["threads"]
    resources:
        mem_mb=config["samtools_sort"]["memory"],
        cpus=config["samtools_sort"]["threads"]
    shell:
        """
        samtools sort -@ {threads} -o {output.sorted_bam} {input} &> {log}
        samtools index -@ {threads} {output.sorted_bam} &> {log}
        """


###############
# QC Metrics for BAM
###############

# Mark duplicate reads in BAM files
rule mark_duplication:
    input: "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.sorted.bam"
    output: 
        bam = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.markDup.bam",
        metric = "output/{sample_id}/{cell_id}/0_qc/markDup/{cell_id}.markDup_metric.txt"
    singularity: "docker://quay.io/biocontainers/picard-slim_2.27.4--hdfd78af_0.sif"
    log: "output/{sample_id}/{cell_id}/log/sort_index_bam.log"
    threads: 1
    resources:
        mem_mb = config["mark_duplication"]["memory"],
        cpus=1
    shell:
        """
        mkdir -p output/{sample_id}/{cell_id}/0_qc/markDup
        
        picard \
            -Xmx${resource.mem_mb}M \
            MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metric} &> {log}
        """

finalBAM = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.markDup.bam"

# Flagstat mtric
rule flagstat:
    input: finalBAM
    output: "output/{sample_id}/{cell_id}/0_qc/flagstat/{cell_id}.flagstat.txt"
    singularity: "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    log: "output/{sample_id}/{cell_id}/log/flagstat.log"
    threads: config["flagstat"]["threads"]
    resources:
        mem_mb = config["flagstat"]["memory"],
        cpus=config["flagstat"]["threads"]
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/flagstat
        samtools flagstat \
            -@ {threads} \
            {input} > {output} 2> {log}
        """

# Insert Size metric
rule CollectInsertSizeMetrics:
    input: finalBAM
    output:
        metric = "output/{sample_id}/{cell_id}/0_qc/insertsize/{cell_id}.insertSize.txt",
        histogram = "output/{sample_id}/{cell_id}/0_qc/insertsize/{cell_id}.insertSize.png",
    singularity: "docker://quay.io/biocontainers/picard-slim_2.27.4--hdfd78af_0.sif"
    log: "output/{sample_id}/{cell_id}/log/CollectInsertSizeMetrics.log"
    threads: 1
    resources:
        mem_mb = config["CollectInsertSizeMetrics"]["memory"],
        cpus=1
    shell:
        """
        picard \
            -Xmx${resource.mem_mb}M \
            CollectInsertSizeMetrics \
            -INPUT {input} \
            -OUTPUT {output.metric} \
            -Histogram_FILE {output.histogram} \
            -VALIDATION_STRINGENCY LENIENT \
            -MAX_RECORDS_IN_RAM 150000 &> {log}
        """

# WGS Metrics
rule CollectWgsMetrics:
    input: 
        bam = finalBAM,
        reference = config["reference_genome"]
    output:
        metric = "output/{sample_id}/{cell_id}/0_qc/WgsMetrics/{cell_id}.wgsMetric.txt",
    singularity: "docker://quay.io/biocontainers/picard-slim_2.27.4--hdfd78af_0.sif"
    log: "output/{sample_id}/{cell_id}/log/CollectWgsMetrics.log"
    threads: 1
    resources:
        mem_mb = config["CollectWgsMetrics"]["memory"],
        cpus=1
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/WgsMetrics

        picard \
            -Xmx${resource.mem_mb}M \
            CollectWgsMetrics \
            INPUT={input.bam} \
            OUTPUT={output.metric} \
            REFERENCE_SEQUENCE={input.reference} \
            COVERAGE_CAP=500 \
            MAX_RECORDS_IN_RAM=150000 &> {log}
        """
# Mitochondria Coverage
rule getMitoBed:
    input: config["reference_genome"]+".fai"
    output: temp("output/{sample_id}/mt.bed")
    localrule: True
    shell: 
        """
        awk -v OFS="\t" '$1 == "chrM" {{print $1, "0", $2}}' {input} > {output}
        """

rule removeDup:
    input: finalBAM
    output: "output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.rmDup.bam"
    singularity: "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
    log: "output/{sample_id}/{cell_id}/log/removeDup.log"
    threads: 1
    resources:
        mem_mb=config["removeDup"]["memory"],
        cpus=1
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/mtDNA/
        samtools view -uF 0x400 {input} > {output} 2> {log}
        """

rule bamToBed:
    input: "output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.rmDup.bam"
    output: "output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.bed"
    singularity: "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
    log: "output/{sample_id}/{cell_id}/log/bamToBed.log"
    threads: 1
    resources:
        mem_mb=config["bamToBed"]["memory"],
        cpus=1
    shell:
        """
        bedtools bamtobed -i {input} > {output} 2> {log} 
        """

rule getMtCoverage:
    input: 
        bed = "output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.bed",
        mt = "output/{sample_id}/mt.bed"
    output: "output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.mtCov.bed"
    singularity: "docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0"
    log: "output/{sample_id}/{cell_id}/log/getMtCoverage.log"
    threads: 1
    resources:
        mem_mb=config["getMtCoverage"]["memory"],
        cpus=1
    shell:
        """
        bedtools coverage \
            -a {input.mt} \
            -b {input.bed} > {output} 2> {log}
        """

# Collect GC Bias metric
rule GCbias:
    input: 
        bam = finalBAM,
        reference = config["reference_genome"]
    output:
        metric: "output/{sample_id}/{cell_id}/0_qc/gcBias/{cell_id}.gcBias_metric.txt"
        chart: "output/{sample_id}/{cell_id}/0_qc/gcBias/{cell_id}.gcBias_metric.pdf"
        summary: "output/{sample_id}/{cell_id}/0_qc/gcBias/{cell_id}.gcBias_summary.txt"\
    singularity: "docker://quay.io/biocontainers/picard-slim_2.27.4--hdfd78af_0.sif"
    log: "output/{sample_id}/{cell_id}/log/GCbias.log"
    threads: 1
    resources:
        mem_mb=config["GCbias"]["memory"],
        cpus=1
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/gcBias

        picard CollectGcBiasMetrics \
            -I {input.bam} \
            -O {output.metrics} \
            --CHART_OUTPUT {output.chart} \
            -S {output.summary} \
            -R {input.reference} \
            --VALIDATION_STRINGENCY LENIENT &> {log}
        """

rule multiQC:
    input: 
        gcBias = "output/{sample_id}/{cell_id}/0_qc/gcBias/{cell_id}.gcBias_metric.txt",
        mtCov = "output/{sample_id}/{cell_id}/0_qc/mtDNA/{cell_id}.mtCov.bed",
        WGSmetric = "output/{sample_id}/{cell_id}/0_qc/WgsMetrics/{cell_id}.wgsMetric.txt",
        insertSize = "output/{sample_id}/{cell_id}/0_qc/insertsize/{cell_id}.insertSize.txt",
        flagstat = "output/{sample_id}/{cell_id}/0_qc/flagstat/{cell_id}.flagstat.txt",
        markDup = "output/{sample_id}/{cell_id}/0_qc/markDup/{cell_id}.markDup_metric.txt",
        fastq_screen1 = "output/{sample_id}/{cell_id}/0_qc/fastq_screen/{cell_id}_1_screen.txt",
        fastq_screen2 = "output/{sample_id}/{cell_id}/0_qc/fastq_screen/{cell_id}_2_screen.txt",
        zip_1 = "output/{sample_id}/{cell_id}/0_qc/fastqc/{cell_id}_1_fastqc.zip",
        zip_2 = "output/{sample_id}/{cell_id}/0_qc/fastqc/{cell_id}_2_fastqc.zip"
    output: "output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/multiqc_report.html"
    singularity: "docker://quay.io/biocontainers/multiqc_1.19--pyhdfd78af_0"
    log: "output/{sample_id}/{cell_id}/log/multiQC.log"
    threads: 1
    resources:
        mem_mb=config["multiQC"]["memory"],
        cpus=1
    shell:
        """
        multiqc \
            -c src/multiqc.yaml \
            --force \
            --outdir output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc \
            output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc &> {log}
        """