## Load config values
configfile: "config/input.yaml"
configfile: "config/parameter.yaml"
configfile: "config/resource.yaml"

# Get ids
input_dir = config["input_dir"]
sample_id = config["sample_id"]
cell_ids = config['cell_ids']

rule all:
    input: expand("output/"+sample_id+"/{cell_id}/2_alignment/{cell_id}.markDup.bam", cell_id = cell_ids)

###############
# Preprocess fastq
###############

# Trims adapter sequence and poly-G tail
rule fastp:
    input:
        fq1 = input_dir + "/{cell_id}_1.fastq.gz",
        fq2 = input_dir + "/{cell_id}_2.fastq.gz"
    output:
        fq1 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_1_trimmed.fastq.gz",
        fq2 = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_2_trimmed.fastq.gz",
        html = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_fastp.html",
        json = "output/{sample_id}/{cell_id}/1_preprocessing/{cell_id}_fastp.json"
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
            TODO
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

# Mark duplicate reads in BAM files
rule mark_duplication:
    input: "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.sorted.bam"
    output: 
        bam = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.markDup.bam",
        metric = "output/{sample_id}/{cell_id}/2_alignment/{cell_id}.markDup_metric.txt"
    singularity: "docker://quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
    log: "output/{sample_id}/{cell_id}/log/sort_index_bam.log"
    threads: 1
    resources:
        mem_mb = config["mark_duplication"]["memory"],
        cpus=1
    shell:
        """
        gatk MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metric} &> {log}
        """

###############
# QC Metrics
###############

# Generate fastq screen config file
rule fastq_screen_config:
    output: "output/fastq.config"
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
    singularity: "docker://quay.io/biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0"
    threads: 20
    log: "output/{sample_id}/{cell_id}/log/fastq_screen.log"
    shell:
        """
        mkdir -p output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastp
        
        fastq_screen \
            --aligner bowtie \
            --conf {input.config} \
            --outdir output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastp \
            --tag \
            --thread {threads} \
            {input.fq1} {input.fq2} &> {log}
        
        mv output/{wildcards.sample_id}/{wildcards.cell_id}/0_qc/fastp/*tagged.fastq.gz \
            output/{wildcards.sample_id}/{wildcards.cell_id}/1_preprocessing
        """