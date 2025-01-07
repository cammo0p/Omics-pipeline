# rules/try_preprocessing.smk

# FastQC rule for raw R1 and R2 files
rule raw_fastqc:
    input:
        # Input is the raw FASTQ file located in rawdata (both R1 and R2)
        fastq_R1 = f"{base_directory}/rawdata/{{id}}_R1_001.fastq.gz",
        fastq_R2 = f"{base_directory}/rawdata/{{id}}_R2_001.fastq.gz"
    output:
        # Output for FastQC: HTML reports for both R1 and R2
        f"{current_directory}/output/fastQC/raw/{id}_R1_001_fastqc.html",
        f"{current_directory}/output/fastQC/raw/{id}_R2_001_fastqc.html"
    threads: 8
    params:
        outdir = f"{current_directory}/output/fastQC/raw"
    shell:
        """
        fastqc {input.fastq_R1} -o {params.outdir}
        fastqc {input.fastq_R2} -o {params.outdir}
        """

# Trimmomatic rule for trimming raw reads
rule trimmomatic_trim:
    input:
        R1 = f"{base_directory}/rawdata/{{id}}_R1_001.fastq.gz",
        R2 = f"{base_directory}/rawdata/{{id}}_R2_001.fastq.gz"
    output:
        # Trimmomatic produces trimmed paired and unpaired files
        R1_paired = "output/fastq/trimmed/{id}_tr_1P.fq.gz",
        R2_paired = "output/fastq/trimmed/{id}_tr_2P.fq.gz",
        R1_unpaired = "output/fastq/trimmed/{id}_tr_1U.fq.gz",
        R2_unpaired = "output/fastq/trimmed/{id}_tr_2U.fq.gz"
    threads: 32  
    shell:
        """
        trimmomatic PE -threads {threads} \
        {input.R1} {input.R2} \
        {output.R1_paired} {output.R1_unpaired} \
        {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:{config[adapters]} MINLEN:{config[minlen]}
        """

# Merging trimmed FASTQ files from L002 and L003 lanes into final files
rule cat_trimmed_fastq:
    input:
        L002_tr_1P = "output/fastq/trimmed/{id_base}_L002_tr_1P.fq.gz",
        L003_tr_1P = "output/fastq/trimmed/{id_base}_L003_tr_1P.fq.gz",
        L002_tr_2P = "output/fastq/trimmed/{id_base}_L002_tr_2P.fq.gz",
        L003_tr_2P = "output/fastq/trimmed/{id_base}_L003_tr_2P.fq.gz"
    output:
        R1_merged = "output/fastq/final/{id_base}_R1.fq.gz",
        R2_merged = "output/fastq/final/{id_base}_R2.fq.gz"
    shell:
        """
        cat {input.L002_tr_1P} {input.L003_tr_1P} > {output.R1_merged}
        cat {input.L002_tr_2P} {input.L003_tr_2P} > {output.R2_merged}
        """

# FastQC rule for merged R1 and R2 files
rule merged_fastqc:
    input:
        # Input is the raw FASTQ file located in rawdata (both R1 and R2)
        fastq_R1_merged = "output/fastq/final/{id_base}_R1.fq.gz",
        fastq_R2_merged = "output/fastq/final/{id_base}_R2.fq.gz"
    output:
        # Output for FastQC: HTML reports for both R1 and R2
        "output/fastQC/final/{id_base}_R1_fastqc.html",
        "output/fastQC/final/{id_base}_R2_fastqc.html"
    threads: 8
    params:
        outdir = "output/fastQC/final"
    shell:
        """
        fastqc {input.fastq_R1_merged} -o {params.outdir}
        fastqc {input.fastq_R2_merged} -o {params.outdir}
        """
