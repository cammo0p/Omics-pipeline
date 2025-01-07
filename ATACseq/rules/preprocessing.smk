# Preprocessing.smk
# Data Quality Control

# rules/preprocessing.smk

# Load the base directory path from config
base_directory = config["base_directory"]

rule raw_fastqc:
    input: 
        # Input is the raw FASTQ file located in rawdata
        fastq = f"{base_directory}/rawdata/{{id}}.fastq.gz"
    output:
        # Output for FastQC: an HTML report and a zip file
        "output/fastQC/raw/{{id}}_fastqc.html"
    threads: 2
    params:
        # Directory where FastQC will save results
        outdir = "output/fastQC/raw"
    shell:
        """
        docker run -it multiqc/multiqc \
        fastqc {input.fastq} -o {params.outdir}
        """

rule raw_multiqc:
    input:
        # MultiQC will use all FastQC HTML reports generated in the previous rule
        expand("output/fastQC/raw/{id}_fastqc.html", id=SAMPLES)
    output:
        # Output is the consolidated MultiQC report
        multiqc_html = "output/fastQC/raw/multiqc/multiqc_report.html"
    shell:
        """
       docker-multiqc .
        """
rule trimmomatic_trim
    input:
        # Paired-end input files for Trimmomatic
        r1 = f"{base_directory}/rawdata/{{id}}_R1_001.fastq.gz",
        r2 = f"{base_directory}/rawdata/{{id}}_R2_001.fastq.gz"
    output:
        # Trimmomatic produces trimmed paired and unpaired files
        tr_1P = "output/fastq/trimmed/{{id}}_tr_1P.fq.gz",
        tr_1U = "output/fastq/trimmed/{{id}}_tr_1U.fq.gz",
        tr_2P = "output/fastq/trimmed/{{id}}_tr_2P.fq.gz",
        tr_2U = "output/fastq/trimmed/{{id}}_tr_2U.fq.gz"
    threads: 32
    params:
        adapters = "/opt/miniconda/envs/aidanfoo/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa"  # Path to adapter file
    shell:
        """
       	trimmomatic PE -threads {threads} \
       	{input.r1} {input.r2} {output.tr_1P} {output.tr_1U} {output.tr_2P} {output.tr_2U}
       	ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

# + fastqc and multiqc for trimmed fastq 


rule cat_trimmed_fastqc
	input:
		L002_tr_1P = "output/fastq/trimmed/{{id}}_L002_tr_1P.fq.gz",
		L003_tr_1P = "output/fastq/trimmed/{{id}}_L003_tr_1P.fq.gz",
		L002_tr_2P = "output/fastq/trimmed/{{id}}_L002_tr_2P.fq.gz",
		L003_tr_2P = "output/fastq/trimmed/{{id}}_L003_tr_2P.fq.gz",
	output:
		R1 = "output/fastq/final/{{id}}_R1.fq.gz"
		R2 = "output/fastq/final/{{id}}_R2.fq.gz"
	shell:
		"""
		cat {input.L002_tr_1P} {input.L003_tr_1P} > {output.R1}
		cat {input.L002_tr_2P} {input.L003_tr_2P} > {output.R2}
		"""
