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
        "output/fastQC/raw/fastqc.html"
    threads: 2
    params:
        # Directory where FastQC will save results
        outdir = "output/fastQC/raw"
    shell:
        """
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
        docker run -it multiqc/multiqc
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

# Merge lanes
rule cat_trimmed_fastq
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

# Fastqc run for merged FastQs
rule merged_fastqc:
    input: 
        # Input is the raw FASTQ file located in rawdata
        merged_fastq_R1 = "output/fastq/final/{{id}}_R1.fq.gz"
        merged_fastq_R2 = "output/fastq/final/{{id}}_R2.fq.gz"
    output:
        # Output for FastQC: an HTML report and a zip file
        "output/fastQC/final/fastqc.html"
    threads: 2
    params:
        # Directory where FastQC will save results
        outdir = "output/fastQC/final"
    shell:
        """
        docker run -it multiqc/multiqc \
        fastqc {input.fastq} -o {params.outdir}
        """

rule merged_multiqc
    input:
        # MultiQC will use all FastQC HTML reports generated in the previous rule
        expand("output/fastQC/final/fastqc.html")
    output:
        # Output is the consolidated MultiQC report
        multiqc_html = "output/fastQC/final/multiqc/multiqc_report.html"
    shell:
        """
       docker-multiqc .
        """

# BWA-MEM to align reads to the reference genome
rule bwa_align:
    input:
        # Input is the final, merged FASTQ files
        fastq_R1 = "output/fastq/final/{id_base}_R1.fq.gz",
        fastq_R2 = "output/fastq/final/{id_base}_R2.fq.gz"
    output:
        # Sorted BAM file and index
        bam_sorted = "output/bam/{id_base}.sorted.bam",
        bam_index = "output/bam/{id_base}.sorted.bam.bai"
    threads: 32
    shell:
        """
        bwa mem {ref_combined_genome} {input.fastq_R1} {input.fastq_R2} > output/bam/{wildcards.id_base}.sam && \
        samtools view -bS output/bam/{wildcards.id_base}.sam > output/bam/{wildcards.id_base}.bam && \
        samtools sort output/bam/{wildcards.id_base}.bam -o {output.bam_sorted} && \
        samtools index {output.bam_sorted} && \
        rm output/bam/{wildcards.id_base}.sam output/bam/{wildcards.id_base}.bam
        """
# docker-gatk to mark duplicates
rule gatk_markdup:
    input:
        # Input is the sorted BAM file
        bam_sorted = "gatk/my_data/output/bam/{id_base}.sorted.bam"
    output:
        # Mark duplicates
        bam_markdup = "gatk/my_data/output/bam/{id_base}.sorted_markdup.bam",
        bam_markdup_txt = "gatk/my_data/output/bam/{id_base}.markdup_metrics.txt"
    shell:
        """
        docker run -v {current_directory}:/gatk/my_data -it broadinstitute/gatk:4.1.3.0 \
        ./gatk MarkDuplicates -I {input.bam_sorted} -O {output.bam_markdup} -M {output.bam_markdup_txt}
        """
# Filter BAM files for peak calling
# Removes mitochondrial reads and duplicate reads

rule samtools_noMT_noDup:
    input:
        # Input is the marked duplicates BAM file
        bam_markdup = "output/bam/{id_base}.sorted_markdup.bam"
    output:
        # Output BAM file with no mitochondrial reads and no duplicates
        bam_noMT_noDup = "output/bam/final{id_base}.sorted_noMT_noDup.bam",
        # Flagstat report for quality control
        flagstat_report = "output/bam/{id_base}.sorted_noMT_noDup.txt"
    shell:
        """
        # Filter out mitochondrial (MT) reads and duplicates, and create the new BAM file
        samtools view -h -F 1024 {input.bam_markdup} | grep -v 'MT' | samtools view -b -o {output.bam_noMT_noDup} && \
        
        # Index the resulting BAM file
        samtools index {output.bam_noMT_noDup} && \
        
        # Generate the flagstat report for the filtered BAM file
        samtools flagstat {output.bam_noMT_noDup} > {output.flagstat_report} \

        """
# Peak files
rule macs2_call_peak:
    input:
        bam_file = "output/bam/final/{id_base}.sorted_noMT_noDup.bam"
    output:
        peak_xls = "output/peak/{id_base}_peaks.xls",
        peak_bed = "output/peak/{id_base}_summits.bed",
        peak_narrowPeaks = "output/peak/{id_base}_peaks.narrowPeak",
        macs2_log = "output/bed/macs2_log_{id_base}.txt"
    shell:
        """
        macs2 callpeak -t {input.bam_file} -f BAMPE -n {wildcards.id_base} -g {config[genome_size]} --keep-dup all 2 > {output.macs2_log} 
        """


# BigWig files
rule deeptools_bw:
    input:
        bam_file = "output/bam/final/{id_base}.sorted_noMT_noDup.bam"
    output:
        bw_file = "output/bigwig/{id_base}_coverage.bw"
    shell:
        """
        conda activate deeptools \
        bamCoverage -b {input.bam_file} -o {output.bw_file}
        """