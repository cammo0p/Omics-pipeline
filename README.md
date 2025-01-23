# Omics-pipeline

A pipeline for DNA-seq, RNA-seq, and ATAC-seq analysis. This repository is currently under development, with pipelines being built for each analysis type.

## Pipeline Structure

- **`Snakefile`**: Acts as the main blueprint, orchestrating the workflow by incorporating rules and commands from each `.smk` file in the `rules` directory.
- **`config/` directory**: Contains YAML files with paths, parameters. Also, includes a `.tsv` file listing sample IDs and a configuration file specifying paths and settings.
- **`envs/` directory**: Contains YAML files specifying dependencies and version requirements for each analysis environment.
- **`rules/` directory**: Contains modular Snakemake (`.smk`) files, each defining a specific step in the workflow.

## Main Workflow Steps for ATAC-seq

The pipeline consists of three primary stages, each requiring a corresponding `.smk` file to be executed:

1. **`preprocessing.smk`**
   - **Purpose**: Preprocesses raw sequencing data.
   - **Description**: Each sample has data split across two lanes, necessitating a merging step after trimming to produce final FASTQ files for downstream analysis.
   - **Key Steps**:
     - **[FastQC]** - Quality control on raw reads.
     - **[Trimmomatic]** - Trimming adapters and low-quality sequences.
     - Merging lanes to create consolidated, trimmed FASTQ files.
     - **(Optional)** **[MultiQC]** - Summarising all samples by aggregating the FastQC HTML outputs across samples (not included in the pipeline).

2. **`alignment.smk`**
   - **Purpose**: Aligns the merged FASTQ files to the reference genome.
   - **Description**: This rule handles the alignment of preprocessed reads to a combined human-virus genome reference, preparing the data for downstream analysis such as peak calling or expression quantification.
   - **Key Steps**:
     - **[BWA MEM]** - Aligns the reads to the reference genome using BWA-MEM.
     - **[GATK]** - (Run with Docker) Marks duplicates in the aligned BAM file to improve downstream analyses.
     
3. **`filtering.smk`**
   - **Purpose**: Removes reads that are not useful for peakcalling. 
   - **Description**:This rule filters out mitochondrial (MT) reads and duplicates, generates a flagstat report for quality control.
   - **Key Steps**:
   -  **[Samtools]** - Filter out mitochondrial (MT) reads (Harvard ATAC-seq module for removeChrom python script), Keep properly paired reads (-f 3), Name sort the BAM for fixmate, Fix mate information, Sort and index the files, Convert BAM to BEDPE format, Generate a flagstat report
   -   **[Picard]** - Remove PCR duplicates
       
4. **`peakcalling.smk`**
   - **Purpose**: Calls peaks from the BAM files, identifying regions of significant enrichment.
   - **Description**: This rule runs **MACS2** to call peaks on filtered BAM files (after removing mitochondrial reads and duplicates). The results can be used for downstream analysis such as differential peak calling or visualisation.
   - **Key Steps**:
     - **[MACS2 Call Peak]** - Uses MACS2 to call peaks from BAM files, specifying the format (`BAMPE` for paired-end), genome size, and keeping duplicate reads.
