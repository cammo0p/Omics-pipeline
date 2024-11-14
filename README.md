# HIV-pipeline

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
   - **Description**: This rule handles the alignment of preprocessed reads to a combined human-virus genome reference, preparing the data for peak calling or expression quantification.

3. **`postprocessing.smk`**
   - **Purpose**: Finalizes data for downstream analysis (e.g., peak calling for ATAC-seq).
   - **Description**: In this stage, sorted and indexed BAM files are generated, along with any additional steps required for visualization or statistical analysis.
