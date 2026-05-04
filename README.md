# Repeat Identification and Analysis Pipeline

## Description
This Nextflow-based pipeline can be used to identify and analyze highly repetitive sequence clusters. Originally designed to analyze repeats belonging to the sequence stuttering pehnomenon and optimized on the genome of *Gallus gallus*. It reports repeat cluster metadata as well as calculated statistics.

## Features


## Input
There are only two required inputs for the pipeline: a FASTA sequence file (example is in test_sequence folder) - most often a whole genome sequence and a KEGG organism mapping file (included in this repo).

## Output
There are several generated outputs:
- SQLite database populated by repeat clusters
- BED file 
It is also possible to generate an HTML visualization using the provided script: html_vis_db_usage.py

## Instalation

### Prerequisites
Before running the pipeline, ensure the following core tools are installed on your system:
* **[Conda](https://docs.conda.io/en/latest/miniconda.html)** (Miniconda or Anaconda) to manage dependencies.
* **[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)** (version 22.0+) to execute the workflow.

### Environment Set Up
```bash
git clone https://github.com/emafialova/Repeat_Pipeline.git
cd YourRepoName

# Create and activate environment
conda env create -f environment.yml
conda activate pipeline_env
```


