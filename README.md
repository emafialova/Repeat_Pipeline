# Repeat Identification adn Analysis Pipeline

## Description
This Nextflow-based pipeline can be used to identify and analyze highly repetitive sequence clusters. Originally designed to analyze repeats belonging to the sequence stuttering pehnomenon and optimized on the genome of *Gallus gallus*. It reports repeat cluster metadata as well as calculated statistics.

## Input
There are only two required inputs for the pipeline: a FASTA sequence file (example is in test_sequence folder) - most often a whole genome sequence and a KEGG organism mapping file (included in this repo).

## Output
There are several generated outputs:
- SQLite database populated by repeat clusters
- BED file 
It is also possible to generate an HTML visualization using the provided script: html_vis_db_usage.py

## Instalation
```bash
git clone [https://github.com/YourUsername/YourRepoName.git](https://github.com/YourUsername/YourRepoName.git)
cd YourRepoName

# Create and activate environment
conda create -n stutter_env python=3.9
conda activate stutter_env
```


## Prerequisites

