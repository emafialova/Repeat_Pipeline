# Repeat Identification and Analysis Pipeline

## Description
This Nextflow-based pipeline can be used to identify and analyze highly repetitive sequence clusters. Originally designed to analyze repeats belonging to the sequence stuttering phenomenon and optimized on the genome of *Gallus gallus*. It reports repeat cluster metadata as well as calculated statistics.

## Features
- localizes imperfect repeats
- thresholds adaptable to a given usecase
- provides output in the form of an SQLite database 
- calculates statistics for each cluster (see table with details below)
- possible HTML visualization for selected repeat cluster `html_vis_db_usage.py`

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
## Usage

### Paralel Run for Whole Genome Analysis
For processing optimization, the pipeline can be executed across multiple chromosomes in paralel. This is the recommended approach for whole-genome analysis and is managed by the provided python wrapper script: `run_all_chr_mp.py`
To run the parallel execution, use the following command:
```bash
python run_all_chr_mp.py --data_source path/to/fasta --chromosome_prefix chr --workers 10
```
If you wish to use the pipeline only on a selected subset of the chromosomes, use the optional chr_list flag:
```bash
python run_all_chr_mp.py --data_source path/to/fasta --chromosome_prefix chr --workers 10 --chr_list "chr1,chr2,chr3"
```

### Single nextlfow pipeline run
If you want to process a single sequence, you can bypass the wrapper and run the Nextflow pipeline directly. The pipeline's default parameters are set in the `nextflow.config` file, however, you can override any of the parameters using the -- flag.
```bash
nextflow run main.nf \
    --sequence path/to/your_sequence.fna \
    --chromosome_id "chr1" \
    --sequence_id "target_gene_name" \
    --outdir Results/custom_run
    --db_path path/to/db
```

## Input
There are only two required inputs for the pipeline: a FASTA sequence file (example is in test_sequence folder) - most often a whole genome sequence and a KEGG organism mapping file (`KEGG_mapping.txt` included in this repo).

## Output
There are several generated outputs:
- SQLite database populated by repeat clusters
- BED file 
It is also possible to generate an HTML visualization using the provided script: `html_vis_db_usage.py`




