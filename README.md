# Repeat Identification and Analysis Pipeline

## Description
This Nextflow-based pipeline can be used to identify and analyze highly repetitive sequence clusters. Originally designed to analyze repeats belonging to the sequence stuttering phenomenon and optimized on the genome of *Gallus gallus*. It reports repeat cluster metadata as well as calculated statistics.

Avian genomes exhibit distinct characteristics, shaped by extensive adaptations during their evolution within the dinosaur lineage. For a long time, a number of genes were believed to be evolutionarily lost in avian genomes. However, recent studies suggest that many of these genes are not truly absent but rather located in regions that are technically difficult to analyze—such as microchromosomes. These regions are characterized by high GC content and a high density of repetitive elements, including the sequence stuttering phenomenon. These characteristics may contribute to genomic instability and their further analysis can provide insights into mechanisms of evolutionary change and selection. To address this need, a specialized computational pipeline capable of the identification and analysis of repetitive features, including the sequence stuttering phenomenon, has been developed. 

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Parameters](#parameters)
- [Example Usage](#example-usage)
- [Citation](#citation)
- [Contact](#contact)

## Features
- localizes imperfect repeats
- thresholds adaptable to a given usecase
- provides output in the form of an SQLite database 
- calculates statistics for each cluster (see table with details below)
- possible HTML visualization for selected repeat cluster `html_vis_db_usage.py`

## Installation

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

### Database Set Up
Since the information generated from the pipeline is being stored in an SQLite3 database, it must be initialized before the pipeline execution using the following command:
```bash
python A1_00_db_prep.py path/to/db
```

### Parallel Run for Whole Genome Analysis
For processing optimization, the pipeline can be executed across multiple chromosomes in parallel. This is the recommended approach for whole-genome analysis and is managed by the provided Python wrapper script: `run_all_chr_mp.py`
To run the parallel execution, use the following command:
```bash
python run_all_chr_mp.py --data_source path/to/fasta --chromosome_prefix chr --workers 10
```
If you wish to use the pipeline only on a selected subset of the chromosomes, use the optional chr_list flag:
```bash
python run_all_chr_mp.py --data_source path/to/fasta --chromosome_prefix chr --workers 10 --chr_list "chr1,chr2,chr3"
```

### Single Nextflow pipeline run
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
- BED file 
- SQLite database populated by repeat clusters
Below, you can find an Entity-Relationship DIagram for the SQLite database:

![DB](Images/DB_setup.png) 

It is also possible to generate an HTML visualization using the provided script: `html_vis_db_usage.py`

## Final Parameter configuration
This pipeline consists of 5 steps, each step has several parameters which can be altered by the user. The results shown below were generated with the following parameter configuration:

| Step | Parameter | Value |
| :--- | :--- | :--- |
| Step 1 |`Min k-mer Length (k)` | 10 | 
| Step 1 |`Min k-mer Frequency (m)` | 10 | 
| Step 1 |`Sliding Window Size (W)` | 10 000 | 
| Step 1 |`Sliding Window Step Size (D)` | 5 000 | 
| Step 2 |`Max Distance Between Occurrences (d)` | 1 000 | 
| Step 2 |`Min Region Length (l)` | 30 | 
| Step 2 |`Min Number of k-mer in a Region (n)` | 10 | 
| Step 3 |`Min Number of k-mer Occurrences (r)` | 7 | 
| Step 3 |`Max k-mer Length (K)` | 2 000 | 
| Step 3 |`Max Hamming Distance (h)` | 15% | 
| Step 3 |`Frequency Threshold (f)` | GC(region) | 
| Step 3 |`Time Limit (t)` | 300 s | 
| Step 5 |`DBSCAN (eps)` | 0.45 | 
| Step 5 |`DBSCAN (MinPts)` | 2 | 
| Additional Filtering |`Region Length` | 70 | 
| Additional Filtering |`Occurrences Median` | 7 | 
| Additional Filtering |`Extension` | YES | 

The configuration was set up on **Gallus gallus genome - assembly GCA_024206055.2_GGswu**.

## Example usage
I will demonstrate the usage of the pipeline on the gene AKT2 which is located on the 32nd chromosome **(Gallus gallus genome, assembly GCA_024206055.2_GGswu: CP100586.2:2596287-2608733)**.
Below is a self dotplot generated using the YASS program [1] with default parameters. It is visible that there are several repeat clusters, identifying and analyzing them is the objective of the pipeline.

![Dotplot AKT2](Images/Dotplot_AKT2.png) 

The sequence of this gene is available in the folder `test_data`, results of the pipeline generated for this gene are available in the folder `test_results`. Final clusters bed file is called: `cluster_islands_C7_new.bed`, database with all clusters is available as `clusters_db_AKT2.db`.

### Command
The results were generated using the command below:
```bash
nextflow run main.nf \                        
    --sequence ./test_data/AKT2_seq.fasta \
    --chromosome_id "CP100586.2" \
    --sequence_id "AKT2" \
    --outdir ./test_results/AKT2 \
    --db_path "$(pwd)/test_results/clusters_db_AKT2.db"
```
### Output 
The information about clusters identified and analyzed by this pipeline is automatically loaded into the prepared database. It can be accessed via command line using sqlite3. The database stores information for all clusters across all chromosomes.

Moreover, for each sequence/chromosome, two output files, which store the information, are generated:
- bed file which stores information about position of a cluster and a cluster ID `cluster_islands_C7_new.bed`
- tsv file with cluster information `cluster_islands_C7_new_details.tsv`
The bed file can be loaded into standard bioinformatic tools such as UCSC Genome Browser as custom track, allowing the user to see the positions of the clusters. 

The statistics calculated for each cluster present in both database and the tsv file include:
| Statistic | Description | 
| :--- | :--- | 
| **Region Size (bp)** | The absolute length of the clustered genomic region, from the start of the first occurrence of a final k-mer to the end of the last. |
| **GC Content (%)** | The percentage of guanine and cytosine bases across the whole region. |
| **Coverage Ratio** | The fraction of the total region size that is covered by the final k-mers. *This demonstrates the density of the repeats within the defined region.* |
| **Number of k-mers** | The absolute count of distinct final k-mers present. *This demonstrates the variability of the k-mers present in the region.*  |
| **Median Length (bp)** | The median length of all final k-mers present in the cluster. Using median rather than mean prevents k-mers of extreme sizes of skewing the descriptor. |
| **Average and Median occurrences** | The mean and median frequency counts of all final k-mers. *These metrics show the repetitive nature of the final k-mers in the cluster.* |
| **Conservation Score** | A measure of how similar the final k-mers are to the representative sequence of the cluster. Measured using Levenshtein distance between the k-mers. |
| **Average Raw and Normalized Smith-Waterman Distance** | The mean similarity scores calculated by local pairwise alignment of all final k-mers. *These metrics show the internal cluster cohesion and quantify how closely related the final k-mer sequences are.* |

### HTML Visualization
To visualize a specific cluster, it is necessary to locate its unique ID (e.g. GGA32-BA46CF24.01) in the fourth column of the output .bed file or in the database.  
HTML visualizations of selected clusters are available in the folder `test_results/HTML`, the following command was used to generate that of cluster GGA32-BA46CF24.01:
```bash
python html_vis_db_usage.py \
    --db ./test_results/clusters_db_AKT2.db \
    --fasta ./test_data/AKT2_seq.fasta \
    --id GGA32-BA46CF24.01 \
    --org "Gallus gallus" \
    --out ./test_results/HTML/report_GGA32-BA46CF24.01.html \
    --chunk 120 
```

The generated HTML is split into five parts:
- (A) header with cluster metadata such as genomic range and organism 
- (B) sequence visualization where positions of each final k-mer are shown 
- (C) cluster metrics including information about region length, coverage and GC content
- (D) phylogenetic tree showing the relations between core k-mers
- (E) multiple sequence alignment of the final k-mers

![HTML vis - part 1 AKT2](Images/HTML_parts_A-B.png)
![HTML vis - part 2 AKT2](Images/HTML_part_C.png)

## Summary
This repository presents a specialized computational pipeline capable of performing a systematic analysis of complex repetitive regions. The pipeline provides a robust framework for repetitive cluster analysis, offering outputs in the form of a database as well as standardized bioinformatic BED file format. Furthermore, the user is also able to generate an HTML report with detailed information about the repetitive cluster. Ultimately, this pipeline serves as a key tool for the analysis of repetitive regions, providing valuable insights that were previously obscured. 

## Citation
If you use this pipeline, please cite:

Fialová, E. (2026). *Repetitive Elements and Stutter Genes in Avian Genomes*. Master's thesis, Univeristy of Chemistry and Technology, Prague, Czech Republic.

## Contact
For any questions or support, please contact:
- ema.fialova@img.cas.cz

This work was carried out with the support of ELIXIR CZ Research Infrastructure (ID LM2023055, MEYS CR)

## References
[1] Noé, L., & Kucherov, G. (2005). YASS: enhancing the sensitivity of DNA similarity search. *Nucleic Acids Research*, 33(suppl_2), W540-W543. https://doi.org/10.1093/nar/gki478
