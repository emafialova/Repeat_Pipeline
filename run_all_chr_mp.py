import subprocess
import re
import argparse
import sys
import os
from multiprocessing import Pool
from functools import partial
from Bio import SeqIO

# This function runs on the worker threads
def process_chromosome(chrom, template_text, base_dir):
    work_dir = os.path.join(base_dir, "work_dirs", chrom)
    os.makedirs(work_dir, exist_ok=True)

    pattern_sequence_id = r'sequence_id\s*=\s*".*"'
    pattern_chromosome_id = r'chromosome_id\s*=\s*".*"'
    pattern_outdir = r'outdir\s*=\s*".*"'

    print(f"Starting analysis for: {chrom}")

    new_outdir = os.path.join(base_dir, "Results", "A1", chrom)
    # Apply substitutions
    new_config = template_text
    new_config = re.sub(pattern_sequence_id,  f'sequence_id = "{chrom}"',       new_config)
    new_config = re.sub(pattern_chromosome_id, f'chromosome_id = "{chrom}"',    new_config)
    new_config = re.sub(pattern_outdir,        f'outdir = "{new_outdir}"',      new_config)

    # Write run-specific config file
    # We use the chromosome name in the filename so workers don't overwrite each other
    new_config = new_config.replace("${projectDir}", base_dir)
    run_config_name = f"run_config_{chrom}.config"
    config_path = os.path.join(work_dir, run_config_name)
    with open(config_path, "w") as f:
        f.write(new_config)

    MAIN_NF = os.path.join(base_dir, "main.nf")

    # Build Nextflow command
    # Added -resume is safer for long parallel runs
    cmd = ["nextflow", "run", MAIN_NF, "-c", run_config_name, "-resume"]

    try:
        # check=True raises an error if the pipeline fails
        subprocess.run(cmd, check=True, cwd=work_dir)
        print(f"SUCCESS: {chrom}")
    except subprocess.CalledProcessError:
        print(f"FAILURE: {chrom} encountered an error.")

def main():
    parser = argparse.ArgumentParser(description="Run Nextflow pipeline for each chromosome")
    parser.add_argument("fasta_file", help="Path to the input FASTA file containing genomic sequences.")
    parser.add_argument("chromosome_prefix", help="Prefix of the chromosome identifiers to search for in the fasta file.")
    parser.add_argument("--chr_list", help="Optional: Comma-separated list of chromosome IDs to process (overrides fasta parsing).", default=None)
    parser.add_argument("--workers", help="Number of parallel processes", type=int, default=30)

    if len(sys.argv) == 1:
            print("\nHi! You can use this script to run Nextflow pipeline for each chromosome found in a FASTA file.")
            print("You need to provide input files and parameters.\n")
            parser.print_help(sys.stderr)
            sys.exit(1)

    args = parser.parse_args()
    
    chr_list = []

    if args.chr_list:
        chr_list = args.chr_list.split(",")
    else:
        # Extract chromosome names from the FASTA file
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            if record.id.startswith(args.chromosome_prefix):
                chr_list.append(record.id)

    print(f"Found {len(chr_list)} chromosomes matching prefix '{args.chromosome_prefix}': {chr_list}")
    print(f"Starting processing with {args.workers} parallel workers...")

    # Path to your base config template
    TEMPLATE_CONFIG = "nextflow.config"

    # Read the template config file once
    with open(TEMPLATE_CONFIG, "r") as f:
        template_text = f.read()

    base_dir = os.getcwd()
    worker_func = partial(process_chromosome, template_text=template_text, base_dir = base_dir)

    # Many workers will pick items from chr_list and run worker_func(item)
    with Pool(processes=args.workers) as pool:
        pool.map(worker_func, chr_list)
        
    print("\nAll jobs finished.")

if __name__ == '__main__':
    main()
