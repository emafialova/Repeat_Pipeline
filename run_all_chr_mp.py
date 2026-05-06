import subprocess
import re
import argparse
import sys
import os
from multiprocessing import Pool
from functools import partial
from Bio import SeqIO

# This function runs on the worker threads
def process_chromosome(chrom, arguments, base_dir):
    work_dir = os.path.join(base_dir, "work_dirs", chrom)
    os.makedirs(work_dir, exist_ok=True)

    new_outdir = os.path.join(base_dir, arguments.output_directory, chrom)
    sequence_id = arguments.sequence_id if arguments.sequence_id else chrom

    MAIN_NF = os.path.join(base_dir, "main.nf")
    print(f"Starting analysis for: {chrom}")

    sequence_path = os.path.join(base_dir, arguments.fasta_file)
    db_path = os.path.join(base_dir,arguments.db_path)
    # Build Nextflow command
    # Added -resume is safer for long parallel runs
    cmd = [
        "nextflow", "run", MAIN_NF, 
        "-resume",
        "--sequence", sequence_path,
        "--outdir", new_outdir,
        "--db_path", db_path,
        "--org", arguments.org,
        "--sequence_id", sequence_id,
        "--chromosome_id", chrom
        ]
    if arguments.mapping:
        cmd.extend(["--mapping", os.path.abspath(arguments.mapping)])
    if arguments.start is not None:
        cmd.extend(["--start", str(arguments.start)])
    if arguments.end is not None:
        cmd.extend(["--end", str(arguments.end)])

    try:
        # check=True raises an error if the pipeline fails
        subprocess.run(cmd, check=True, cwd=work_dir)
        print(f"SUCCESS: {chrom}")
    except subprocess.CalledProcessError:
        print(f"FAILURE: {chrom} encountered an error.")

def main():
    parser = argparse.ArgumentParser(description="Run Nextflow pipeline for each chromosome")
    parser.add_argument("fasta_file", help="Path to the input FASTA file containing genomic sequences.")
    parser.add_argument("chromosome_prefix", help="Prefix of the chromosome identifiers to search for in the fasta file (e.g. NC_).")
    # Required Arguments (to change it for given organism)
    parser.add_argument("--output_directory", required=True,help="Output directory for results.")
    parser.add_argument("--db_path", required=True, help="Path to the SQLite database.")
    parser.add_argument("--org", required=True, help="Organism name (e.g., 'Gallus gallus').")
    # Optional Arguments
    parser.add_argument("--start", help="Start position for analysis (default: 0).", type=int, default=0)
    parser.add_argument("--end", help="End position for analysis (default: entire sequence).", type=int, default=None)
    parser.add_argument("--mapping", help="Path to KEGG mapping file (defaults to KEGG_mapping.txt in repo).", default=None)
    parser.add_argument("--sequence_id", help="Name of analyzed sequence. Default is chromosome ID.", default=None)
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

    base_dir = os.getcwd()
    worker_func = partial(process_chromosome, arguments=args, base_dir = base_dir)

    # Many workers will pick items from chr_list and run worker_func(item)
    with Pool(processes=args.workers) as pool:
        pool.map(worker_func, chr_list)
        
    print("\nAll jobs finished.")

if __name__ == '__main__':
    main()
