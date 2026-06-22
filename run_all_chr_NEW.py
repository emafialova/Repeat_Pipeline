#!/usr/bin/env python3
import subprocess
import re
import argparse
import sys
import os
from multiprocessing import Pool
from functools import partial
from Bio import SeqIO

def process_chromosome(chrom, arguments, base_dir):
    work_dir = os.path.join(base_dir, "work_dirs", chrom)
    os.makedirs(work_dir, exist_ok=True)

    new_outdir = os.path.join(base_dir, arguments.output_directory, chrom)
    sequence_id = arguments.sequence_id if arguments.sequence_id else chrom

    MAIN_NF = os.path.join(base_dir, "main.nf")
    print(f"Starting Nextflow analysis for: {chrom}")

    sequence_path = os.path.join(base_dir, arguments.fasta_file)
    db_path = os.path.join(base_dir, arguments.db_path)
    
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
        subprocess.run(cmd, check=True, cwd=work_dir)
        print(f"SUCCESS Nextflow: {chrom}")
    except subprocess.CalledProcessError:
        print(f"FAILURE: {chrom} encountered an error during Nextflow cluster discovery.")

def main():
    parser = argparse.ArgumentParser(description="Run complete pipeline suite with post-processing optimizations.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file containing genomic sequences.")
    parser.add_argument("chromosome_prefix", help="Prefix of the chromosome identifiers to search for (e.g. NC_).")
    parser.add_argument("--output_directory", required=True, help="Output directory for results.")
    parser.add_argument("--db_path", required=True, help="Path to the SQLite database.")
    parser.add_argument("--org", required=True, help="Organism name (e.g., 'Gallus gallus').")
    
    # Custom post-processing flags
    parser.add_argument("--genome_id", default="hg38", help="UCSC Genome reference string (default: hg38)")
    parser.add_argument("--filter_threshold", type=float, default=0.80, help="Overlap tolerance threshold for TRF/Exons (default: 0.80)")
    
    parser.add_argument("--start", help="Start position for analysis.", type=int, default=0)
    parser.add_argument("--end", help="End position for analysis.", type=int, default=None)
    parser.add_argument("--mapping", help="Path to KEGG mapping file.", default=None)
    parser.add_argument("--sequence_id", help="Name of analyzed sequence.", default=None)
    parser.add_argument("--chr_list", help="Optional: Comma-separated list of chromosome IDs to process.", default=None)
    parser.add_argument("--workers", help="Number of parallel processes", type=int, default=30)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    chr_list = []

    if args.chr_list:
        chr_list = args.chr_list.split(",")
    else:
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            if record.id.startswith(args.chromosome_prefix):
                chr_list.append(record.id)

    print(f"Found {len(chr_list)} chromosomes matching prefix '{args.chromosome_prefix}'")
    print(f"Starting parallel cluster generation using {args.workers} workers...")

    base_dir = os.getcwd()
    worker_func = partial(process_chromosome, arguments=args, base_dir=base_dir)

    # 1. RUN PARALLEL NEXTFLOW JOBS
    with Pool(processes=args.workers) as pool:
        pool.map(worker_func, chr_list)
        
    print("\nNextflow jobs finished for all chromosomes. Transitioning to global filtering steps...")

    # Define targets relative to your global execution paths
    target_output_root = os.path.join(base_dir, args.output_directory)
    #global_db_path = os.path.join(base_dir, args.db_path)

    # 2. RUN STEP 6: ADDITIONAL DATABASE FILTERING (Executes once globally)
    print("\n=========================================================")
    print("RUNNING GLOBAL STEP 6: Additional Database Filtering")
    print("=========================================================")
    
    step6_cmd = [
        "python", "scripts/Additional_filtering.py",
        "--db", args.db_path,
        "--outdir", target_output_root,
        #"--bed", "FINAL_clusters.bed"
    ]
    try:
        subprocess.run(step6_cmd, check=True)
        print("Step 6 execution complete.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing Additional_filtering.py: {e}")
        sys.exit(1)

    # 3. RUN STEP 7: GENOMIC REPEAT MASKING VIA UCSC (Executes once globally)
    print("\n=========================================================")
    print("RUNNING GLOBAL STEP 7: Centromere & Satellite Filtration")
    print("=========================================================")
    
    step7_cmd = [
        "python", "scripts/6_Filtering.py",
        "-g", args.genome_id,
        "-d", target_output_root,   # Root folder where chromosome folders sit
        "-o", "final_after_filt.bed",# Name of file to save in EACH folder
        "-r", "filtered_out.csv",    # Name of audit log in EACH folder
        "-t", str(args.filter_threshold)
    ]
    try:
        subprocess.run(step7_cmd, check=True)
        print("\nAll pipeline tasks successfully complete!")
        print(f"Chromosome-specific custom tracks are available as 'final_after_filt.bed' inside your output directories.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing UCSC genomic filtering: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()