#!/usr/bin/env python3
import subprocess
import re
import argparse
import sys
import os
import shutil  # Added for safe work directory cleanup
from multiprocessing import Pool
from functools import partial
from Bio import SeqIO

def process_chromosome(chrom, arguments, base_dir):
    work_dir = os.path.join(base_dir, "work_dirs", chrom)
    os.makedirs(work_dir, exist_ok=True)

    new_outdir = os.path.join(base_dir, arguments.output_directory, chrom)
    sequence_id = arguments.sequence_id if arguments.sequence_id else chrom

    MAIN_NF = os.path.join(base_dir, "main2.nf")
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
        return (chrom, True) # Return success status
    except subprocess.CalledProcessError:
        print(f"FAILURE: {chrom} encountered an error during Nextflow cluster discovery.")
        return (chrom, False) # Return failure status

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
    
    base_dir = os.getcwd()

    # --- ADDITION 1: Create Results and Output Parent Directories Natively ---
    print("\nPreparing target output filesystem paths...")
    target_output_root = os.path.join(base_dir, args.output_directory)
    os.makedirs(target_output_root, exist_ok=True)

    # --- ADDITION 2: Instantiate and Initialize SQLite Schema Globally ---
    print(f"Initializing central pipeline database sequence: {args.db_path}")
    # Ensure directory containing the database actually exists
    db_dir = os.path.dirname(os.path.abspath(args.db_path))
    if db_dir:
        os.makedirs(db_dir, exist_ok=True)
        
    db_init_cmd = ["python", "./scripts/00_db_prep.py", os.path.abspath(args.db_path)]
    try:
        subprocess.run(db_init_cmd, check=True)
        print("Database initialized successfully.")
    except subprocess.CalledProcessError as e:
        print(f"CRITICAL ERROR: Could not create base SQLite database schema ({e}). Exiting.")
        sys.exit(1)

    print(f"\nStarting parallel cluster generation using {args.workers} workers...")
    worker_func = partial(process_chromosome, arguments=args, base_dir=base_dir)

    # 1. RUN PARALLEL NEXTFLOW JOBS
    with Pool(processes=args.workers) as pool:
        # Collect execution results to verify if any individual chromosome crashed
        execution_results = pool.map(worker_func, chr_list)
        
    print("\nNextflow jobs finished for all chromosomes. Transitioning to global filtering steps...")

    # 2. RUN STEP 6: ADDITIONAL DATABASE FILTERING (Executes once globally)
    print("\n=========================================================")
    print("RUNNING GLOBAL STEP 6: Additional Database Filtering")
    print("=========================================================")
    
    step6_cmd = [
        "python", "scripts/Additional_filtering.py",
        "--db", args.db_path,
        "--outdir", target_output_root,
    ]
    
    step6_success = False
    try:
        subprocess.run(step6_cmd, check=True)
        print("Step 6 execution complete.")
        step6_success = True
    except subprocess.CalledProcessError as e:
        print(f"Error executing Additional_filtering.py: {e}")
        sys.exit(1)

    # --- ADDITION 3: Post-Pipeline Verification and Selective Work Dir Erasure ---
    # Check if Step 6 succeeded, and make sure no Nextflow worker failed
    nextflow_failures = [chrom for chrom, success in execution_results if not success]
    
    if step6_success and len(nextflow_failures) == 0:
        print("\n=========================================================")
        print("SUCCESS AUDIT PASSED: Clearing internal temporary work data.")
        print("=========================================================")
        for chrom in chr_list:
            # Targets the local runtime work cache directory structure
            target_work_path = os.path.join(base_dir, "work_dirs", chrom, "work")
            if os.path.exists(target_work_path):
                try:
                    shutil.rmtree(target_work_path)
                    print(f"Cleaned up caching directory for: {chrom}")
                except Exception as cleanup_error:
                    print(f"Warning: Could not remove directory {target_work_path}: {cleanup_error}")
        print("Server directory clean-up successfully finalized.")
    else:
        print("\n=========================================================")
        print("WARNING: Pipeline optimization cleanup skipped due to execution errors.")
        if nextflow_failures:
            print(f"The following chromosomes encountered failures: {', '.join(nextflow_failures)}")
        if not step6_success:
            print("Step 6 Post-processing engine failed to resolve successfully.")
        print("=========================================================")


if __name__ == '__main__':
    main()