#!/usr/bin/env python3
import os
import re
import argparse
import sys
import subprocess
from multiprocessing import Pool
from functools import partial
from Bio import SeqIO

def parse_fasta_metadata(fasta_path):
    """
    Scans the FASTA file to find the first genuine CHROMOSOME entry.
    Extracts the organism name and the chromosome prefix dynamically, 
    completely skipping unplaced scaffolds or contigs. Supports both
    NCBI RefSeq (brackets) and GenBank metadata styles.
    """
    organism = "Unknown Species"
    prefix = ""
    
    try:
        with open(fasta_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                desc_lower = record.description.lower()
                
                # SECURITY FILTER: Only parse records explicitly labeled as a chromosome
                if "chromosome" not in desc_lower:
                    continue  # Skip unplaced scaffolds/linkers and keep looking
                
                # --- 1. DYNAMIC CHROMOSOME PREFIX EXTRACTION ---
                # e.g., 'CP100555.1' -> 'CP' or 'NC_000001.11' -> 'NC_'
                match_pref = re.match(r'^([a-zA-Z_]+)', record.id)
                if match_pref:
                    prefix = match_pref.group(1)
                else:
                    prefix = record.id[:2] # Fallback to first two characters
                
                # --- 2. DYNAMIC ORGANISM NAME EXTRACTION ---
                desc = record.description
                
                # Scenario A: RefSeq Brackets check (e.g., [Homo sapiens])
                match_bracket = re.search(r'\[(.*?)\]', desc)
                if match_bracket:
                    organism = match_bracket.group(1)
                
                # Scenario B: GenBank Format Check (e.g., ">CP100555.1 Gallus gallus breed Huxu...")
                else:
                    # Strip out the ID at the front
                    words = desc.replace(record.id, "").strip().split()
                    if len(words) >= 2:
                        # Standard binomial nomenclature is always the first two words (Genus species)
                        organism = f"{words[0]} {words[1]}"
                
                # Break immediately once the first genuine chromosome format is successfully parsed
                break
                
    except Exception as e:
        print(f"Error parsing metadata stream for {fasta_path}: {e}")
        
    return organism, prefix

def run_single_genome(fasta_path, output_root, workers_per_genome):
    """Builds and executes the execution statement for an isolated genome profile."""
    file_name = os.path.basename(fasta_path)
    
    print(f"\n[ORCHESTRATOR] Analyzing incoming target file: {file_name}")
    org_name, chr_prefix = parse_fasta_metadata(fasta_path)
    
    if chr_prefix == "":
        print(f"[SKIP] Warning: No chromosome entries found in {file_name}. Skipping pipeline calculation.")
        return

    # --- DYNAMIC INITIALS EXTRACTION ---
    # Split "Gallus gallus" into ["Gallus", "gallus"]
    name_parts = org_name.split()
    if len(name_parts) >= 2:
        # Take 'G' from Gallus and 'g' from gallus, then capitalize to make 'GG'
        org_initials = (name_parts[0][0] + name_parts[1][0]).upper()
    else:
        # Fallback if the name is just one word
        org_initials = org_name[:2].upper()
        
    print(f" -> Detected Organism:   '{org_name}'")
    print(f" -> Generated Initials:  '{org_initials}'")
    print(f" -> Detected Prefix:     '{chr_prefix}'")
    
    # Standardize paths using your clean two-letter shortcut convention
    db_name = f"clusters_db_{org_initials}.db"  # e.g., clusters_db_GG.db
    db_path = os.path.join(output_root, db_name)
    out_dir = os.path.join(output_root, org_initials)  # e.g., ./Results/GG
    
    # Build your command string targeting pipeline_wrapper.py
    cmd = [
        "python", "pipeline_wrapper.py",
        fasta_path,
        chr_prefix,
        "--workers", str(workers_per_genome),
        "--org", org_name,
        "--genome_id", org_initials, 
        "--db_path", db_path,
        "--output_directory", out_dir
    ]
    
    print(f"[RUNNING] Processing genome pipeline command for {org_initials}...")
    try:
        subprocess.run(cmd, check=True)
        print(f"[SUCCESS] Finished full pipeline routine for genome: {org_initials}")
    except subprocess.CalledProcessError as e:
        print(f"[FAILURE] Master loop error detected while processing {org_initials}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Parallel master wrapper to stream genomes sequentially using robust biological discovery flags.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing your input .fna/.fasta files")
    parser.add_argument("-o", "--output_dir", default="./Results", help="Base directory to build databases and output directories")
    parser.add_argument("--slots", type=int, default=10, help="Number of genomes to compute simultaneously (default: 10)")
    parser.add_argument("--threads_per_job", type=int, default=10, help="Workers parameter passed internally to pipeline_wrapper.py (default: 10)")
    
    args = parser.parse_args()
    
    # Collect all FASTA targets inside directory
    extensions = ('.fna', '.fasta', '.fa')
    fasta_files = [
        os.path.abspath(os.path.join(args.input_dir, f)) 
        for f in os.listdir(args.input_dir) 
        if f.lower().endswith(extensions)
    ]
    
    print(f"=========================================================")
    print(f"MASTER ORCHESTRATOR INITIATED: Found {len(fasta_files)} genomes to map.")
    print(f"Parallel Execution Target Configuration: {args.slots} slots simultaneously")
    print(f"Compute Footprint per Slot: {args.threads_per_job} threads")
    print(f"=========================================================")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up the processing pool worker footprint
    worker_func = partial(run_single_genome, output_root=args.output_dir, workers_per_genome=args.threads_per_job)
    
    with Pool(processes=args.slots) as pool:
        pool.map(worker_func, fasta_files)
        
    print("\n=========================================================")
    print("GLOBAL BATCH PROCESSING MIGRATION ROUTINE COMPLETE.")
    print("=========================================================")

if __name__ == "__main__":
    main()