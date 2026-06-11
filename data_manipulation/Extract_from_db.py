#!/usr/bin/env python3
import argparse
import sqlite3
import sys
import pandas as pd

# Define the definitive mapping of user-friendly names to actual DB table columns
COLUMN_MAPPING = {
    # Instances table
    'id': ('instances.human_id', 'instances'),
    'internal_id': ('instances.internal_id', 'instances'),
    'chrom': ('instances.chrom', 'instances'),
    'start': ('instances.start_pos', 'instances'),
    'stop': ('instances.end_pos', 'instances'),
    
    # Cluster Details table
    'rep_sequence': ('cluster_details.rep_sequence', 'cluster_details'),
    'gc': ('cluster_details.gc', 'cluster_details'),
    'cons': ('cluster_details.cons', 'cluster_details'),
    'med_len': ('cluster_details.med_len', 'cluster_details'),
    'kmer_num': ('cluster_details.kmer_num', 'cluster_details'),
    'region_size': ('cluster_details.region_size', 'cluster_details'),
    'cov': ('cluster_details.cov', 'cluster_details'),
    'avg_occ': ('cluster_details.avg_occ', 'cluster_details'),
    'med_occ': ('cluster_details.med_occ', 'cluster_details'),
    'avg_raw_lev': ('cluster_details.avg_raw_lev', 'cluster_details'),
    'avg_norm_lev': ('cluster_details.avg_norm_lev', 'cluster_details'),
    
    # Hashes table
    'sequence': ('hashes.sequence', 'hashes'),
    'short_hash': ('hashes.short_hash', 'hashes'),
    'full_hash': ('hashes.full_hash', 'hashes'),
    
    # Core Kmer Info table
    'core_seq': ('core_kmer_info.core_seq', 'core_kmer_info'),
    'core_positions': ('core_kmer_info.positions', 'core_kmer_info'),
    
    # Final Kmer Info table
    'kmer_seq': ('final_kmer_info.kmer_seq', 'final_kmer_info'),
    'final_positions': ('final_kmer_info.positions', 'final_kmer_info'),
    'is_rep': ('final_kmer_info.is_rep', 'final_kmer_info'),
    'core_kmer': ('final_kmer_info.core_kmer', 'final_kmer_info')
}

def build_dynamic_query(requested_columns):
    """Determines required tables and dynamically compiles an SQL JOIN string."""
    select_clauses = []
    required_tables = set()
    
    for col in requested_columns:
        if col not in COLUMN_MAPPING:
            print(f"Error: Column '{col}' is invalid.")
            print(f"Available options are: {', '.join(COLUMN_MAPPING.keys())}")
            sys.exit(1)
        
        db_col, table_name = COLUMN_MAPPING[col]
        # Use 'AS' to ensure the output CSV column headers match the user's requested names
        select_clauses.append(f"{db_col} AS {col}")
        required_tables.add(table_name)
        
    # 'instances' is our core dimensional entity table, so it serves as the base
    query = f"SELECT {', '.join(select_clauses)} FROM instances"
    
    # Left join tables conditionally only if fields from them were requested
    if 'cluster_details' in required_tables:
        query += " LEFT JOIN cluster_details ON instances.internal_id = cluster_details.internal_id"
        
    if 'hashes' in required_tables:
        query += " LEFT JOIN hashes ON instances.full_hash = hashes.full_hash"
        
    if 'core_kmer_info' in required_tables:
        query += " LEFT JOIN core_kmer_info ON instances.internal_id = core_kmer_info.internal_id"
        
    if 'final_kmer_info' in required_tables:
        query += " LEFT JOIN final_kmer_info ON instances.internal_id = final_kmer_info.internal_id"
        
    query += ";"
    return query

def main():
    parser = argparse.ArgumentParser(description="Custom SQLite pipeline data extractor to CSV.")
    parser.add_argument("-d", "--database", required=True, help="Path to your SQLite database file")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output CSV file")
    parser.add_argument("-c", "--columns", required=True, nargs='+', 
                        help="Space-separated list of desired columns (e.g., chrom id start stop gc core_seq region_size)")
    
    args = parser.parse_args()
    
    # 1. Generate the customized SQL statement based on command-line arguments
    sql_query = build_dynamic_query(args.columns)
    
    print(f"Connecting to database: {args.database}")
    print(f"Executing dynamic query to compile requested fields: {', '.join(args.columns)}")
    
    try:
        # 2. Execute query via standard library and hand off table matrix to Pandas
        conn = sqlite3.connect(args.database)
        df = pd.read_sql_query(sql_query, conn)
        conn.close()
        
        # 3. Save report output
        df.to_csv(args.output, index=False)
        print(f"Success! Extracted {len(df)} rows. Saved to: {args.output}\n")
        
    except Exception as e:
        print(f"Database extraction failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()