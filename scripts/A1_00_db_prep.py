import sqlite3
import os
import argparse

def initialize_global_database(db_path):
    if os.path.exists(db_path):
        print(f"Database {db_path} already exists. Skipping initialization.")
        return

    print(f"Initializing extended genomic database: {db_path}")
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Enable WAL for parallel safety
    c.execute('PRAGMA journal_mode=WAL;')

    # --- TABLE 1: INSTANCES ---
    c.execute('''CREATE TABLE IF NOT EXISTS instances (
        internal_id TEXT PRIMARY KEY, 
        chrom TEXT, 
        start_pos INTEGER, 
        end_pos INTEGER, 
        human_id TEXT,
        short_hash TEXT,
        full_hash TEXT
    )''')

    # --- TABLE 2: CLUSTER DETAILS ---
    c.execute('''CREATE TABLE IF NOT EXISTS cluster_details (
        internal_id TEXT PRIMARY KEY, 
        rep_sequence TEXT, 
        gc REAL, 
        cons REAL, 
        med_len REAL, 
        avg_raw_lev REAL, 
        avg_norm_lev REAL, 
        cov REAL, 
        region_size INTEGER, 
        kmer_num INTEGER,
        avg_occ REAL,
        med_occ REAL
    )''')
    
    # --- TABLE 3: KMERS ---
    c.execute('''CREATE TABLE IF NOT EXISTS kmers (
        internal_id TEXT PRIMARY KEY, 
        kmer_list TEXT
    )''')

    # --- TABLE 4: FINAL KMER INFO ---
    c.execute('''CREATE TABLE IF NOT EXISTS final_kmer_info (
        internal_id TEXT, 
        kmer_seq TEXT, 
        positions TEXT, 
        is_rep INTEGER,
        core_kmer TEXT,
        FOREIGN KEY(internal_id) REFERENCES instances(internal_id)
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS hashes (
        sequence TEXT PRIMARY KEY, 
        short_hash TEXT, 
        full_hash TEXT)''')
    
    c.execute('''CREATE TABLE IF NOT EXISTS core_kmer_info (
        internal_id TEXT, 
        core_seq TEXT, 
        positions TEXT, 
        FOREIGN KEY(internal_id) REFERENCES instances(internal_id)
    )''')
    
    # Essential Indexes 
    c.execute("CREATE INDEX IF NOT EXISTS idx_human_id ON instances(human_id)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_full_hash ON instances(full_hash)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_kmer_lookup ON final_kmer_info(internal_id)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_hash_sequence ON hashes(sequence)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_core_lookup ON core_kmer_info(internal_id)")

    conn.commit()
    conn.close()
    print("✅ Database prepared and extended successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Initialize the genomic cluster database with extended schema.")
    parser.add_argument("db_path", default="genome_clusters.db", help="Path to the SQLite database file to create or initialize.")
    args = parser.parse_args()
    initialize_global_database(args.db_path)