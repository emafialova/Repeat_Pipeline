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

    c.execute('''CREATE TABLE IF NOT EXISTS Hashes (
        sequence TEXT , 
        short_hash TEXT, 
        full_hash TEXT PRIMARY KEY)''')
    
    c.execute('''CREATE TABLE IF NOT EXISTS Instances (
        internal_id TEXT PRIMARY KEY, 
        chrom TEXT, 
        start_pos INTEGER, 
        end_pos INTEGER, 
        human_id TEXT,
        full_hash TEXT,
        FOREIGN KEY(full_hash) REFERENCES Hashes(full_hash)
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS Cluster_details (
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
        med_occ REAL,
        FOREIGN KEY(internal_id) REFERENCES Instances(internal_id)
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS Core_kmer_info (
        internal_id TEXT, 
        core_seq TEXT, 
        positions TEXT, 
        PRIMARY KEY (internal_id, core_seq),
        FOREIGN KEY(internal_id) REFERENCES Instances(internal_id)
    )''')

    c.execute('''CREATE TABLE IF NOT EXISTS Final_kmer_info (
        internal_id TEXT, 
        kmer_seq TEXT, 
        positions TEXT, 
        is_rep INTEGER,
        core_kmer TEXT,
        PRIMARY KEY (internal_id, kmer_seq),
        FOREIGN KEY(internal_id) REFERENCES Instances(internal_id),
        FOREIGN KEY(internal_id, core_kmer) REFERENCES Core_kmer_info(internal_id, core_seq)
    )''')
 
    # Essential Indexes 
    c.execute("CREATE INDEX IF NOT EXISTS idx_human_id ON Instances(human_id)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_full_hash ON Instances(full_hash)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_kmer_lookup ON Final_kmer_info(internal_id)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_hash_sequence ON Hashes(sequence)")
    c.execute("CREATE INDEX IF NOT EXISTS idx_core_lookup ON Core_kmer_info(internal_id)")

    conn.commit()
    conn.close()
    print("Database prepared and extended successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Initialize the genomic cluster database with extended schema.")
    parser.add_argument("db_path", default="genome_clusters.db", help="Path to the SQLite database file to create or initialize.")
    args = parser.parse_args()
    initialize_global_database(args.db_path)