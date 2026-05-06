#!/usr/bin/env python3
"""
Cluster Database Filter
Applies strict biological thresholds to genomic clusters and outputs:
1. A valid BED file of passing clusters.
2. A CSV report of all failed clusters and the specific reasons why they failed.
"""

import sqlite3
import pandas as pd
import argparse
import os

def filter_database(db_path, out_bed, out_csv, target_chrom, outdir):
    print(f"📂 Connecting to database: {db_path}")
    conn = sqlite3.connect(db_path)
    
    # 1. Extract main instance and cluster details
    query_main = """
        SELECT i.internal_id, i.human_id, i.chrom, i.start_pos, i.end_pos,
               d.med_occ, d.region_size, d.kmer_num, d.med_len, d.rep_sequence,
               d.gc, d.cons, d.avg_raw_lev, d.avg_norm_lev, d.cov, d.avg_occ
        FROM instances i
        JOIN cluster_details d ON i.internal_id = d.internal_id
    """
    df = pd.read_sql_query(query_main, conn)
    if target_chrom:
        df = df[df['chrom'] == target_chrom].copy()
    # 2. Extract the Core Branching metric
    # Groups by internal_id and core_kmer, counts the final kmers per core, 
    # and then finds the MAX of those counts for each cluster.
    query_cores = """
        SELECT internal_id, 
               MAX(CASE WHEN LENGTH(kmer_seq) > LENGTH(core_kmer) THEN 1 ELSE 0 END) as has_extended_kmer,
               GROUP_CONCAT(kmer_seq) as kmers
        FROM final_kmer_info
        GROUP BY internal_id
    """
    core_df = pd.read_sql_query(query_cores, conn)
    conn.close()

    # 3. Merge the datasets
    df = df.merge(core_df, on='internal_id', how='left')
    
    # Fill NaNs with 0 just in case a cluster has no records in final_kmer_info
    #df['max_kmers_per_core'] = df['max_kmers_per_core'].fillna(0)
    df['has_extended_kmer'] = df['has_extended_kmer'].fillna(0)
    df['kmers'] = df['kmers'].fillna("")
    print(f"📊 Total raw clusters loaded: {len(df)}")

    # 4. Define the filtering logic
    # We will build a list of failure reasons for every row
    def check_filters(row):
        reasons = []
        
        # Rule 1: Median occurrence >= 7
        if row['med_occ'] < 7:
            reasons.append(f"med_occ < 7 ({row['med_occ']})")
            
        # Rule 3: Region length >= 70
        if row['region_size'] < 70:
            reasons.append(f"region_size < 70 ({row['region_size']})")
            
        # Rule 4: At least one core must have > 1 final kmer
        #if row['max_kmers_per_core'] <= 1:
        #    reasons.append(f"max kmers per core <= 1 ({row['max_kmers_per_core']})")
        if row['has_extended_kmer'] == 0:
            reasons.append("no final kmer extended beyond its core")
            
        return " | ".join(reasons) if reasons else "PASS"

    # Apply the logic
    df['filter_status'] = df.apply(check_filters, axis=1)

    # 5. Split into Passed and Failed DataFrames
    passed_df = df[df['filter_status'] == "PASS"].copy()
    failed_df = df[df['filter_status'] != "PASS"].copy()

    print(f"✅ Clusters passed: {len(passed_df)}")
    print(f"❌ Clusters failed: {len(failed_df)}")

    # 6. Prep BED File format for Passed Clusters
    bed_df = passed_df[['chrom', 'start_pos', 'end_pos', 'human_id']].copy()
    bed_df['score'] = 0
    bed_df['strand'] = '+'
    bed_df = bed_df[['chrom', 'start_pos', 'end_pos', 'human_id', 'score', 'strand']]
    
    cols = ['human_id', 'filter_status', 'chrom', 'start_pos', 'end_pos', 
            'region_size', 'med_occ', 'kmer_num', 'med_len', 'has_extended_kmer', 'rep_sequence']

    # 7. Export files into corresponding chromosome folders
    #for chrom in df['chrom'].unique():
    #    # Ensure the chromosome directory exists
    #    chrom_folder = os.path.join(outdir, chrom)
    #    os.makedirs(chrom_folder, exist_ok=True)
    #    
    #    # Build the specific paths for this chromosome
    #    chrom_bed_out = os.path.join(chrom_folder, os.path.basename(out_bed))
    #    chrom_csv_out = os.path.join(chrom_folder, os.path.basename(out_csv))
    #    
    #    # Filter data for the current chromosome in the loop
    #    chrom_passed = bed_df[bed_df['chrom'] == chrom]
    #    chrom_failed = failed_df[failed_df['chrom'] == chrom]
    #    
    #    if not chrom_passed.empty:
    #        chrom_passed.to_csv(chrom_bed_out, sep='\t', header=False, index=False)
    #        print(f"💾 Saved passing clusters to BED: {chrom_bed_out}")
    #        
    #    if not chrom_failed.empty:
    #        chrom_failed[cols].to_csv(chrom_csv_out, index=False)
    #        print(f"💾 Saved filtered-out audit report to CSV: {chrom_csv_out}")
    for chrom in df['chrom'].unique():
        chrom_folder = os.path.join(outdir, chrom)
        os.makedirs(chrom_folder, exist_ok=True)
        
        chrom_passed = passed_df[passed_df['chrom'] == chrom]
        chrom_failed = failed_df[failed_df['chrom'] == chrom]

        # Save BED
        if not chrom_passed.empty:
            bed_out = os.path.join(chrom_folder, out_bed)
            bed_cols = ['chrom', 'start_pos', 'end_pos', 'human_id']
            chrom_passed[bed_cols].assign(score=0, strand='+').to_csv(bed_out, sep='\t', header=False, index=False)

            # --- NEW: Save Detailed TSV (Two lines per cluster) ---
            details_out = os.path.join(chrom_folder, "filtered_clusters_details.tsv")
            with open(details_out, 'w') as f:
                for _, row in chrom_passed.iterrows():
                    # Line 1: Metadata
                    meta = (f"rep={row['rep_sequence']};gc={row['gc']};cons={row['cons']};"
                            f"med_len={row['med_len']};avg_raw_lev={row['avg_raw_lev']};"
                            f"avg_norm_lev={row['avg_norm_lev']};cov={row['cov']};"
                            f"region_size={row['region_size']};kmer_num={row['kmer_num']};"
                            f"avg_occ={row['avg_occ']};med_occ={row['med_occ']}")
                    f.write(f"{row['human_id']}\t{meta}\n")
                    # Line 2: K-mers
                    f.write(f"{row['human_id']}\t{row['kmers']}\n")
            
            print(f"💾 Saved details to TSV: {details_out}")

        # Save Rejection Audit
        if not chrom_failed.empty:
            csv_out = os.path.join(chrom_folder, out_csv)
            cols = ['human_id', 'filter_status', 'chrom', 'start_pos', 'end_pos', 
                    'region_size', 'med_occ', 'kmer_num', 'med_len', 'has_extended_kmer', 'rep_sequence']
            chrom_failed[cols].to_csv(csv_out, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter cluster database and output a BED file and Rejection CSV.")
    parser.add_argument("--db", required=True, help="Path to the SQLite database")
    parser.add_argument("--bed", default="filtered_clusters.bed", help="Output path for the clean BED file")
    parser.add_argument("--csv", default="rejected_clusters_audit.csv", help="Output path for the rejected clusters CSV")
    parser.add_argument("--chromosome", default=None, help="Specific chromosome to filter by (optional)")
    parser.add_argument("--outdir", required = True, help="Directory to save output files (will be organized by chromosome)")
    
    args = parser.parse_args()
    filter_database(args.db, args.bed, args.csv, args.chromosome, args.outdir)