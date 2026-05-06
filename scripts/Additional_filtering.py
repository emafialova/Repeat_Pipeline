#!/usr/bin/env python3
import sqlite3
import pandas as pd
import argparse
import os

def filter_database(db_path, out_bed, out_csv, target_chrom, outdir):
    print(f"Connecting to database: {db_path}")
    conn = sqlite3.connect(db_path)
    
    # Extract main instance and cluster details
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

    # Group by internal_id and core_kmer, count the final kmers per core and find the MAX of those counts for each cluster.
    query_cores = """
        SELECT internal_id, 
               MAX(CASE WHEN LENGTH(kmer_seq) > LENGTH(core_kmer) THEN 1 ELSE 0 END) as has_extended_kmer,
               GROUP_CONCAT(kmer_seq) as kmers
        FROM final_kmer_info
        GROUP BY internal_id
    """
    core_df = pd.read_sql_query(query_cores, conn)
    conn.close()

    df = df.merge(core_df, on='internal_id', how='left')
    
    # Fill NaNs with 0 just in case a cluster has no records in final_kmer_info
    df['has_extended_kmer'] = df['has_extended_kmer'].fillna(0)
    df['kmers'] = df['kmers'].fillna("")
    print(f"Total raw clusters loaded: {len(df)}")

    # Build a list of failure reasons for every row
    def check_filters(row):
        reasons = []
        # Median occurrence >= 7
        if row['med_occ'] < 7:
            reasons.append(f"med_occ < 7 ({row['med_occ']})")
            
        # Region length >= 70
        if row['region_size'] < 70:
            reasons.append(f"region_size < 70 ({row['region_size']})")
            
        # At least one core must have > 1 final kmer
        if row['has_extended_kmer'] == 0:
            reasons.append("no final kmer extended beyond its core")
            
        return " | ".join(reasons) if reasons else "PASS"

    # Apply the logic
    df['filter_status'] = df.apply(check_filters, axis=1)

    # Split into Passed and Failed DataFrames
    passed_df = df[df['filter_status'] == "PASS"].copy()
    failed_df = df[df['filter_status'] != "PASS"].copy()

    print(f"Clusters passed: {len(passed_df)}")
    print(f"Clusters failed: {len(failed_df)}")

    # Prep BED File format for Passed Clusters
    bed_df = passed_df[['chrom', 'start_pos', 'end_pos', 'human_id']].copy()
    bed_df['score'] = 0
    bed_df['strand'] = '+'
    bed_df = bed_df[['chrom', 'start_pos', 'end_pos', 'human_id', 'score', 'strand']]
    
    cols = ['human_id', 'filter_status', 'chrom', 'start_pos', 'end_pos', 
            'region_size', 'med_occ', 'kmer_num', 'med_len', 'has_extended_kmer', 'rep_sequence']

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

            # Save Detailed TSV 
            details_out = os.path.join(chrom_folder, "FINAL_clusters_details.tsv")
            with open(details_out, 'w') as f:
                for _, row in chrom_passed.iterrows():
                    meta = (f"rep={row['rep_sequence']};gc={row['gc']};cons={row['cons']};"
                            f"med_len={row['med_len']};avg_raw_lev={row['avg_raw_lev']};"
                            f"avg_norm_lev={row['avg_norm_lev']};cov={row['cov']};"
                            f"region_size={row['region_size']};kmer_num={row['kmer_num']};"
                            f"avg_occ={row['avg_occ']};med_occ={row['med_occ']}")
                    f.write(f"{row['human_id']}\t{meta}\n")
                    f.write(f"{row['human_id']}\t{row['kmers']}\n")
            
            print(f"Saved details to TSV: {details_out}")

        # Save Rejection Audit
        if not chrom_failed.empty:
            csv_out = os.path.join(chrom_folder, out_csv)
            cols = ['human_id', 'filter_status', 'chrom', 'start_pos', 'end_pos', 
                    'region_size', 'med_occ', 'kmer_num', 'med_len', 'has_extended_kmer', 'rep_sequence']
            chrom_failed[cols].to_csv(csv_out, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter cluster database and output a BED file and Rejection CSV.")
    parser.add_argument("--db", required=True, help="Path to the SQLite database")
    parser.add_argument("--bed", default="FINAL_clusters.bed", help="Output path for the clean BED file")
    parser.add_argument("--csv", default="rejected_clusters.csv", help="Output path for the rejected clusters CSV")
    parser.add_argument("--chromosome", default=None, help="Specific chromosome to filter by (optional)")
    parser.add_argument("--outdir", required = True, help="Directory to save output files (will be organized by chromosome)")
    
    args = parser.parse_args()
    filter_database(args.db, args.bed, args.csv, args.chromosome, args.outdir)