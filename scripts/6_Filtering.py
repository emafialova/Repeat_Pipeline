#!/usr/bin/env python3
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd
from sqlalchemy import create_engine

def get_ucsc_engine(genome):
    """Creates a connection engine to the public UCSC Euro-Mirror MySQL server."""
    return create_engine(f"mysql+pymysql://genome@genome-euro-mysql.soe.ucsc.edu/{genome}")

def fetch_refseq_mapping(engine):
    """Fetches the UCSC-to-RefSeq chromosome conversion dictionary."""
    alias_query = "SELECT chrom AS ucsc_chrom, alias AS refseq_chrom FROM chromAlias WHERE source='refSeq';"
    try:
        df_alias = pd.read_sql(alias_query, engine)
        return dict(zip(df_alias['ucsc_chrom'], df_alias['refseq_chrom']))
    except Exception as e:
        print(f"Warning: Could not fetch chromAlias mapping ({e}). Falling back to UCSC naming.")
        return {}

def fetch_ucsc_tracks(genome):
    """Downloads all 3 reference tracks from UCSC and maps names to RefSeq format."""
    print(f"Connecting to UCSC Euro-Mirror for assembly: {genome}...")
    engine = get_ucsc_engine(genome)
    mapping = fetch_refseq_mapping(engine)
    
    # --- 1. Centromeres & Satellites Track ---
    print("-> Downloading Centromere & Satellite data...")
    sat_query = "SELECT genoName AS chrom, genoStart AS start, genoEnd AS end, repName AS name FROM rmsk WHERE repClass='Satellite';"
    cen_query = "SELECT chrom, chromStart AS start, chromEnd AS end, name FROM centromeres;"
    
    df_sat = pd.read_sql(sat_query, engine)
    try:
        df_cen = pd.read_sql(cen_query, engine)
    except Exception:
        df_cen = pd.DataFrame(columns=['chrom', 'start', 'end', 'name'])
    df_blacklist = pd.concat([df_sat, df_cen], ignore_index=True)
    
    # --- 2. TRF (Simple Repeats) Track ---
    print("-> Downloading TRF Simple Repeats track...")
    trf_query = "SELECT chrom, chromStart AS start, chromEnd AS end, sequence AS name FROM simpleRepeat;"
    df_trf = pd.read_sql(trf_query, engine)
    
    # --- 3. Coding Exons Track ---
    print("-> Downloading ncbiRefSeq Gene track and parsing coding exons...")
    exon_query = "SELECT chrom, exonStarts, exonEnds FROM ncbiRefSeq WHERE name LIKE 'NM_%%';" 
    # NM_ filters for validated coding transcripts, ignoring non-coding RNAs (NR_)
    df_raw_exons = pd.read_sql(exon_query, engine)
    
    # Parse the comma-separated exon coordinate blocks from SQL row format
    exon_rows = []
    for _, row in df_raw_exons.iterrows():
        # UCSC strings end with trailing commas, strip and split them
        starts = [int(x) for x in row['exonStarts'].decode('utf-8').strip(',').split(',') if x]
        ends = [int(x) for x in row['exonEnds'].decode('utf-8').strip(',').split(',') if x]
        for s, e in zip(starts, ends):
            exon_rows.append({'chrom': row['chrom'], 'start': s, 'end': e, 'name': 'exon'})
    df_exons = pd.DataFrame(exon_rows).drop_duplicates().reset_index(drop=True)
    
    # --- Format Types and Apply RefSeq Mapping ---
    for df in [df_blacklist, df_trf, df_exons]:
        if not df.empty:
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            if mapping:
                df['chrom'] = df['chrom'].map(mapping).fillna(df['chrom'])
                
    return df_blacklist, df_trf, df_exons

def merge_intervals(starts, ends):
    """Flattens overlapping genomic intervals to calculate true physical base-pair coverage."""
    if len(starts) == 0:
        return 0
    
    # Pair up and sort intervals strictly by their start positions
    intervals = sorted(zip(starts, ends), key=lambda x: x[0])
    
    merged_len = 0
    curr_start, curr_end = intervals[0]
    
    for next_start, next_end in intervals[1:]:
        if next_start < curr_end:
            # There is an overlap, push the current end boundary forward if needed
            curr_end = max(curr_end, next_end)
        else:
            # No overlap, count the accumulated length of the completed interval
            merged_len += (curr_end - curr_start)
            curr_start, curr_end = next_start, next_end
            
    # Add the final remaining interval
    merged_len += (curr_end - curr_start)
    return merged_len

def process_filtering(pipeline_df, df_blacklist, df_trf, df_exons, threshold):
    keep_rows = []
    dropped_report_rows = []
    padding = 100
    
    for chrom, group in pipeline_df.groupby('chrom'):
        b_chrom = df_blacklist[df_blacklist['chrom'] == chrom]
        t_chrom = df_trf[df_trf['chrom'] == chrom]
        e_chrom = df_exons[df_exons['chrom'] == chrom]
        
        for _, row in group.iterrows():
            p_start, p_end, p_id = row['start'], row['end'], row['id']
            p_len = p_end - p_start
            
            # --- FILTER 1: Strict Centromere / Satellite Check (With Padding) ---
            sat_overlaps = ((p_start - padding) < b_chrom['end']) & ((p_end + padding) > b_chrom['start'])
            if sat_overlaps.any():
                hit_names = ",".join(set(b_chrom['name'][sat_overlaps]))
                dropped_report_rows.append({
                    'chrom': chrom, 'start': p_start, 'end': p_end, 'cluster_id': p_id,
                    'filter_reason': 'Centromere/Satellite', 'matched_features': hit_names, 'overlap_fraction': 1.0
                })
                continue
                
            # --- FILTER 2: TRF Simple Repeat Fraction Check (with Interval Merging) ---
            trf_overlaps = (p_start < t_chrom['end']) & (p_end > t_chrom['start'])
            if trf_overlaps.any():
                o_starts = np.maximum(p_start, t_chrom['start'][trf_overlaps].values)
                o_ends = np.minimum(p_end, t_chrom['end'][trf_overlaps].values)
                
                # MERGE STACKED OVERLAPS HERE
                total_trf_overlap = merge_intervals(o_starts, o_ends)
                trf_fraction = total_trf_overlap / p_len
                
                if trf_fraction >= threshold:
                    raw_motifs = set(t_chrom['name'][trf_overlaps].dropna())
                    decoded_motifs = {m.decode('utf-8') if isinstance(m, bytes) else str(m) for m in raw_motifs}
                    hit_motifs = ",".join(decoded_motifs)[:50]
                    
                    dropped_report_rows.append({
                        'chrom': chrom, 'start': p_start, 'end': p_end, 'cluster_id': p_id,
                        'filter_reason': 'TRF_SimpleRepeat', 'matched_features': hit_motifs, 'overlap_fraction': round(trf_fraction, 3)
                    })
                    continue

            # --- FILTER 3: Coding Exon Fraction Check (with Interval Merging) ---
            exon_overlaps = (p_start < e_chrom['end']) & (p_end > e_chrom['start'])
            if exon_overlaps.any():
                o_starts = np.maximum(p_start, e_chrom['start'][exon_overlaps].values)
                o_ends = np.minimum(p_end, e_chrom['end'][exon_overlaps].values)
                
                # MERGE STACKED OVERLAPS HERE (In case transcripts overlap each other)
                total_exon_overlap = merge_intervals(o_starts, o_ends)
                exon_fraction = total_exon_overlap / p_len
                
                if exon_fraction >= threshold:
                    dropped_report_rows.append({
                        'chrom': chrom, 'start': p_start, 'end': p_end, 'cluster_id': p_id,
                        'filter_reason': 'Pure_Coding_Exon', 'matched_features': 'RefSeq_Exon', 'overlap_fraction': round(exon_fraction, 3)
                    })
                    continue
            
            # Eligible Cluster survives all filters
            keep_rows.append(row)
            
    return pd.DataFrame(keep_rows), pd.DataFrame(dropped_report_rows)

def process_filtering2(pipeline_df, df_blacklist, df_trf, df_exons, threshold):
    """
    Applies strict filter for satellites/centromeres, 
    and fractional overlap filtering for TRF and Exons.
    Returns: (dataframe of clean rows, list of records containing skipped reasons)
    """
    keep_rows = []
    dropped_report_rows = []
    padding = 100
    
    for chrom, group in pipeline_df.groupby('chrom'):
        b_chrom = df_blacklist[df_blacklist['chrom'] == chrom]
        t_chrom = df_trf[df_trf['chrom'] == chrom]
        e_chrom = df_exons[df_exons['chrom'] == chrom]
        
        for _, row in group.iterrows():
            p_start, p_end, p_id = row['start'], row['end'], row['id']
            p_len = p_end - p_start
            
            # --- FILTER 1: Strict Centromere / Satellite Check (With Padding) ---
            sat_overlaps = ((p_start - padding) < b_chrom['end']) & ((p_end + padding) > b_chrom['start'])
            if sat_overlaps.any():
                hit_names = ",".join(set(b_chrom['name'][sat_overlaps]))
                dropped_report_rows.append({
                    'chrom': chrom, 'start': p_start, 'end': p_end, 'cluster_id': p_id,
                    'filter_reason': 'Centromere/Satellite', 'matched_features': hit_names, 'overlap_fraction': 1.0
                })
                continue
                
            # --- FILTER 2: TRF Simple Repeat Fraction Check ---
            trf_overlaps = (p_start < t_chrom['end']) & (p_end > t_chrom['start'])
            if trf_overlaps.any():
                o_starts = np.maximum(p_start, t_chrom['start'][trf_overlaps].values)
                o_ends = np.minimum(p_end, t_chrom['end'][trf_overlaps].values)
                total_trf_overlap = (o_ends - o_starts).sum()
                trf_fraction = total_trf_overlap / p_len
                
                if trf_fraction >= threshold:
                    #hit_motifs = ",".join(set(t_chrom['name'][trf_overlaps].dropna()))[:50] # truncated for readability
                    # Convert bytes to strings dynamically if they are bytes objects
                    raw_motifs = set(t_chrom['name'][trf_overlaps].dropna())
                    decoded_motifs = {m.decode('utf-8') if isinstance(m, bytes) else str(m) for m in raw_motifs}
                    hit_motifs = ",".join(decoded_motifs)[:50] # truncated for readability
                    dropped_report_rows.append({
                        'chrom': chrom, 'start': p_start, 'end': p_end, 'cluster_id': p_id,
                        'filter_reason': 'TRF_SimpleRepeat', 'matched_features': hit_motifs, 'overlap_fraction': round(trf_fraction, 3)
                    })
                    continue

            # --- FILTER 3: Coding Exon Fraction Check (Protects Exon-Intron border expansions!) ---
            exon_overlaps = (p_start < e_chrom['end']) & (p_end > e_chrom['start'])
            if exon_overlaps.any():
                o_starts = np.maximum(p_start, e_chrom['start'][exon_overlaps].values)
                o_ends = np.minimum(p_end, e_chrom['end'][exon_overlaps].values)
                total_exon_overlap = (o_ends - o_starts).sum()
                exon_fraction = total_exon_overlap / p_len
                
                if exon_fraction >= threshold:
                    dropped_report_rows.append({
                        'chrom': chrom, 'start': p_start, 'end': p_end, 'cluster_id': p_id,
                        'filter_reason': 'Pure_Coding_Exon', 'matched_features': 'RefSeq_Exon', 'overlap_fraction': round(exon_fraction, 3)
                    })
                    continue
            
            # Elgible Cluster survives all filters
            keep_rows.append(row)
            
    return pd.DataFrame(keep_rows), pd.DataFrame(dropped_report_rows)

def main():
    parser = argparse.ArgumentParser(description="Advanced Nextflow-ready filtration of pipeline cluster files.")
    parser.add_argument("-g", "--genome", default="hg38", help="UCSC Genome ID (e.g., hg38, mm39)")
    parser.add_argument("-d", "--organism_dir", required=True, help="Path to your directory containing nested chromosome folders (e.g., path/to/HS/)")
    parser.add_argument("-o", "--clean_bed_output", required=True, help="Path to save the final clean, filtered BED file for coworkers")
    parser.add_argument("-r", "--report_output", required=True, help="Path to save the audit report of what was dropped and why (CSV)")
    parser.add_argument("-t", "--threshold", type=float, default=0.80, help="Overlap fraction threshold to drop simple repeats or exons (default: 0.80)")
    
    args = parser.parse_args()
    
    # 1. Download and compile database references from UCSC
    df_blacklist, df_trf, df_exons = fetch_ucsc_tracks(args.genome)
    
    # 2. Target pipeline 'filtered_clusters.tsv' instances across directories
    search_pattern = os.path.join(args.organism_dir, "**", "filtered_clusters.tsv")
    tsv_files = glob.glob(search_pattern, recursive=True)
    
    print(f"\nFound {len(tsv_files)} files across your chromosome directories to process.")
    if not tsv_files:
        print("Error: No 'filtered_clusters.tsv' instances recovered.")
        sys.exit(1)
        
    all_clean_clusters = []
    all_dropped_records = []
    
    # 3. Stream through individual chromosome files
    for tsv_path in tsv_files:
        if os.stat(tsv_path).st_size == 0:
            continue
            
        try:
            df_pipe = pd.read_csv(
                tsv_path, sep='\t', header=None, 
                names=['chrom', 'start', 'end', 'id'], usecols=[0, 1, 2, 3],
                dtype={'chrom': str, 'start': int, 'end': int, 'id': str}
            )
        except Exception as e:
            print(f"   Skipping file due to parsing errors: {tsv_path} ({e})")
            continue
            
        clean_df, dropped_df = process_filtering(df_pipe, df_blacklist, df_trf, df_exons, args.threshold)
        
        if not clean_df.empty:
            all_clean_clusters.append(clean_df)
        if not dropped_df.empty:
            all_dropped_records.append(dropped_df)
            
    # 4. Generate Final Outputs
    print("\n--- Compiling Outputs ---")
    
    # Save clean BED file for your coworkers
    if all_clean_clusters:
        final_clean_df = pd.concat(all_clean_clusters, ignore_index=True)
        final_clean_df = final_clean_df.sort_values(by=['chrom', 'start'])
        final_clean_df.to_csv(args.clean_bed_output, sep='\t', header=False, index=False)
        print(f"Success: Cleaned BED tracks written to: {args.clean_bed_output}")
        print(f"         Saved {len(final_clean_df)} high-confidence clusters.")
    else:
        print("Notice: No clusters survived the filtering criteria. Clean file is empty.")
        open(args.clean_bed_output, 'w').close() # touch empty file
        
    # Save audit drop report for your thesis
    if all_dropped_records:
        final_dropped_df = pd.concat(all_dropped_records, ignore_index=True)
        final_dropped_df.to_csv(args.report_output, index=False)
        print(f"Success: Detailed filtration audit log saved to: {args.report_output}")
        print(f"         Flagged and excluded {len(final_dropped_df)} regions of genomic background noise.")
    else:
        print("No regions were filtered out.")

if __name__ == "__main__":
    main()