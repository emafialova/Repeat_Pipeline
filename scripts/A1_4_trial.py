#!/usr/bin/env python3
import numpy as np
from Levenshtein import distance as lev_distance
import argparse
import ast
from sklearn.cluster import DBSCAN
import sqlite3
import json
import hashlib
import re
import os
import sys
import time
from collections import Counter, defaultdict
from Bio import SeqIO
import multiprocessing as mp
from Bio import Align

# ------------------ ALIGNER SETUP ------------------
ALIGNER = Align.PairwiseAligner()
ALIGNER.mode = 'local' 
ALIGNER.match_score = 2
ALIGNER.mismatch_score = -1
ALIGNER.open_gap_score = -2
ALIGNER.extend_gap_score = -1

# ------------------ ID GENERATION HELPERS ------------------

def load_ready_mapping(mapping_path, target_org):
    target_org = target_org.lower()
    try:
        with open(mapping_path, 'r') as f:
            next(f) 
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                if parts[2].lower() == target_org:
                    return parts[1]
    except Exception as e:
        print(f"Error loading organism mapping: {e}")
        sys.exit(1)
    sys.exit(1)

def to_base36(n):
    chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if n < 36: return f"0{chars[n]}"
    high, low = n // 36, n % 36
    if high > 35: return "ZZ" 
    return f"{chars[high]}{chars[low]}"

def parse_fasta_chroms(fasta_path, pattern):
    info_regex = re.compile(r"chromosome\s+([^\s,;]+)", re.IGNORECASE)
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>') and pattern in line:
                header = line[1:].strip()
                match = info_regex.search(header)
                if match:
                    val = match.group(1).strip(',.;').upper()
                    if val.isdigit():
                        return f"{int(val):02d}"
                    elif val in ['X', 'Y', 'W', 'Z', 'M', 'C', 'MT', 'CP']:
                        return f"0{val}" if len(val) == 1 else val[:2]
                    else:
                        return val[:2].zfill(2)
                break
    return "00"

def get_or_update_hash_map(map_path, rep_seq):
    seq = rep_seq.strip().upper()
    # Generate the actual hashes for the current sequence
    current_short_hash = hashlib.md5(seq.encode()).hexdigest()[:8].upper()
    current_full_hash = hashlib.md5(seq.encode()).hexdigest().upper()
    
    # sequence -> (short, full)
    seq_to_hashes = {}
    # short_hash -> sequence (for collision detection)
    short_hash_to_seq = {}

    if os.path.exists(map_path):
        with open(map_path, 'r') as f:
            next(f) # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                s, h, fh = parts
                seq_to_hashes[s] = (h, fh)
                short_hash_to_seq[h] = s

    # 1. If sequence already exists, return the stored TUPLE
    if seq in seq_to_hashes:
        return seq_to_hashes[seq]

    # 2. Collision Check: Does this short hash belong to a DIFFERENT sequence?
    if current_short_hash in short_hash_to_seq:
        print(f"CRITICAL ERROR: Hash collision for {current_short_hash}!")
        sys.exit(1)

    # 3. New entry: Append to file
    file_mode = 'a' if os.path.exists(map_path) else 'w'
    with open(map_path, file_mode) as f:
        if file_mode == 'w':
            f.write("sequence\thash\tfull_hash\n")
        f.write(f"{seq}\t{current_short_hash}\t{current_full_hash}\n")
    
    return current_short_hash, current_full_hash

def get_db_hashes(cursor, rep_seq):
    """
    Checks database for existing hash. If not found, generates and saves it.
    This is parallel-safe when wrapped in a transaction.
    """
    seq = rep_seq.strip().upper()
    cursor.execute("SELECT short_hash, full_hash FROM hashes WHERE sequence = ?", (seq,))
    row = cursor.fetchone()
    
    if row:
        return row[0], row[1]
    
    # Generate new hashes
    full_h = hashlib.md5(seq.encode()).hexdigest().upper()
    short_h = full_h[:8]
    
    cursor.execute("INSERT INTO hashes (sequence, short_hash, full_hash) VALUES (?, ?, ?)", 
                   (seq, short_h, full_h))
    return short_h, full_h

# ------------------ DATA LOADING & ANALYSIS ------------------

def extract_all(final_nodes_file):
    """
    Reads final_nodes.tsv-like file: header then lines of 'start<TAB>end<TAB>dict'.
    Returns:
      all_regions: dict {region_name: (region_start:int, region_end:int, core_dict:dict)}
      region_list: list of region names
    """
    all_regions, region_list = {}, []
    with open(final_nodes_file, mode='r') as f:
        next(f)
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                start, end, kmers = line.split("\t")
            except ValueError:
                continue
            start, end, kmers = line.split("\t")
            region_name = f"region_{i}"
            all_regions[region_name] = (int(start), int(end), ast.literal_eval(kmers))
            region_list.append(region_name)
    return all_regions, region_list

def parse_tsv_info(region_tsv_file):
    """
    Reads a region_Y.tsv produced by your program:
    header line (gene|start|end|region_id)
    header2 (column names)
    then lines: kmer\tcount\tlength\tcompactness\tpositions\tstarting_kmer_flag
    Returns dict: kmer -> (count:int, positions: List[int])
    """
    kmer_tsv_dict = {}
    with open(region_tsv_file, mode='r') as f:
        try:
            next(f); next(f)
        except StopIteration: return kmer_tsv_dict
        for ln in f:
            parts = ln.strip().split("\t")
            if len(parts) < 5: continue
            kmer = parts[0]
            try: count = int(float(parts[1]))
            except: count = 0
            # Clean positions field
            positions_field = parts[4]
            # positions field expected like: "[12,34,56]" or "[]"
            positions = []
            p = positions_field.strip()
            if p.startswith("[") and p.endswith("]"):
                inner = p[1:-1].strip()
                if inner:
                    try:
                        positions = [int(x) for x in inner.split(",") if x.strip() != ""]
                    except Exception:
                        # fallback: try to split on non-numeric
                        import re
                        nums = re.findall(r"\d+", inner)
                        positions = [int(x) for x in nums]
            else:
                # If it's not bracketed, try to parse numbers inside
                import re
                nums = re.findall(r"\d+", p)
                positions = [int(x) for x in nums]
            kmer_tsv_dict[kmer] = (count, positions)
    return kmer_tsv_dict

def filter_kmers_by_abc(region_kmer_dict, abc_path):
    """
    Filters out kmers that appear in the FIRST column of the abc file.
    
    Args:
        region_kmer_dict (dict): {'core_kmer': [list_of_final_kmers]}
        abc_path (str): path to region_*.abc file (3-column TSV)
    
    Returns:
        dict: same structure as input, but with kmers removed if found in abc file's first column
    """
    if not os.path.exists(abc_path): return region_kmer_dict
    abc_firstcol = set()
    with open(abc_path, 'r') as f:
        for line in f:
            if line.strip():
                first = line.split('\t', 1)[0]  # only first column
                abc_firstcol.add(first)
    filtered_dict = {}
    for core, final_list in region_kmer_dict.items():
        new_list = [k for k in final_list if k not in abc_firstcol]
        if new_list:  # only keep cores that still have at least one final
            filtered_dict[core] = new_list

    return filtered_dict

def create_distance_matrix(sequences):
    N = len(sequences)
    dist_matrix = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(i+1, N):
            s1, s2 = sequences[i], sequences[j]
            sw_score = ALIGNER.score(s1, s2)
            max_possible = min(len(s1), len(s2)) * ALIGNER.match_score
            dist_matrix[i, j] = dist_matrix[j, i] = 1.0 - (max(0.0, sw_score) / max_possible) if max_possible > 0 else 1.0
    return dist_matrix

def find_representative(sequence_list, tsv_kmer_dict, dist_matrix, idxs):
    if not idxs: return None
    cluster_seqs = [sequence_list[i] for i in idxs]
    freqs = []
    for seq in cluster_seqs:
        val = tsv_kmer_dict.get(seq)
        if val is None:
            freqs.append(0)
        else:
            try:
                freqs.append(int(val[0]))
            except Exception:
                freqs.append(0)
    if freqs:
        max_freq = max(freqs)
        if freqs.count(max_freq) == 1 and max_freq > 0:
            rep_idx_local = freqs.index(max_freq)
            return cluster_seqs[rep_idx_local]
    # medoid fallback: compute average distances within cluster
    if dist_matrix is None:
        return cluster_seqs[0]
    sub = dist_matrix[np.ix_(idxs, idxs)]
    avg = sub.mean(axis=1)
    medoid_local_idx = int(np.argmin(avg))
    return cluster_seqs[medoid_local_idx]

def run_dbscan(dist_matrix, eps, min_samples):
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    clustering.fit(dist_matrix)
    return clustering.labels_

def gc_content_region(sequence):
    if not sequence:
        return 0.0
    gc_abs = sum(1 for c in sequence.upper() if c in "GC")
    return 100.0 * gc_abs / len(sequence)

def avg_sw_stats(sequence_list, idxs, dist_matrix):
    if len(idxs) < 2: return 0.0, 0.0
    sub = dist_matrix[np.ix_(idxs, idxs)]
    triu_vals = sub[np.triu_indices(len(idxs), k=1)]
    avg_norm = float(triu_vals.mean()) if triu_vals.size > 0 else 0.0

    total_raw_sw = 0
    count = 0
    for i in range(len(idxs)):
        for j in range(i+1, len(idxs)):
            s1 = sequence_list[idxs[i]]
            s2 = sequence_list[idxs[j]]
            raw_sw_score = ALIGNER.score(s1, s2) 
            total_raw_sw += raw_sw_score
            count += 1
            
    avg_raw_sw = total_raw_sw / count if count > 0 else 0.0
    return (avg_raw_sw, avg_norm)

def merge_regions(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for interval in sorted_intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    return merged

def calc_coverage(tsv_kmers, sequence_list, idxs):
    all_cluster_seqs = [sequence_list[i] for i in idxs]
    coverage = 0
    covered_regions = []
    for sequence in all_cluster_seqs:
        vals = tsv_kmers.get(sequence)
        #print(sequence, vals)
        length = len(sequence)
        positions = vals[1]
        for i in positions:
            covered_regions.append((i, i+length-1))
    merged_regions = merge_regions(covered_regions)
    for region in merged_regions:
        start = region[0]
        stop = region[1]
        coverage += stop - start + 1
    return coverage, merged_regions[0][0], merged_regions[-1][1]

def calc_conservation(sequence_list, idxs, rep):
    single_conservations = []   
    for i in idxs:
        s = sequence_list[i]
        d = lev_distance(s, rep)
        denom = max(len(s), 1)
        d_norm = d / denom
        cons = 1.0 - d_norm
        single_conservations.append(cons)
    avg_cons = np.mean(single_conservations) if single_conservations else 0
    return avg_cons

def merge_clusters_by_overlap(cluster_info):
    if not cluster_info: return []
    intervals = sorted([{'s': v['start'], 'e': v['end'], 'c': k, 'i': v['idxs']} for k, v in cluster_info.items()], key=lambda x: x['s'])
    islands = []
    curr = {'start': intervals[0]['s'], 'end': intervals[0]['e'], 'kmer_idxs': list(intervals[0]['i']), 'merged_cids': [intervals[0]['c']]}
    for iv in intervals[1:]:
        if iv['s'] <= curr['end']:
            curr['end'] = max(curr['end'], iv['e'])
            curr['kmer_idxs'].extend(iv['i'])
            curr['merged_cids'].append(iv['c'])
        else:
            islands.append(curr)
            curr = {'start': iv['s'], 'end': iv['e'], 'kmer_idxs': list(iv['i']), 'merged_cids': [iv['c']]}
    islands.append(curr)
    return islands

# ------------------ ANALYSIS ENGINE ------------------

def analyze_region(region_index, region_tsv, region_directory, region, complete_sequence, sequence_id, region_name, region_start, region_end, region_kmer_dict, eps, min_s):
    print(f"  Started Processing region file: {region_index}")
    tsv_kmer_dict = parse_tsv_info(region_tsv)
    region_abc = os.path.join(region_directory, region + ".abc")
    region_kmer_dict = filter_kmers_by_abc(region_kmer_dict, region_abc)
    
    sequence_list, core_of_seq = [], []
    for core, finals in region_kmer_dict.items():
        for s in finals:
            sequence_list.append(s); core_of_seq.append(core)
    if not sequence_list: return None

    dist_matrix = create_distance_matrix(sequence_list)
    labels = run_dbscan(dist_matrix, eps, min_s)
    
    initial_clusters = {}
    for cid in set(labels):
        idxs = [i for i, l in enumerate(labels) if l == cid]
        if not idxs:
            continue
        _, m_start, m_end = calc_coverage(tsv_kmer_dict, sequence_list, idxs)
        initial_clusters[cid] = {'start': m_start, 'end': m_end, 'idxs': idxs}

    islands_raw = merge_clusters_by_overlap(initial_clusters)
    final_island_info = {}

    for i, island in enumerate(islands_raw, start=1):
        idxs = island['kmer_idxs']
        rep = find_representative(sequence_list, tsv_kmer_dict, dist_matrix, idxs)
        cov_nom, m_start, m_end = calc_coverage(tsv_kmer_dict, sequence_list, idxs)
        avg_raw, avg_norm = avg_sw_stats(sequence_list, idxs, dist_matrix)
        counts = [tsv_kmer_dict.get(sequence_list[idx], (0,))[0] for idx in idxs]

        island_offset = m_start
        # Extract Kmer Data for DB
        kmer_data = []
        unique_cores_in_island = set()
        for idx in idxs:
            seq = sequence_list[idx]
            core_seq = core_of_seq[idx]
            unique_cores_in_island.add(core_seq)
            original_positions = tsv_kmer_dict.get(seq, (0, []))[1]
            updated_positions = [pos - island_offset for pos in original_positions]
            kmer_data.append({
                'seq': seq,
                'pos': updated_positions,
                'is_rep': 1 if seq == rep else 0,
                'core': core_seq
            })

        core_data = []
        for c_seq in unique_cores_in_island:
            core_original_positions = tsv_kmer_dict.get(c_seq, (0, []))[1]
            core_data.append({
                'seq': c_seq,
                'pos': [pos - island_offset for pos in core_original_positions]
            })

        r_size = m_end - m_start + 1
        final_island_info[f"I{i}"] = {
            'rep': rep, 'start': m_start, 'end': m_end, 'region_size': r_size,
            'gc': gc_content_region(complete_sequence[region_start:region_end+1]),
            'avg_raw_sw': avg_raw, 'avg_norm_sw': avg_norm,
            'coverage': cov_nom / max(r_size, 1), 'avg_occ': np.mean(counts),
            'median_occ': np.median(counts), 'kmer_num': len(idxs),
            'median_length': np.median([len(sequence_list[x]) for x in idxs]),
            'conservation': calc_conservation(sequence_list, idxs, rep),
            'all_cluster_seqs': [sequence_list[x] for x in idxs],
            'kmer_data': kmer_data,
            'core_data': core_data
        }
    return region_index, region_start, region_end, final_island_info, region_index

# ------------------ MAIN EXECUTION ------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("final_nodes")
    parser.add_argument("region_directory")
    parser.add_argument("sequence_fasta")
    parser.add_argument("chromosome")
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--org", required=True)
    parser.add_argument("--db_path", default="genome_clusters.db")
    parser.add_argument("--output", default="clusters_islands.tsv")
    parser.add_argument("--dbscan_eps", type=float, default=0.45)
    parser.add_argument("--dbscan_min_samples", type=int, default=2)
    parser.add_argument("--workers", type=int, default=4)
    
    if len(sys.argv) == 1:
        print("\nWelcome to the Kmer Clustering Tool!")
        print("You need to provide input files and parameters.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    # Initialize Context
    org_code = load_ready_mapping(args.mapping, args.org)
    c_code = parse_fasta_chroms(args.sequence_fasta, args.chromosome)
    #c_code = chrom_map.get(args.chromosome, "00")
    suffix_tracker = Counter()

    # Run Multi-Processing
    all_regions, region_list = extract_all(args.final_nodes)
    #record = SeqIO.read(args.sequence_fasta, "fasta")
    complete_seq = None
    for record in SeqIO.parse(args.sequence_fasta, "fasta"):
        if args.chromosome in record.id or args.chromosome in record.description:
            complete_seq = str(record.seq); break
    if complete_seq is None: print(f"Error: {args.chromosome} not found"); sys.exit(1)
    
    tasks = []
    for idx, region in enumerate(region_list, start=1):
        r_start, r_end, r_dict = all_regions[region]
        r_tsv = os.path.join(args.region_directory, f"{region}.tsv")
        if not os.path.exists(r_tsv):
            print(f"Warning: region tsv not found for {r_tsv}!")
            continue
        tasks.append((idx, r_tsv, args.region_directory, region, complete_seq, args.chromosome, region, r_start, r_end, r_dict, args.dbscan_eps, args.dbscan_min_samples))

    with mp.Pool(processes=args.workers) as pool:
        results = pool.starmap(analyze_region, tasks)
    
    results = sorted([r for r in results if r], key=lambda x: x[0])

    # ------------------ DB POPULATION & TRANSFORMED OUTPUT ------------------
    # SQL Transaction with 60s timeout for parallel chromosome runs
    conn = sqlite3.connect(args.db_path, timeout=60.0)
    cursor = conn.cursor()
    
    with open(args.output, "w") as out_f, \
         open(args.output.replace(".tsv", "_details.tsv"), "w") as out_d:
        
        cursor.execute("BEGIN TRANSACTION")
        for _, r_start, r_end, island_info, r_id in results:
            for i_id, info in island_info.items():
                abs_start, abs_end = r_start + info['start'], r_start + info['end']
                #abs_start, abs_end = r_start, r_end
                
                # ID Logic: GAL01-A1B2C3D4.01
                #short_h, full_h = get_or_update_hash_map(args.hash_map, info['rep'])
                short_h, full_h = get_db_hashes(cursor, info['rep'])
                base_id = f"{org_code}{c_code}-{short_h}"
                suffix_tracker[base_id] += 1
                human_id = f"{base_id}.{to_base36(suffix_tracker[base_id])}"
                internal_id = f"{org_code}{c_code}-{full_h}.{to_base36(suffix_tracker[base_id])}"

                # Write Files
                out_f.write(f"{args.chromosome}\t{abs_start}\t{abs_end}\t{human_id}\t0\t+\n")
                
                metric_str = (f"rep={info['rep']};gc={info['gc']:.1f};cons={info['conservation']:.3f};"
                              f"med_len={info['median_length']};avg_raw_lev={info['avg_raw_sw']:.2f};"
                              f"avg_norm_lev={info['avg_norm_sw']:.3f};cov={info['coverage']:.3f};"
                              f"region_size={info['region_size']};kmer_num={info['kmer_num']};"
                              f"avg_occ={info['avg_occ']:.2f};med_occ={info['median_occ']:.2f}")
                
                out_d.write(f"{human_id}\t{metric_str}\n{human_id}\t{','.join(info['all_cluster_seqs'])}\n")

                # DB Insertion
                cursor.execute("INSERT OR REPLACE INTO instances VALUES (?, ?, ?, ?, ?, ?, ?)",
                               (internal_id, args.chromosome, abs_start, abs_end, human_id, short_h, full_h))
                
                cursor.execute("INSERT OR REPLACE INTO cluster_details VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                               (internal_id, info['rep'], info['gc'], info['conservation'], info['median_length'], 
                                info['avg_raw_sw'], info['avg_norm_sw'], info['coverage'], info['region_size'], 
                                info['kmer_num'], info['avg_occ'], info['median_occ']))
                
                cursor.execute("INSERT OR REPLACE INTO kmers VALUES (?, ?)", (internal_id, ','.join(info['all_cluster_seqs'])))

                for k in info['kmer_data']:
                    cursor.execute("INSERT INTO final_kmer_info VALUES (?, ?, ?, ?, ?)",
                                   (internal_id, k['seq'], json.dumps(k['pos']), k['is_rep'], k['core']))
                for c in info['core_data']:
                    cursor.execute("INSERT INTO core_kmer_info VALUES (?, ?, ?)",
                                   (internal_id, c['seq'], json.dumps(c['pos'])))

        cursor.execute("COMMIT")
    conn.close()
    print(f"✅ Processed {args.chromosome} and updated {args.db_path}")

if __name__ == "__main__":
    main()
