import argparse
import os
import csv
from collections import defaultdict
import sys
from itertools import combinations
import time
import multiprocessing as mp

csv.field_size_limit(sys.maxsize)

def find_reachable_leaves(start_node, edge_file):
    """Iterative DFS on a directed graph loaded from edge_file, return leaf nodes reachable from start_node."""
    graph = defaultdict(list)
    outdegree = defaultdict(int)

    with open(edge_file) as f:
        next(f)  # skip header 1
        next(f)  # skip header 2
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue
            src, label, dst = parts
            if label == "SIMILAR":
                continue
            graph[src].append(dst)
            outdegree[src] += 1
            if dst not in outdegree:
                outdegree[dst] = 0

    visited = set()
    leaf_nodes = set()
    stack = [start_node]

    while stack:
        node = stack.pop()
        if node in visited:
            continue
        visited.add(node)
        if outdegree[node] == 0:
            leaf_nodes.add(node)
        else:
            for neighbor in graph[node]:
                if neighbor not in visited:
                    stack.append(neighbor)
    return leaf_nodes

def filter_leaf_nodes(leaf_node_set):
    """Remove kmers that are strict subsequences of others."""
    leaf_node_list = list(leaf_node_set)
    filtered_list = []
    for i in leaf_node_list:
        is_subseq = False
        for j in leaf_node_list:
            if i != j and i in j:
                is_subseq = True
                break
        if not is_subseq:
            filtered_list.append(i)
    return filtered_list

def parse_region_from_tsv(tsv_file):
    """Parse first header line 'gene|start|end|region_num' to extract region_start and region_end."""
    with open(tsv_file) as f:
        header_line = f.readline().strip()
    parts = header_line.split('|')
    if len(parts) >= 3:
        try:
            region_start = int(parts[1])
            region_end = int(parts[2])
            return region_start, region_end
        except ValueError:
            pass
    return None, None

def extract_starting_kmers(tsv_file):
    """Extract kmers flagged as starting kmers."""
    starting_kmers = []
    with open(tsv_file) as f:
        lines = f.readlines()[1:]  # skip first header line
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        if row.get("starting_kmer_flag", "").strip() == "1":
            starting_kmers.append(row["kmer"].strip())
    unique_kmers = set(starting_kmers)
    starting_kmers = list(unique_kmers)
    return starting_kmers

def extract_kmer_info(tsv_file):
    """
    Extract all kmers and their info: positions, length, base_rel.freq.
    Returns dict: {kmer: {"length": int, "positions": [int,...], "base_rel_freq": float}}
    """
    kmer_info = {}
    with open(tsv_file) as f:
        lines = f.readlines()[1:]  # skip first header line
    reader = csv.DictReader(lines, delimiter="\t")
    required_cols = ["kmer", "length", "positions", "base_rel.freq"]
    for col in required_cols:
        if col not in reader.fieldnames:
            raise ValueError(f"'{col}' column not found in {tsv_file}")
    for row in reader:
        kmer = row["kmer"].strip()
        length = int(row["length"])
        base_rel_freq = float(row["base_rel.freq"])
        positions = [int(p) for p in row["positions"].split(",") if p.strip().isdigit()]
        kmer_info[kmer] = {
            "length": length,
            "positions": positions,
            "base_rel_freq": base_rel_freq
        }
    return kmer_info

def filter_true_final_kmers(leaf_nodes, abc_file):
    """
    Remove leaf nodes that appear in first column of the .abc file.
    Only truly final kmers are kept.
    """
    abc_kmers = set()
    with open(abc_file) as f:
        next(f)
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                abc_kmers.add(parts[0].strip())
    return [kmer for kmer in leaf_nodes if kmer not in abc_kmers]

def group_and_select_kmers(final_kmers, kmer_info):
    """
    Group kmers with same number of positions AND same absolute difference across positions.
    Pick one representative per group: longest kmer, then highest base_rel.freq.
    """
    groups = []
    for kmer1, kmer2 in combinations(final_kmers, 2):
        pos1 = kmer_info[kmer1]["positions"]
        pos2 = kmer_info[kmer2]["positions"]
        if len(pos1) != len(pos2):
            continue
        diffs = [abs(a - b) for a, b in zip(pos1, pos2)]
        if len(set(diffs)) == 1:
            merged = False
            for group in groups:
                if kmer1 in group or kmer2 in group:
                    group.update([kmer1, kmer2])
                    merged = True
                    break
            if not merged:
                groups.append(set([kmer1, kmer2]))

    grouped = set().union(*groups) if groups else set()
    for kmer in final_kmers:
        if kmer not in grouped:
            groups.append(set([kmer]))

    selected = []
    for group in groups:
        best = max(group, key=lambda k: (kmer_info[k]["length"], kmer_info[k]["base_rel_freq"]))
        selected.append(best)
    return selected

def filter_redundant(leaf_nodes, kmer_info):
    """
    Check each kmer against the original set of others.
    Safeguards against mutual redundancy by ensuring we always keep one of a redundant pair.
    """
    nodes = list(leaf_nodes)
    intervals = {}
    for k in nodes:
        length = kmer_info[k]["length"]
        positions = sorted(kmer_info[k]["positions"])
        intervals[k] = [(p, p + length) for p in positions]

    redundant_to_remove = set()
    
    for target in nodes:
        target_intervals = intervals[target]
        
        # Build merged intervals from all other kmers NOT yet marked for removal
        # This safeguard ensures we don't drop mutually redundant kmers
        others = [k for k in nodes if k != target and k not in redundant_to_remove]
        if not others:
            continue
            
        other_intervals = []
        for o in others:
            other_intervals.extend(intervals[o])
        
        other_intervals.sort()
        merged_blocks = []
        for curr in other_intervals:
            if not merged_blocks:
                merged_blocks.append(list(curr))
            else:
                last = merged_blocks[-1]
                if curr[0] <= last[1]:
                    last[1] = max(last[1], curr[1])
                else:
                    merged_blocks.append(list(curr))
        
        # Coverage check
        is_fully_covered = True
        for (t_start, t_end) in target_intervals:
            occ_covered = False
            for (m_start, m_end) in merged_blocks:
                if m_start <= t_start and m_end >= t_end:
                    occ_covered = True
                    break
            if not occ_covered:
                is_fully_covered = False
                break
        
        if is_fully_covered:
            redundant_to_remove.add(target)

    return [k for k in nodes if k not in redundant_to_remove]

def process_region(region_id, region_tsv_file, region_abc_file):
    print(f"  Started Processing region file: {region_id}")
    region_start, region_end = parse_region_from_tsv(region_tsv_file)
    if region_start is None or region_end is None:
        return None
    region_starting_kmers = extract_starting_kmers(region_tsv_file)
    if not region_starting_kmers:
        return None
    region_kmer_info = extract_kmer_info(region_tsv_file)
    kmer_dict = {}
    for kmer in region_starting_kmers:
        leaf_nodes = find_reachable_leaves(kmer, region_abc_file)
        leaf_nodes_filtered = filter_leaf_nodes(leaf_nodes)
        leaf_nodes_final = filter_true_final_kmers(leaf_nodes_filtered, region_abc_file)
        
        # Apply filters
        leaf_nodes_selected = group_and_select_kmers(leaf_nodes_final, region_kmer_info)
        leaf_nodes_final_set = filter_redundant(leaf_nodes_selected, region_kmer_info)
        kmer_dict[kmer] = leaf_nodes_final_set
    # check if all final kmers are unqiue and dont appear for multiple starting kmers    
    # if any final kmer does, keep it for the shorter starting kmer and remove it from the longer starting kmer (since they are in a parent-child relationship)
    unique_final_kmers = set()
    for start_kmer, final_kmers in kmer_dict.items():
        unique_final_kmers.update(final_kmers)
    if len(unique_final_kmers) < sum(len(fk) for fk in kmer_dict.values()):
        kmer_dict_unique = {}
        unique_list = list(unique_final_kmers) #[kmer1, kmer2]
        for unique_kmer in unique_list:
            found_in = []
            for start_kmer, final_kmers in kmer_dict.items():
                if unique_kmer in final_kmers:
                    found_in.append(start_kmer)
            if len(found_in) == 1:
                kmer_dict_unique[found_in[0]] = kmer_dict_unique.get(found_in[0], []) + [unique_kmer]
            else:
                best_start = min(found_in, key=lambda sk: region_kmer_info[sk]["length"])
                kmer_dict_unique[best_start] = kmer_dict_unique.get(best_start, []) + [unique_kmer]
        kmer_dict = kmer_dict_unique        
    print(f"  Finished Processing region file: {region_id}")
    return region_start, region_end, kmer_dict
    
def main():
    parser = argparse.ArgumentParser(description="Find final reachable nodes.")
    parser.add_argument("input_dir")
    parser.add_argument("output_file")
    parser.add_argument("--workers", type=int, default=4)

    if len(sys.argv) == 1:
        print("\nWelcome to the fourth step focused on final kmer extraction!")
        print("You need to provide input files and parameters.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    files = os.listdir(args.input_dir)
    region_ids = {f.split(".")[0] for f in files if f.startswith("region_") and (f.endswith(".tsv") or f.endswith(".abc"))}

    def extract_region_number(region_id):
        try: return int(region_id.split("_")[1])
        except: return float('inf')

    region_ids = sorted(region_ids, key=extract_region_number)
    tasks = []

    for region_id in region_ids:
        tsv_file = os.path.join(args.input_dir, region_id + ".tsv")
        abc_file = os.path.join(args.input_dir, region_id + ".abc")
        if os.path.exists(tsv_file) and os.path.exists(abc_file):
            tasks.append((region_id, tsv_file, abc_file))
    
    with mp.Pool(processes=args.workers) as pool:
        all_results = pool.starmap(process_region, tasks)

    with open(args.output_file, "w") as out_f:
        out_f.write("region_start\tregion_end\tkmer_dict\n")
        for res in all_results:
            if res:
                out_f.write(f"{res[0]}\t{res[1]}\t{res[2]}\n")
    
    print(f"\nProcessing completed. Output saved to: {args.output_file}")

if __name__ == "__main__":
    main()