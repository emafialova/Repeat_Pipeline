import argparse
import time
import ast  # safer literal evaluation
from Bio import SeqIO
import sys

def is_low_complexity_kmer(kmer: str) -> bool:
    """
    Return True if kmer is of low complexity (all the same character).
    """
    return len(set(kmer)) == 1

def split_kmers(input_row: str, position_dict, threshold=1000):
    kmer, _, _, _, all_positions = input_row.strip().split("\t")
    
    if is_low_complexity_kmer(kmer):
        return position_dict  # skip this kmer entirely
    # Safely convert string to list of integers
    positions_list = ast.literal_eval(all_positions)
    
    groups = []
    current_group = [positions_list[0]]

    for i in range(1, len(positions_list)):
        if positions_list[i] - positions_list[i - 1] <= threshold:
            current_group.append(positions_list[i])
        else:
            # Gap > threshold, save the current group
            groups.append(current_group)
            current_group = [positions_list[i]]

    # Save the final group
    if current_group:
        groups.append(current_group)

    position_dict[kmer] = groups
    return position_dict

def merge_split_kmers(position_dict, min_distance=30, min_occurrences=10):
    """
    Merge overlapping or adjacent intervals from the position_dict,
    filter based on length (min_distance),
    and keep only those with at least min_occurrences positions inside.
    """
    all_positions_flat = []
    for kmer, groups in position_dict.items():
        for group in groups:
            all_positions_flat.extend(group)

    all_positions_flat.sort()

    # Build intervals for merging
    all_groups = []
    for kmer, groups in position_dict.items():
        for group in groups:
            all_groups.append((min(group), max(group)))
    # Sort intervals by start position
    all_groups.sort(key=lambda x: x[0])

    merged = []
    for interval in all_groups:
        if not merged:
            merged.append(list(interval))
        else:
            last = merged[-1]
            # If overlapping or touching, merge
            if last[1] >= interval[0] - 1:
                last[1] = max(last[1], interval[1])
            else:
                merged.append(list(interval))

    filtered_merged = []
    for m in merged:
        start, end = m
        # Filter on length
        if end - start < min_distance:
            continue
        # Count the number of positions within the interval
        count = sum(start <= pos <= end for pos in all_positions_flat)
        if count >= min_occurrences:
            filtered_merged.append((start, end))
    return filtered_merged

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze the frequent kmers to determine interesting regions.")
    parser.add_argument("kmer_file", help="Path to the input file containing kmers.")
    parser.add_argument("output_file", help="Base path for the output (no extension), e.g., ./A1/chr31/A1_kmers_full_gene_02")
    parser.add_argument("analyzed_sequence", help = "Path to the original analyzed sequence (FASTA format!).")
    parser.add_argument("--threshold", type=int, default=1000, help="Threshold for position gap splitting (default: 1000)")
    
    if len(sys.argv) == 1:
        print("\nWelcome to the second step used for the creation of regions of interest!")
        print("You need to provide input files and parameters.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    start_time = time.time()

    with open(args.kmer_file, "r") as f:
        kmer_records = f.readlines()[1:]  # skip header

    position_dictionary = {}
    for line in kmer_records:
        split_kmers(line, position_dictionary, threshold=args.threshold)

    merged_intervals = merge_split_kmers(position_dictionary)
    print(merged_intervals)

    with open(f"{args.output_file}.tsv", "w") as out_f:
        out_f.write(f"{args.kmer_file} analysis of found kmers - transformation into areas of interest\n")
        out_f.write("start\tstop\tkmers\n")
        for interval in merged_intervals:
            start, end = interval
            contributing_kmers = []
            for kmer, groups in position_dictionary.items():
                for group in groups:
                    if any(start <= pos <= end for pos in group):
                        contributing_kmers.append(kmer)
                        break  # avoid duplicates if multiple groups from same kmer
            out_f.write(f"{start}\t{end}\t{contributing_kmers}\n")

    # Uncomment the following block if you want to visualize the merged intervals on the original sequence (requires matplotlib)
    
    #sequence = SeqIO.read(args.analyzed_sequence, "fasta")
    #seq = str(sequence.seq)
    
    #plt.figure(figsize=(15, 2))
    #plt.hlines(y=0, xmin=0, xmax=len(seq), color='gray', linewidth=2)

    #for interval in merged_intervals:
        #y_level = 20
        #plt.plot([interval[0], interval[1]], [y_level, y_level], lw=6, color="blue")

    #plt.xlabel("Genomic position (bp)")
    #plt.title("Merged Repeat Intervals")
    #plt.yticks(range(len(merged_intervals) + 1))
    #plt.tight_layout()
    #plt.savefig(args.output_file + "_merged_intervals.png")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Total runtime: {minutes} minutes {seconds} seconds ({elapsed_time:.2f} seconds)")
    print(f"Outputs saved to {args.output_file}.tsv")
