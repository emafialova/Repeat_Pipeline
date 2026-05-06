import numpy as np
import pickle
import time
from suffix_array import SuffixArray
from Bio import SeqIO
import argparse
import re
from collections import defaultdict
import random
import sys

def generate_symbols():
    for i in range(5, 1000000):
        yield i

substitutions = {}
string_levels = []

def find_all_interesting_areas(LCP: list, k: int, m: int):
    possible = []
    real_areas = {}
    area_num = 1
    for index, i in enumerate(LCP):
        if i >= k:
            possible.append(index - 1)
        elif i < k and len(possible) < m:
            possible = []
        elif i < k and len(possible) >= m:
            real_areas[area_num] = possible
            possible = []
            area_num += 1
    return real_areas

def inner_substitute(most_freq_kmer: str, new_symbol: int, input_list: list):
    i = 0
    output_list = []
    while i < len(input_list):
        if input_list[i:i + len(most_freq_kmer)] == list(most_freq_kmer):
            output_list.append(new_symbol)
            i += len(most_freq_kmer)
        else:
            output_list.append(input_list[i])
            i += 1
    return output_list

def substitute(areas_dic: dict, input_string: list, substitutions: dict, k: int, def_generator):
    highest_freq = max([len(value) for value in areas_dic.values()])
    highest_key = [key for key, v in areas_dic.items() if len(v) == highest_freq]
    positions = areas_dic[highest_key[0]]
    most_freq_kmer = tuple(input_string[positions[0]:positions[0] + k])

    if most_freq_kmer not in substitutions:
        substitutions[most_freq_kmer] = next(def_generator)

    new_symbol = substitutions[most_freq_kmer]
    new_string = inner_substitute(most_freq_kmer, new_symbol, input_string)
    string_levels.append(new_string)
    return new_string

def outer_fce(string, k, i, def_generator):
    if len(string) <= 1:
        return string

    sa = SuffixArray(string)
    output = sa.suffix_array()
    lcp = sa.longest_common_prefix()
  
    areas_of_interest = find_all_interesting_areas(lcp, k, i)

    areas = {}
    for key, value in areas_of_interest.items():
        areas[key] = output[value[0]:value[-1] + 1]

    if len(areas) == 0:
        return string

    new_string = substitute(areas, string, substitutions, k, def_generator)
    return outer_fce(new_string, k, i, def_generator)

def get_corresponding_letters(number: int, substitutions_dictionary: dict, rules_dictionary: dict):
    corresponding_list = substitutions_dictionary[number]
    result = []
    for i in corresponding_list:
        if i in rules_dictionary:
            result.append(rules_dictionary[i])
        elif i in substitutions_dictionary:
            interm_result = get_corresponding_letters(i, substitutions_dictionary, rules_dictionary)
            result.extend(interm_result)
    return result

NT = [1, 2, 3, 4, 0]
starting_kmers = []

def simple_processing(input_string, rules, some_substitutions, output_file, window_start_absolute, kmer_dictionary):
    k = 1
    for j in some_substitutions.values():
        full_str = "".join(rules[m] if m in NT else "".join(get_corresponding_letters(m, some_substitutions, rules)) for m in j)
        starting_kmers.append(full_str)
        positions = [m.start() + window_start_absolute for m in re.finditer(f"(?={re.escape(full_str)})", input_string)]
        for i in positions:
            #output_handle.write(f"{full_str}\t{len(positions)}\t{len(full_str)}\t{positions[0]}\t{positions}\n")
            kmer_dictionary[full_str].append(i)

# =============== MAIN EXECUTION ===============

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Search for kmers in a FASTA file using suffix arrays with sliding windows.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file (sequence to be analyzed)")
    parser.add_argument("output_file", help="Base output path without extension")
    parser.add_argument("--min_kmer_length", type=int, default=10, help="Minimum k-mer length (default: 10)")
    parser.add_argument("--min_occurrences", type=int, default=10, help="Minimum occurrences to consider a k-mer (default: 10)")
    parser.add_argument("--window_size", type=int, default=10000, help="Size of sliding window in bp (default: 10000)")
    parser.add_argument("--overlap", type=int, default=5000, help="Size of window overlap in bp (default: 5000)")

    if len(sys.argv) == 1:
        print("\nWelcome to the first step used for searching kmers!")
        print("You need to provide input files and parameters.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    record = SeqIO.read(args.fasta_file, "fasta")
    intron = str(record.seq)
    print(f"Sequence loaded from {args.fasta_file}")
    start_time = time.time()

    step_size = args.window_size - args.overlap
    total_length = len(intron)

    rules = {"A": 1, "T": 2, "C": 3, "G": 4, "N": 0}
    rev_rules = {value: key for key, value in rules.items()}

    kmers = defaultdict(list)

    output_tsv_path = args.output_file
    with open(output_tsv_path, "w") as output_handle:
        output_handle.write("kmer\tcount\tlength\tfirst_position\tall_positions\n")

        for window_start in range(0, total_length, step_size):
            window_end = min(window_start + args.window_size, total_length)
            window_seq = intron[window_start:window_end]
            window_seq = re.sub(r'[^ATCG]', 'N', window_seq)
            window_seq_num = [rules[base] for base in window_seq]
            print(f"Processing window {window_start}-{window_end} ({window_end - window_start} nt)")
            generator = generate_symbols()
            substitutions.clear()
            string_levels.clear()
            processed_string = outer_fce(window_seq_num, args.min_kmer_length, args.min_occurrences, generator)

            rev_substitutions = {value: list(key) for key, value in substitutions.items()}
            simple_processing(window_seq, rev_rules, rev_substitutions, args.output_file, window_start, kmers)
        #print(kmers)
        kmers_set = {key:set(value) for key,value in kmers.items()}
        for kmer, positions in kmers_set.items():
            position_list = sorted(list(positions))
            output_handle.write(f"{kmer}\t{len(position_list)}\t{len(kmer)}\t{position_list[0]}\t{position_list}\n")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    print(f"Done. Total runtime: {minutes} min {seconds} sec ({elapsed_time:.2f} sec).")
    print(f"Results saved to {output_tsv_path}")
