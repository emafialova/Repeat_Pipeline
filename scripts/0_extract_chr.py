import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract whole chromosome or specific sequence from FASTA based on GFF coordinates or chromosome name.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file containing genomic sequences.")
    parser.add_argument("output_file", help="Path to the output FASTA file to save the extracted sequence.")
    parser.add_argument("chromosome_id", help="Chromosome ID to extract from the FASTA file.")
    parser.add_argument("--start", type=int, default=None, help="Start position (1-based) for subsequence extraction.")
    parser.add_argument("--end", type=int, default=None, help="End position (1-based, inclusive) for subsequence extraction.")

    if len(sys.argv) == 1:
        print("\nWelcome to a step used to extract whole or partial chromosomal sequence from a fasta file!")
        print("You need to provide input files and parameters.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    start = None
    end = None
    if args.start and args.end:
        start = args.start - 1      # Convert to 0-based index
        end = args.end          # end is exclusive in Python slicing

    chromosome_id = args.chromosome_id

    with open(f"{args.fasta_file}", mode='r') as soubor:
        data = soubor.read()
    entries = data.strip().split(">")

    chr_sequence = ""

    for entry in entries:
        if entry.startswith(chromosome_id):
            lines = entry.strip().split("\n")
            sequence = ''.join(lines[1:])  # remove header and join sequence lines
            chr_sequence = sequence
            break

    if start and end:
        subsequence = chr_sequence[start:end].upper()

    output_file = args.output_file
    output_name = output_file.strip().split("/")[-1].split(".")[0]
    
    with open(f"{output_file}", mode = "w") as soubor:
        if args.start and args.end:
            header = f">{output_name}|{chromosome_id},{args.start}-{args.end}\n"
            seq_to_record = subsequence
        else:
            header = f">{output_name}|{chromosome_id}\n"
            seq_to_record = chr_sequence.upper()
        soubor.write(header)
        soubor.write(seq_to_record)
