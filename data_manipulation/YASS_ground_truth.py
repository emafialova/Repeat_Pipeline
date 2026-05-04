import os
import subprocess
import argparse

def run_yass(fasta_file, output_file, yass_path, diag_tol):
    """
    Runs YASS on a fasta file and filters for diagonal matches.
    
    diag_tol: max distance from diagonal to keep a hit (default 50 bp)
    """
    tmp_output = output_file + ".tmp"
    
    # Run YASS with tab-separated output (-d 3)
    cmd = [
        yass_path,
        fasta_file,
        fasta_file,
        "-o", tmp_output,  # correct way to specify output
        "-d", "3",          # tab-separated output
        "-w", "0"           # disable post-processing (recommended for long sequences)
    ]
    
    print(f"Running YASS: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    print("Parsing YASS output...")
    intervals = []
    with open(tmp_output) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            # Assuming YASS TSV format: query_start, query_end, subject_start, subject_end, score, evalue
            q_start, q_end, s_start, s_end = map(int, parts[:4])
            e_value = float(parts[8])
            
            # Keep only matches close to diagonal (query vs subject)
            if e_value > 1e-5:
                continue
            if q_start == s_start and q_end == s_end:
                continue
            if abs(q_start - s_start) <= diag_tol and abs(q_end - s_end) <= diag_tol:
                intervals.append((min(q_start, q_end), max(q_start, q_end)))
    
    if not intervals:
        print("No diagonal matches found.")
        return
    
    # Merge overlapping/adjacent intervals
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        # Only merge if intervals are overlapping OR within max_gap bases
        max_gap = 10  # adjustable threshold between hits
        if start <= last_end + max_gap:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

        #if start <= last_end + 1:  # merge if overlapping or adjacent
        #    merged[-1] = (last_start, max(last_end, end))
        #else:
        #    merged.append((start, end))
    
    # Write merged intervals to output
    with open(output_file, "w") as out:
        for start, end in merged:
            out.write(f"{start}\t{end}\n")
    
    print(f"Finished. Merged intervals saved to {output_file}")
    os.remove(tmp_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run YASS and extract diagonal repeat regions")
    parser.add_argument("fasta", help="Input fasta file")
    parser.add_argument("output", help="Output file for merged intervals")
    parser.add_argument("--yass", default="/Users/emafialova/Projects/Dependencies/yass/src/yass", help="Path to YASS executable")
    parser.add_argument("--diag_tol", type=int, default=100, help="Max distance from diagonal to keep hits")
    args = parser.parse_args()
    
    run_yass(args.fasta, args.output, args.yass, args.diag_tol)
