import argparse
import sys
from Bio import SeqIO
import numpy as np

def confusion_matrix2(reference_file, pipeline_file, sequence_length):
    TP = 0
    FP = 0
    FN = 0      
    TN = 0
    # Load reference intervals
    reference_intervals = []
    with open(reference_file) as ref_f:
        for line in ref_f:
            start, end = line.strip().split("\t")
            reference_intervals.append((int(start), int(end)))
    # Load pipeline intervals from bed file
    pipeline_intervals = []
    with open(pipeline_file) as pipe_f:
        for line in pipe_f:
            _, start, end, _, _, _ = line.strip().split("\t")
            pipeline_intervals.append((int(start), int(end)))
    # Calculate TP, FP, FN, TN
    previous_end = 0 # prepare to later take care of gaps between intervals
    for i, (ref_start, ref_end) in enumerate(reference_intervals):
        if i == 0:  # take care of bases before first reference interval
            TN += ref_start  # assuming genome starts at position 0
        elif i < len(reference_intervals) - 1:
            covered_bases = 0  # track how many bases of this reference are covered
            for pipe_start, pipe_end in pipeline_intervals:
                # Case 1: pipeline fully within reference
                if pipe_start >= ref_start and pipe_end <= ref_end:
                    overlap = pipe_end - pipe_start + 1
                    TP += overlap
                    covered_bases += overlap
                    # pipeline_intervals.remove((pipe_start, pipe_end))  # optional
                # Case 2: partial overlap on right
                elif pipe_start >= ref_start and pipe_start <= ref_end and pipe_end > ref_end:
                    overlap = ref_end - pipe_start + 1
                    TP += overlap
                    FP += pipe_end - ref_end
                    covered_bases += overlap
                # Case 3: partial overlap on left
                elif pipe_start < ref_start and pipe_end >= ref_start and pipe_end <= ref_end:
                    overlap = pipe_end - ref_start + 1
                    TP += overlap
                    FP += ref_start - pipe_start
                    covered_bases += overlap
                elif pipe_start < ref_start and pipe_end > ref_end:
                    TP += ref_end - ref_start + 1
                    FP += (ref_start - pipe_start) + (pipe_end - ref_end)
                    covered_bases += ref_end - ref_start + 1
            # Now compute FN: reference length minus covered bases
            ref_length = ref_end - ref_start + 1
            FN += max(0, ref_length - covered_bases)
            TN += ref_start - previous_end - 1
            previous_end = ref_end
        elif i == len(reference_intervals) - 1:  # take care of remaining bases
            TN += max(0, sequence_length - ref_end)

    #for i, (ref_start, ref_end) in enumerate(reference_intervals):
    #    if i == 0: # take care of bases before first reference interval
    #        TN += ref_start  # assuming genome starts at position 0
    #    elif i < len(reference_intervals) - 1:
    #        for pipe_start, pipe_end in pipeline_intervals:
    #            if pipe_start >= ref_start and pipe_end <= ref_end:
    #                TP += pipe_end-pipe_start + 1
    #                #pipeline_intervals.remove((pipe_start, pipe_end)) # since the pipeline interval is fully used, remove it
    #            elif pipe_start >= ref_start and pipe_start <= ref_end and pipe_end > ref_end:
    #                TP += ref_end - pipe_start + 1
    #                FP += pipe_end - ref_end
    #            elif pipe_start < ref_start and pipe_end <= ref_end and pipe_end >= ref_start:
    #                TP += pipe_end - ref_start + 1
    #                FP += ref_start - pipe_start
    #        TN += ref_start - previous_end - 1
    #        previous_end = ref_end
    #    elif i < len(reference_intervals) - 1:
    #        # add all bases that are not covered by pipeline intervals BUT are covered by reference as FN
    #        covered = 0
    #        for pipe_start, pipe_end in pipeline_intervals:
    #            overlap_start = max(ref_start, pipe_start)
    #            overlap_end = min(ref_end, pipe_end)
    #            if overlap_start <= overlap_end:
    #                covered += overlap_end - overlap_start + 1
    #        FN += (ref_end - ref_start + 1) - covered
    #    elif i == len(reference_intervals) - 1: #take care of remaining bases
    #        TN += max(0, sequence_length - ref_end)  
    return TP, FP, FN, TN

def confusion_matrix(reference_file, pipeline_file, sequence_length):
    # Create arrays to mark which bases are reference repeats and pipeline repeats
    ref_mask = np.zeros(sequence_length, dtype=bool)
    pipe_mask = np.zeros(sequence_length, dtype=bool)

    # Load reference intervals
    with open(reference_file) as ref_f:
        for line in ref_f:
            start, end = map(int, line.strip().split("\t"))
            ref_mask[start:end+1] = True  # mark reference repeats

    # Load pipeline intervals
    with open(pipeline_file) as pipe_f:
        for line in pipe_f:
            _, start, end, _, _, _ = line.strip().split("\t")
            start, end = int(start), int(end)
            pipe_mask[start:end+1] = True  # mark detected repeats

    # Now calculate TP, FP, FN, TN per base
    TP = np.sum(ref_mask & pipe_mask)
    FP = np.sum(~ref_mask & pipe_mask)
    FN = np.sum(ref_mask & ~pipe_mask)
    TN = np.sum(~ref_mask & ~pipe_mask)

    return TP, FP, FN, TN


def calculate_sensitivity(reference_file, pipeline_file, sequence_length, output_file):
    TP, FP, FN, TN = confusion_matrix(reference_file, pipeline_file, sequence_length)
    print(f"TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    print(f"Sensitivity: {sensitivity:.4f}")
    with open(output_file, "w") as out_f:
        out_f.write(f"TP\tFP\tFN\tTN\n{TP}\t{FP}\t{FN}\t{TN}\n")
        out_f.write(f"Sensitivity: {sensitivity:.4f}\n")
    return sensitivity


def calculate_specificity(reference_file, pipeline_file, sequence_length, output_file):
    TP, FP, FN, TN = confusion_matrix(reference_file, pipeline_file, sequence_length)
    print(f"TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    print(f"Specificity: {specificity:.4f}")
    with open(output_file, "a") as out_f:
        #out_f.write(f"TP\tFP\tFN\tTN\n{TP}\t{FP}\t{FN}\t{TN}\n")
        out_f.write(f"Sensitivity: {specificity:.4f}\n")
    return specificity

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Caulcuate sensitivity and specificity of pipeline (vs. YASS-detected repeats)")
    parser.add_argument("reference", help="Input identified reference regions file (usually YASS)")
    parser.add_argument("pipeline_output_file", help="File with pipeline-detected regions")
    parser.add_argument("analyzed_sequence", help="File with sequence")
    parser.add_argument("output_file", help="File to write the confusion matrix to")

    if len(sys.argv) == 1:
        print("\nWelcome to a tool used to calculate sensitivity and specificity of the created pipeline!")
        print("You need to provide input files and parameters.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    sequence = SeqIO.read(args.analyzed_sequence, "fasta")
    seq_length = len(str(sequence.seq))

    calculate_sensitivity (args.reference, args.pipeline_output_file, seq_length, args.output_file)
    calculate_specificity (args.reference, args.pipeline_output_file, seq_length, args.output_file)