#!/usr/bin/env python3
"""
Check whether custom repeat regions overlap with CDS features
in a GFF genome annotation.

Usage:
    python analyze_repeat_coding_overlap.py --bed repeats.bed --gff genome.gff
"""

import argparse
from pybedtools import BedTool
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="Check overlap between repeats and CDS regions in GFF.")
    parser.add_argument("bed", help="Input BED file with repeat regions")
    parser.add_argument("gff",help="Genome annotation GFF file")
    parser.add_argument("--prefix", default="repeats", help="Prefix for output files")
    parser.add_argument("--chromosome", default=None, help="Specific chromosome to analyze (use if you want to limit analysis)")
    args = parser.parse_args()

    # --- Load input files ---
    print("Loading BED and GFF files...")
    repeats = BedTool(args.bed)
    annotation = BedTool(args.gff)

    # --- Filter by chromosome if specified ---
    if args.chromosome:
        print(f"Filtering repeats and annotation for chromosome: {args.chromosome}...")
        repeats = repeats.filter(lambda x: x.chrom == args.chromosome).saveas()
        annotation = annotation.filter(lambda x: x.chrom == args.chromosome).saveas()

    #print(annotation)
    gene_dict = defaultdict(list)
    for line in annotation:
        _, _, info_type, start, stop, _, _, _, info = str(line).strip().split("\t")
        if info_type == "transcript":
            gene_id = info.strip().split(";")[0]
            transcript_id =  info.strip().split(";")[1]
            gene_dict[gene_id].append(transcript_id)
    print(gene_dict)
    many_transcripts = 0
    for value_list in gene_dict.values():
        if len(value_list) > 1:
            many_transcripts += 1
    print(many_transcripts)
    #print(info_type)


    # --- Extract CDS features only ---
    print("Extracting CDS features from GFF...")
    cds = annotation.filter(lambda x: x[2] == "transcript")
    #.saveas(f"{args.prefix}_CDS_only.gff")    
    # --- Find repeats that overlap transcript intervals ---
    overlaps = repeats.intersect(cds, wa=True, wb=True)

    # identify unique repeat intervals that overlap
    overlapping_repeats = set()
    for x in overlaps:
        # unique key representing the interval
        key = f"{x.chrom}:{x.start}-{x.end}"
        overlapping_repeats.add(key)

    # --- Create list of non-overlapping repeats ---
    noncoding_repeats = []
    for r in repeats:
        key = f"{r.chrom}:{r.start}-{r.end}"
        if key not in overlapping_repeats:
            noncoding_repeats.append(r)

    # save both
    BedTool(overlaps).saveas(f"{args.prefix}_in_transcript.bed")
    BedTool(noncoding_repeats).saveas(f"{args.prefix}_non_transcript.bed")

    # --- Find overlaps ---
    #print("Finding overlaps between repeats and CDS...")
    #overlaps = repeats.intersect(cds, wa=True, wb=True)
    #overlaps.saveas(f"{args.prefix}_in_CDS.bed")
#
    ## --- Find repeats NOT overlapping CDS ---
    #print("Finding non-coding (non-CDS) repeats...")
    #noncoding = repeats.intersect(cds, v=True)
    #noncoding.saveas(f"{args.prefix}_nonCDS.bed")

    # --- Summary statistics ---
    #n_cds = len(overlaps)
    #n_non = len(noncoding)
    n_total = len(repeats)
    n_overlap = len(overlapping_repeats)
    n_non = len(noncoding_repeats)

    print("\n=== SUMMARY ===")
    print(f"Total repeats:            {n_total}")
    print(f"Repeats overlapping TX:   {n_overlap}")
    print(f"Repeats non-overlapping:  {n_non}")
    #print(f"Percent overlapping:      {100 * n_overlap / n_total:.2f}%")


    #print(f"Total repeats:         {n_total}")
    #print(f"Repeats in CDS:        {n_cds}")
    #print(f"Repeats in non-CDS:    {n_non}")
    if n_total > 0:
        print(f"Percent overlapping:        {100 * n_overlap / n_total:.2f}%")

    print(f"\nOutput files:")
    print(f"  {args.prefix}_CDS_only.gff")
    print(f"  {args.prefix}_in_CDS.bed")
    print(f"  {args.prefix}_nonCDS.bed")
    print("\nDone!")

if __name__ == "__main__":
    main()
