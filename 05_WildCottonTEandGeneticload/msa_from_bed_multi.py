#!/usr/bin/env python3
"""
Extract MSA columns per chromosome based on BED regions.

Usage:
    python msa_from_bed_multi.py --bed sites.bed --msa_dir msa_folder/ --out_dir sub_alignment/
"""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Extract MSA columns from BED regions for multiple chromosomes")
    parser.add_argument("--bed", required=True, help="BED file with regions (0-based, end-exclusive)")
    parser.add_argument("--msa_dir", required=True, help="Directory containing per-chromosome MSA FASTA files (named chrom.fasta)")
    parser.add_argument("--out_dir", required=True, help="Directory to save sub-alignments")
    return parser.parse_args()

def read_bed(bed_file):
    """Return dictionary: chrom -> list of (start, end) tuples"""
    bed_dict = {}
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start, end = int(parts[1]), int(parts[2])
            bed_dict.setdefault(chrom, []).append((start, end))
    return bed_dict

def build_ref_pos_map(ref_seq):
    """Map genomic positions to MSA column indices (skip gaps)."""
    pos_to_col = {}
    ref_pos = 0
    for col, base in enumerate(ref_seq):
        if base != '-':
            pos_to_col[ref_pos] = col
            ref_pos += 1
    return pos_to_col

def extract_sub_alignment(alignment, regions):
    ref_seq = alignment[0].seq
    pos_to_col = build_ref_pos_map(ref_seq)

    columns_to_extract = []
    for start, end in regions:
        for pos in range(start, end):
            if pos in pos_to_col:
                columns_to_extract.append(pos_to_col[pos])

    sub_alignment = MultipleSeqAlignment([])
    for record in alignment:
        new_seq = "".join(record.seq[col] for col in columns_to_extract)
        sub_alignment.append(SeqRecord(Seq(new_seq), id=record.id, description=""))
    return sub_alignment

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    bed_dict = read_bed(args.bed)

    for chrom, regions in bed_dict.items():
        msa_file = os.path.join(args.msa_dir, f"{chrom}.clean.fa")
        if not os.path.exists(msa_file):
            print(f"[Warning] MSA file for {chrom} not found, skipping...")
            continue

        alignment = AlignIO.read(msa_file, "fasta")
        sub_alignment = extract_sub_alignment(alignment, regions)
        out_file = os.path.join(args.out_dir, f"{chrom}_sub_alignment.fasta")
        AlignIO.write(sub_alignment, out_file, "fasta")
        print(f"{chrom}: {len(sub_alignment[0])} columns extracted -> {out_file}")

if __name__ == "__main__":
    main()
