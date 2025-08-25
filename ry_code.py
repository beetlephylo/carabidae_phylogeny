from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re

def parse_args():

    parser = argparse.ArgumentParser(description="RY coding")
    parser.add_argument("-i", "--input", type=str, help="Input fasta")
    parser.add_argument("-o", "--output", type=str, help="Output fasta")
    parser.add_argument("-r", "--ry", action='store_true', help="R/Y output instead of default numeric 1/0")
    parser.add_argument("-t", "--third", action='store_true', help="RY code only third codon positions")
    parser.add_argument("-d", "--data", choices=['nuc', 'mito'], default='mito', help="Specify nuclear or invertebrate mitochondiral DNA if translating only 3rd codon positions. Default is mitochondrial.")

    return parser.parse_args()

# 1/0 code
def replace_numeric(base):
    base = base.upper()
    return (
        '1' if base in {'A', 'G', 'R'} else
        '0' if base in {'C', 'T', 'Y'} else
        '-'
    )


# R/Y code
def replace_ry(base):
    base = base.upper()
    return (
        'R' if base in {'A', 'G', 'R'} else
        'Y' if base in {'C', 'T', 'Y'} else
        '-'
    )


def find_reading_frame(records, data):
    stop_codons = ['TAA', 'TAG'] if data=='mito' else ['TAA', 'TAG', 'TGA']    # Add extra stop codon for nuclear DNA
    frame_totals = [0, 0, 0]
    for rec in records:
        # Count stop codons in each frame
        for f in range(3):
            x = 0
            codons = [rec.seq[i:i+3] for i in range(f, len(rec.seq), 3)]
            for codon in codons:
                if codon in stop_codons:
                    x += 1
            frame_totals[f] += x
    # Choose reading frame with fewest stop codons
    frame = frame_totals.index(min(frame_totals))
    print('Reading frame:', frame)
    length = len(records[0].seq)
    return frame, length


def code_third_codon_position(records, data, ry=False):
    # Get reading frame
    records = list(records)
    frame, length = find_reading_frame(records, data)
    for rec in records:
        seq = list(Seq(rec.seq).upper())
        # Replace third base of each codon
        for i in range(frame + 2, len(seq), 3):
            base = seq[i]
            if ry:
                seq[i] = replace_ry(base)
            else:
                seq[i] = replace_numeric(base)

        rec.seq = Seq(''.join(seq))
    return records

def code_all_bases(records, ry=False):
    records = list(records)
    for rec in records:
        seq = list(Seq(rec.seq).upper())
        for i in range(len(seq)):
            base = seq[i]
            if ry:
                seq[i] = replace_ry(base)
            else:
                seq[i] = replace_numeric(base)
        rec.seq = Seq(''.join(seq))
    return records


def main():
    args = parse_args()
    records = SeqIO.parse(args.input, "fasta")
    if args.third:
        records = code_third_codon_position(records, args.data, args.ry)
    else:
        records = code_all_bases(records, args.ry)
    with open(args.output, "w") as output:
        SeqIO.write(records, output, "fasta")
        print(f'{len(list(records))} records written to {args.output}')

if __name__ == "__main__":
    main()