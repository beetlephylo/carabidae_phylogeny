#!/usr/bin/env python3

from Bio import SeqIO
import csv
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Select longest sequences for PTP delimited species")
    parser.add_argument("-i", "--input", type=str, help="mPTP output txt file")
    parser.add_argument("-s", "--supermatrix", type=str, help="Fasta with the tree's supermatrix")
    parser.add_argument("-k", "--keep", type=str, help="File with list of taxa to keep (eg constraint, outgroup)")
    parser.add_argument("-n", "--named", action='store_true', help="Keep all taxa with binomial (string after last underscore all lower case))")
    parser.add_argument("-o", "--output", type=str, help="Output file with selected fasta IDs")

    return parser.parse_args()


def count_bases(input_fasta):
    count = {}
    records = SeqIO.parse(input_fasta, "fasta")
    x = 0
    for rec in records:
        x += 1
        seq = rec.seq.upper()
        gaps = seq.count('-') + seq.count('N') + seq.count('X') + seq.count('*')
        length = len(seq) - gaps
        count[rec.id] = length
    print(f'{x} taxa in supermatrix')
    return count


def process_mptp_file(mptp_file):
    species_lists = {}
    temp = []
    start_processing = False
    with open(mptp_file) as file:
        lines = file.readlines()
        for line in lines:
            if 'Number of delimited species' in line:
                print(line.strip())
            elif 'Species ' in line:
                species_no = line.strip().replace('Species ', '').replace(':', '')
                species_no = re.findall(r'\d+', line)[0]
                species_lists[species_no] = []
                start_processing = True
            else:
                if start_processing and line.strip() != '':
                    species_lists[species_no].append(line.strip())
    print(f'{len(list(species_lists.keys()))} delimited species groups in {mptp_file}')
    return species_lists


def find_longest(count, species_lists, keep_ids, named):
    keep = set()
    # Save IDs from keep list
    if keep_ids:
        file = open(keep_ids)
        lines = file.readlines()
        keep.add(line.strip() for line in lines)

    for species_no, species_list in species_lists.items():
        nt = 0
        for spec in species_list:
            # Save IDs with species/subspecies at end of id string.
            if named:
                parts = spec.rsplit('_', 1)
                try:
                    if re.fullmatch(r'[a-z]+', parts[1]):
                        keep.add(spec)
                except IndexError:
                    pass            
            # Find longest seequence for each species group
            if count[spec] > nt:
                nt = count[spec]
                longest = spec
        keep.add(longest)

    return keep

def main():
    args = parse_args()
    count = count_bases(args.supermatrix)
    species_lists = process_mptp_file(args.input)
    selected = find_longest(count, species_lists, args.keep, args.named)

    with open(args.output, 'w') as output:
        for species_id in selected:
            output.write(f'{species_id}\n')

    print(f'{len(selected)} taxa saved to {args.output}')
    
if __name__ == "__main__":
    main()