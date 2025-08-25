#!/usr/bin/env python3

from Bio import SeqIO, Phylo
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Filter fasta file from list of IDs")
    parser.add_argument("-i", "--input", type=str, help="Input fasta to be filtered")
    parser.add_argument("-f", "--found", type=str, help="Output file for sequences found in search filter")
    parser.add_argument("-n", "--notfound", type=str, help="Output file for input sequences not present in search filter")
    parser.add_argument("-s", "--search", type=str, help="File with IDs to search for. Specify file type with -t.")

    parser.add_argument('-t', '--type', type=str, choices=['list', 'fasta', 'tree'],
                        help="File with IDs to searhc for: list with one ID per line (default), fasta file or newick tree file")

    return parser.parse_args()


def read_id_file(id_file, format=list):
    filter = []
    with open(id_file) as infile:
        if format == 'list':
            lines = infile.readlines()
            for line in lines:
                if line != '\n':
                    filter.append(line.strip())

        elif format == 'fasta':
            recs = SeqIO.parse(args.search, "fasta")
            for rec in recs:
                filter.append(rec.id)

        elif format == 'tree':
            tree = Phylo.read(infile, "newick")
            taxa = tree.get_terminals()
            filter = [taxon.name for taxon in taxa]
    filter = [re.sub(r';frame=\d', '', f) for f in filter]
    print(f'{len(filter)} IDs in {id_file}')
    return filter


def filter_fasta(id_file, filter, input_fasta, found, not_found):
    records = SeqIO.parse(input_fasta, "fasta")
    print(f'{len(list(records))} records in {input_fasta}')
    frametag = re.compile(r';frame=\d')

    # For each record, check if it is present in filter file (without frame tag)
    if found:
        found_records = (r for r in SeqIO.parse(input_fasta, "fasta") if frametag.sub('', r.id) in filter)
        count = SeqIO.write(found_records, found, "fasta")
        print(f"{count} of {len(filter)} records from {id_file} found in {input_fasta}, saved to {found}")
    if not_found:
        not_found_records = (r for r in SeqIO.parse(input_fasta, "fasta") if frametag.sub('', r.id) not in filter)
        count = SeqIO.write(not_found_records, not_found, "fasta")
        print(f"{count} records from {input_fasta} not present in {id_file}, saved to {not_found}")


def main():
    args = parse_args()
    filter = read_id_file(args.search, args.type)
    filter_fasta(args.search, filter, args.input, found=args.found, not_found=args.notfound)
    

if __name__ == "__main__":
    main()