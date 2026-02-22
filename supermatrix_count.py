#!/usr/bin/env python3

import argparse
import csv
from Bio import SeqIO

# Argument parser
parser = argparse.ArgumentParser(description="Count number of genes and number of bases per gene in supermatrix")
parser.add_argument("-f", "--fasta", type=str, help="Input supermatrix fasta")
parser.add_argument("-p", "--partitions", type=str, help="catfasta2phyml output partition file")
parser.add_argument("-o", "--output", type=str, help="Output csv file")

args = parser.parse_args()

genes_dict = {"12S":    {"code": "A"},
              "16S":    {"code": "B"},
              "18S":    {"code": "C"},
              "28S":    {"code": "D"},
              "ATP6":   {"code": "E"},
              "ATP8":   {"code": "F"},
              "COX1a":  {"code": "G"},
              "COX1b":  {"code": "H"},
              "COX2":   {"code": "I"},
              "COX3":   {"code": "J"},
              "CYTB":   {"code": "K"},
              "ND1":    {"code": "L"},
              "ND2":    {"code": "M"},
              "ND3":    {"code": "N"},
              "ND4.":   {"code": "O"},
              "ND4L.":  {"code": "P"},
              "ND5":    {"code": "Q"},
              "ND6":    {"code": "R"},
              "AK":     {"code": "S"},
              "CAD":    {"code": "T"},
              "EF1A":   {"code": "U"},
              "Wg":     {"code": "V"}}


genes_list = ['12S', '16S', '18S', '28S', 'ATP6', 'ATP8', 'COX1a', 'COX1b', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4.', 'ND4L.', 'ND5', 'ND6', 'AK', 'CAD', 'EF1A', 'Wg']

# Read partition file
if args.partitions:
    with open(args.partitions) as file:
        parts = {}
        genes = []
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            part, crds = line.split(' = ')
            # Get partition location and length
            crds = crds.split('-')
            start = int(crds[0]) - 1
            end = int(crds[1])
            length = end - start
            # Get partition name
            part = part.split('/')[-1]
            # Add to dict list for each gene [start pos, end pos, length]        
            parts[part] = [start, end, length]
            genes.append(part)
        others = [part for part in genes if all(gene not in part for gene in genes_list)]

    print(f'{len(genes)} partitions in partition file')
else:
    print('No partition file provided')


# Write CSV metadata file
with open(args.output, "w") as file:
    writer = csv.writer(file)
    records = SeqIO.parse(args.fasta, "fasta")
    x = 0
    rows = []
    # Without partition file
    if not args.partitions:
        writer.writerow(["taxon", "total_nucleotides"])
        for rec in records:
            x += 1
            gaps = rec.seq.count('-') + rec.seq.count('N') + rec.seq.count('X')
            row = [rec.id, len(rec.seq) - gaps]
            rows.append(row)
    # With partition file
    else:
        genes_header = []
        for g in genes_list:
            for gene in genes:
                if g in gene:
                    genes_header.append(gene)

        writer.writerow(["gene_presence", "taxon", "total_genes", "total_nucleotides"] + genes_header + others)
        for rec in records:
            x += 1
            row = [rec.id, '', '']
            total = 0
            g = 0
            prec = [] # gene prescence/absence annotation for tree
            for gen in genes_list:
                for gene in genes:
                    if gen in gene:
                        # Count gaps in each gene
                        rec.seq = rec.seq.upper()
                        gaps = rec.seq.count('-', parts[gene][0], parts[gene][1]) + rec.seq.count('N', parts[gene][0], parts[gene][1]) + rec.seq.count('X', parts[gene][0], parts[gene][1])
                        length = parts[gene][2] - gaps
                        if length > 0:
                            g += 1
                            prec.append(genes_dict[gen]['code'])
                        else:
                            prec.append('-')
                        row.append(length)
                        total = total + length
            if gene in others:
                # Count gaps in each gene
                rec.seq = rec.seq.upper()
                gaps = rec.seq.count('-', parts[gene][0], parts[gene][1]) + rec.seq.count('N', parts[gene][0], parts[gene][1]) + rec.seq.count('X', parts[gene][0], parts[gene][1])
                length = parts[gene][2] - gaps
                if length > 0:
                    g += 1
                    prec.append(genes_dict[gen]['code'])
                else:
                    prec.append('-')
                row.append(length)
                total = total + length
            row[1] = g
            row[2] = total
            row = [''.join(prec)] + row
            rows.append(row)
    for row in rows:
        writer.writerow(row)
print(f'{x} taxa written to {args.output}')

