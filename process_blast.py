#!/usr/bin/env python3

import argparse
from Bio import Entrez, SeqIO
import time

def parse_args():
    # Paths to the files and directories
    parser = argparse.ArgumentParser(description="Filter BLAST file to get max and min sequence position from hits."
                                    "Choose to keep longest sequence for each NCBI TaxID, or for each genbank accession." 
                                    "Compatible with BLAST outfmt '6 staxids sacc sstart send'"
                                    "Result can be used to extract full sequences from BLAST database.")
    parser.add_argument("-i", "--input", type=str, help="BLAST output file")
    parser.add_argument("-o", "--output", type=str, help="Output file with match coordinates")
    parser.add_argument("-l", "--longest", action='store_true', help="Only keep longest sequence for each txid")
    parser.add_argument("-e", "--email", type=str, help="Your email address for NCBI")

    return parser.parse_args()

def process_blast_results(input):
    blast_hits = {}

    with open(input, "r") as file:
        data = file.readlines()
        for txid, gbid, start, end in (line.strip().split("\t") for line in data):
            start = int(start)
            end = int(end)
                    
            # Check for multiple TXIDs
            if ';' in txid:
                txid = search_genbank(gbid)

            # Swap the values if the alignment is on the reverse strand
            if start > end:
                start, end = end, start

            # Get match coordinates. Add to other matches for same record.    
            if txid in blast_hits:
                if gbid in blast_hits[txid]:
                    # Check for too long matches (errors from genomes)
                    if max([end, blast_hits[txid][gbid]['end']]) - min([start, blast_hits[txid][gbid]['start']]) < 4000:
                        # Check if start is lower than current start
                        if start < int(blast_hits[txid][gbid]['start']):
                            blast_hits[txid][gbid]['start'] = int(start)
                        # Check if end is higher than current end
                        if end > int(blast_hits[txid][gbid]['end']):
                            blast_hits[txid][gbid]['end'] = int(end)
                    else:
                        # If max coordinated of current and previous match are > 4000bp, save as separate sequence
                        # Merge with current _2 sequence if possible
                        if f'{gbid}_2' in blast_hits[txid]:
                            if max([end, blast_hits[txid][f'{gbid}_2']['end']]) - min([start, blast_hits[txid][gbid]['start']]) < 4000:
                                # Find first start and last end match for each gbid
                                if start < int(blast_hits[txid][f'{gbid}_2']['start']):
                                    blast_hits[txid][f'{gbid}_2']['start'] = int(start)
                                if end > int(blast_hits[txid][f'{gbid}_2']['end']):
                                    blast_hits[txid][f'{gbid}_2']['end'] = int(end)
                        else:
                            # Save to dict
                            blast_hits[txid][f'{gbid}_2'] = {'start': int(start), 'end': int(end)}
                else:
                    blast_hits[txid][gbid] = {'start': int(start), 'end': int(end)}
            else:
                blast_hits[txid] = {gbid: {'start': int(start), 'end': int(end)}}

    return blast_hits


# Search GenBank, wait and try again if 'too many requests' error
def search_genbank(ids, retries=60, delay=40):
    for attempt in range(0, retries):
        try:
            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
            results = SeqIO.parse(handle, "gb")
            try:
                for r in results:
                    if "db_xref" in r.features[0].qualifiers:
                        db_xref = r.features[0].qualifiers["db_xref"]
                        for ref in db_xref:
                            if "taxon" in ref:  # Get NCBI taxon, rather than BOLD cross ref
                                txid = "".join(filter(str.isdigit, ref))  # Extract numbers from NCBI taxon value
                                return txid
                    else:
                        print(f'Error getting TXID for {gbid}')
            except TypeError:
                print(f'Error getting TXID for {gbid}')
                continue

        except Entrez.HTTPError:
            print(f"HTTP error fetching records; retry in {delay} seconds")
            time.sleep(delay)
    print(f"Failed to retrieve records after {retries} attempts.")
    return ids

def keep_longest(blast_hits):
    longest = {}
    for txid, gbid in blast_hits.items():
        max_len = 200
        max_acc = ''
        for gb, rec in gbid.items():
            # Get longest match for each TXID
            rec_len = int(rec['end']) - int(rec['start'])
            if rec_len > max_len:
                max_len = rec_len
                max_acc = gb
        longest[txid] = {gbid: {'start': int(rec['start']), 'end': int(rec['end'])}}
    return longest



def main():
    args = parse_args()
    Entrez.email = args.email
    blast_hits = process_blast_results(args.input)
    if args.longest:
        blast_hits = keep_longest(blast_hits)

    with open(args.output,'w') as output:
        x = 0
        for txid, gbid in blast_hits.items():
            for gb, rec in gbid.items():
                gb = gb.replace('_2', '')
                rec_len = int(rec['end']) - int(rec['start'])
                if rec_len > 200:
                    output.write(f'{gb}\t{rec["start"]}\t{rec["end"]}\n')
                    x += 1

    print(f'{x} records written to {args.output}')


if __name__ == "__main__":
    main()





