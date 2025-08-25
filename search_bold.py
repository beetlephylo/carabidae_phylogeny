# Script to download sequences and associated metadata from BOLD v4

import argparse
import requests

def parse_args():
    parser = argparse.ArgumentParser(description="Download barcodes and metadata from BOLD Systems.")
    parser.add_argument("-t", "--taxon", help="Specify taxon (single or comma delimited list)")
    parser.add_argument("-i", "--ids", help="File with list of BOLD IDs, one per line")
    parser.add_argument("-o", "--output", help="Output TSV file path")

    return parser.parse_args()


def search_bold_taxon(taxa):
    # API URL
    taxa = '|'.join(taxa)
    combined_url = f"http://www.boldsystems.org/index.php/API_Public/combined?taxon={taxa}&format=tsv"
    # Download metadata
    print(f"Downloading data for taxon: {taxa}")
    response = requests.get(combined_url)
    tsv_data = response.text
    return tsv_data


def search_bold_ids(ids, chunk_size=100):
    # Split the list of IDs into chunks
    chunks = [ids[i:i + chunk_size] for i in range(0, len(ids), chunk_size)]

    # Initialize an empty list to store the results
    results = []

    # Make API calls for each chunk
    for i, chunk in enumerate(chunks):
        # API URL
        combined_url = f"http://www.boldsystems.org/index.php/API_Public/combined?ids={'|'.join(chunk)}&format=tsv"
        print(combined_url)
        # Download metadata
        print(f"Downloading data for IDs in chunk {i}")
        response = requests.get(combined_url)
        tsv_data = response.text

        # If this is not the first chunk, skip the header line
        if i > 0:
            tsv_data = '\n'.join(tsv_data.split('\n')[1:])

        results.append(tsv_data)

    # Join the results together
    combined_results = '\n'.join(results)

    return combined_results


def main():
    args = parse_args()

    if args.taxon:
        taxa = args.taxon.split(',')
        data = search_bold_taxon(taxa)
        
    elif args.ids:
        with open(args.ids, 'r') as file: 
            ids = file.read().splitlines()
            #ids = '|'.join(ids)
        data = search_bold_ids(ids)

    # Write to file
    if args.output:
        output = args.output
    else:
        output = f"{args.taxon}_metadata.tsv"

    with open(output, "w") as file:
        file.write(data)

    if args.output:
        print(f"Metadata saved to {args.output}")
    else:
        print(f"Metadata saved to {args.taxon}_metadata.tsv")

if __name__ == "__main__":
    main()
