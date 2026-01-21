# carabidae_phylogeny

Code and scripts used for 'Beetles, Barcodes, and Big Data: A Deep Dive into Harpalinae Phylogeny'.

bold_filter.R
Filters output from BOLD v4 and saves a fasta file for each gene.

carabidae_pipeline
Outlines workflow used in 'Beetles, Barcodes, and Big Data: A Deep Dive into Harpalinae Phylogeny'.

filter_fasta.py
Filter a fasta file to separate those with IDs specificed in an input file.

merge_renamed_fasta_csv.R
Join rename_fasta.py output files to show gene representation for each fasta ID.

partitions.py
Change catfasta2phyml.pl partition file to RAxML format.

process_blast.py
Filter BLAST output file to get coordinates of chosen sequences. Compatible with BLAST outfmt '6 staxids sacc sstart send'.

ptp_filter_output.py
Select representative sequence for each delimited group from mPTP output file.

rename_fasta.py
Change fasta IDs, from CSV file.

rename_tree.R
Change nexus or newick tip IDs (eg add taxonomy) from CSV file.

ry_code.py
Change nucleotide fasta file to purine/pyrimidine binary data (RY or 01).

search_bold.py
Retrieve sequences and metadata from BOLD v4 API.

supermatrix_count.py
Get CSV of sequence length for each taxon and each partition in a fasta file (partition file in catfasta2phyml.pl format).