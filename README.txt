Code and scripts used for 'Beetles, Barcodes, and Big Data: A Deep Dive into Harpalinae Phylogeny'.

bold_v4_filter.R
Filters output from BOLD v4 and saves a fasta file for each gene.

bold_v4_search.py
Retrieve sequences and metadata with the BOLD v4 API.

bold_v5_filter.R
Filters output from BOLD v5 and saves a fasta file for each gene.

bold_v5_search.py
Retrieve sequences and metadata using the BOLD v5 API.

clean_rna.py
Aligns RNA sequences using MAFFT, and filters out those that do not meet the similarity threshold.

fasta_length_filter.py
Remove sequences shorter than specified length threshold.

filter_fasta.py
Filter a fasta file to separate those with IDs specificed in an input file.

partitions.py
Change catfasta2phyml.pl partition file to RAxML format.

process_blast.py
Filter BLAST output file to get coordinates of chosen sequences. Compatible with BLAST outfmt '6 staxids sacc sstart send'.

ptp_filter_output.py
Select representative sequence for each delimited group from mPTP output file.

rename_fasta.py
Change fasta IDs using CSV file.

rename_tree.R
Change nexus or newick tip IDs (eg add taxonomy) from CSV file.

ry_code.py
Change nucleotide fasta file to purine/pyrimidine binary data (RY or 01).

supermatrix_count.py
Produce CSV of sequence length for each taxon and each partition in a fasta file (partition file in catfasta2phyml.pl format).