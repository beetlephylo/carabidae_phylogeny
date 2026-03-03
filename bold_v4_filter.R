# Extract sequences from BOLD download. Input existing metadata/sequence download from BOLD v4.

library('tidyverse')
library('getopt')


spec <- matrix(c(
  'tsv',        'c', 2, 'character', 'Path to existing BOLD v4 TSV file to be filtered',
  'genbank',    'g', 2, 'logical',   'Remove sequences also on GenBank',
  'barcode',    'b', 2, 'character', 'Save only barcodes, delete other genes',
  'raw',        'r', 2, 'logical',   'Save raw metadata, before any processing.',
  'bins',       'f', 2, 'logical',   'Keep only one sequence per BIN',
  'species',    's', 2, 'logical',   'Keep only one sequence per named species'
), byrow = T, ncol = 5)

opt <- getopt(spec)


# Read in data
meta <- read.csv(opt$tsv, header = TRUE, sep = "\t", quote = "")


# Filter out GenBank sequences

if (!is.null(opt$genbank)) {
  # Filter for non-empty and non-missing genbank_accession
  f_meta <- meta %>% filter(is.na(genbank_accession) | genbank_accession == "")
  print(paste(nrow(f_meta), 'records remaining after removing those with GenBank accession'))
} else {
  f_meta <- meta
}

# Keep only COI sequences
if ( !is.null(opt$barcode)) {
  f_meta <- f_meta %>% filter(markercode=="COI-5P")
  print(paste(nrow(f_meta), 'records with COI-5P sequences'))
}

# Remove records without sequences or named genes
f_meta <- f_meta %>% filter(markercode!="")
f_meta <- f_meta %>% filter(nucleotides!="")

# Remove gaps from sequences
f_meta$nucleotides <- gsub("-","",as.character(f_meta$nucleotides))
# Get sequence lengths
f_meta$sequence_length <- nchar(f_meta$nucleotides)

# Keep longest sequence for each bin, for each gene
if ( !is.null(opt$bins)) {
  na_or_empty <- f_meta %>% filter(is.na(bin_uri) | bin_uri == "") 
  unique_bin_uri <- f_meta %>% 
    arrange(bin_uri, markercode, desc(sequence_length)) %>% 
    distinct(bin_uri, .keep_all = TRUE) 
  f_meta <- bind_rows(na_or_empty, unique_bin_uri)
  cat(paste(nrow(unique_bin_uri), 'unique BINs\nSaved one sequence for each BIN\n'))
  cat(paste(nrow(na_or_empty), 'records without a BIN\n'))
  cat(paste('Saved', nrow(f_meta), 'records total\n'))
}

# Keep longest sequence for each bin, for each gene
if ( !is.null(opt$species)) {
  na_or_empty <- f_meta %>% filter(is.na(species_name) | species_name == "") 
  unique_species <- f_meta %>% 
    arrange(species_name, markercode, desc(sequence_length)) %>% 
    distinct(species_name, .keep_all = TRUE) 
  f_meta <- bind_rows(na_or_empty, unique_species)
  cat(paste(nrow(unique_species), 'unique speciess\nSaved one sequence for each species\n'))
  cat(paste(nrow(na_or_empty), 'records without a species name\n'))
  cat(paste('Saved', nrow(f_meta), 'records total\n'))
}

#Edit dataframe to match lab metadata
empty <- c('lab_id', 'subgenus', 'subtribe',	'tribe', 'superfamily',	'infraorder',	'suborder',	'genbank_accession')
f_meta[ , empty] <- ''
csv <- f_meta %>% select(genbank_accession,	processid, bin_uri,	lab_id, suborder,	infraorder, superfamily, family_name, 
                         subfamily_name, tribe, subtribe, genus_name, subgenus, species_name, subspecies_name, country,	lat,	lon)

new_names <- c("genbank_accession",	"bold_id",	"bold_bin",	"lab_id",	"suborder", "infraorder",	"superfamily",	
               "family", "subfamily",	"tribe", "subtribe",	"genus",	"subgenus",	"species", "subspecies", "country",	"latitude",	"longitude")
names(csv) <- new_names

csv <- csv %>% unique()


# Write metadata to CSV
write.csv(csv, 'metadata.csv', row.names = FALSE)
print('Metadata saved to metadata.csv')

# Get list of gene names
genes <-c(unique(f_meta$markercode))
print('Genes found:')
print(genes)

# Write fastas
for (gene in genes) {
  df <- f_meta %>% filter(markercode==gene)
  file_out <- file(paste(gene, ".fasta", sep = ''), "w")
  for (i in 1:nrow(df)) {
    cat(">", df$processid[i], "\n", file = file_out, sep = '')
    cat(df$nucleotides[i], "\n", file = file_out, sep = '')
  }
  close(file_out)
  cat(length(df$nucleotides), ' sequences writen to ', gene, '.fasta\n', sep = '')
}