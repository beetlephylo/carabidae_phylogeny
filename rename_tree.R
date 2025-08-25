library(getopt)
library(ape)
library(dplyr)
library(stringr)

# Define the command-line arguments
spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input tree file',
  'csv',      'c', 2, 'character', 'Metadata CSV. Specify custom labels with -l flag, or default assumes new names in first column, old names in second',
  'tips',     't', 2, 'character', 'Column name with original tip names (if not first column)',
  'label',    'l', 2, 'character', 'Custom labels: specify columns to use in labels: comma separated column names in order
                                    Default without this option is new names in first column, old names in second column',
  'output',   'o', 1, 'character', 'Output tree file',
  'renamed',  'r', 2, 'logical',   'CSV output with old and new tips names',
  'drop_old', 'd', 2, 'logical',   'Drop original tip names (default keep old name at start of new name)',
  'prefix',   'p', 2, 'character', 'Add specified string to the start of all labels'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# Read tree and strip quotes from tip names
read_tree <- function(input_tree) {
  nexus <- FALSE
  first_line <- readLines(input_tree, n = 1)
  # Check for the #NEXUS header
  if (grepl("#NEXUS", first_line, ignore.case = TRUE)) {
    nexus <- TRUE
    tree <- (read.nexus(input_tree))
  } else {
    tree <- (read.tree(input_tree))
  }
  cat("Tree has", length(tree$tip.label), "tips\n")
  tree$tip.label <- gsub("(^['\"]|['\"]$)", '', tree$tip.label)
  return(list(tree = tree, nexus = nexus))
}


# Read metadata and process
process_metadata <- function(metadata, tree, tips, label, drop_old, prefix) {
  df <- read.csv(metadata)
  df[df == ''] <- NA
  df[] <- lapply(df, as.character)

  # Get old_names column
  old_names_col <- if (is.null(tips)) names(df)[2] else tips
  df <- df %>% rename(old_names_col = !!sym(old_names_col)) %>%
    filter(old_names_col %in% tree$tip.label)
  
  # If --label is not provided, use the first column as the new name
  if (is.null(label)) {
    cat("Tree will be renamed with the", names(df[1]), "column from", metadata, '\n')
    new_names_col <- names(df)[1] 
    df <- df %>%
      rename(new_names_col = !!sym(new_names_col))
  } else {
    if (!is.null(tips)) {
      label <- str_replace(label, tips, 'old_names_col')
    }
    # Write custom labels
    cat("Tree will be renamed with the following column(s) from", metadata, ": ", label, '\n')
    column_names <- strsplit(opt$label, ",")[[1]]
    df <- df %>%
      rowwise() %>% 
      mutate(new_names_col = ifelse(is.null(drop_old),
        paste(c(old_names_col, c_across(all_of(column_names))), collapse = "_"), 
        paste(c(c_across(all_of(column_names))), collapse = "_"))) %>%
      ungroup()
    }
    # Add prefix
    if (!is.null(prefix)) {
      df$new_names_col <- paste0(prefix, "_", df$new_names_col)
    }
    df <- df %>%
      mutate(new_names_col = gsub("NA", "", new_names_col)) %>%
      mutate(new_names_col = gsub("_+$|\\s", "_", new_names_col)) %>%
      select(old_names_col, new_names_col) %>%
      unique()
  
  return(df)
}

# Rename tree
rename_tips <- function(tree, df) {
  matches <- match(tree$tip.label, df$old_names_col)
  tree$tip.label[!is.na(matches)] <- df$new_names_col[matches[!is.na(matches)]]
  return(tree)
}


# Write tree
write_tree <- function(tree, output, nexus) {
  # tree <- rename_tips(tree, df)
  if (nexus == FALSE) {
    write.tree(tree, file=output)
    cat("Tree written to", output, '\n')
  } else {
    write.nexus(tree, file=output)
    cat("Tree written to", output, '\n')
  }
}


main <- function(opt) {
  input <- read_tree(opt$input)
  df <- process_metadata(opt$csv, input$tree, opt$tips, opt$label, opt$drop_old, opt$prefix)
  tree <- rename_tips(input$tree, df)
  write_tree(tree, opt$output, input$nexus)
}


main(opt)