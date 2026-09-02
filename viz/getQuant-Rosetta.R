# requires:
# results/<gene>/<gene>_protein_haplotypes_matched_to_igsr.tsv -> produced by viz/merge_haplo_with_igsrmd.R
# results/<gene>/folded/ <gene>_protein_path_variant_*_Repair.sc files -> produced in this case by Rosetta FastRelax (on FoldX repaired PDBs)

library(dplyr)
library(readr)
library(stringr)
library(purrr)

gene <- "TTR"

# load the sample data
setwd("~/git/MORPHA/")
samples <- read.delim(
  paste0("results/", gene, "/", gene, "_protein_haplotypes_matched_to_igsr.tsv"), 
  sep = "\t", 
  header = TRUE, 
  check.names = FALSE
)

# locate Rosetta score files (.sc, .fasc, or text score files)
target_path <- paste0("results/", gene, "/folded/")
rosetta_files <- list.files(
  pattern = ".*_variant_\\d+.*\\.sc$", # Adjust pattern to your file extension (.sc / .fasc)
  path = target_path, 
  full.names = TRUE
)

# read Rosetta score files and extract the minimum total_score per variant
rosetta_data <- rosetta_files %>%
  map_dfr(function(file_path) {
    # Read score file lines
    lines <- read_lines(file_path)
    
    # Identify the header line starting with "SCORE:" containing column names
    header_line <- lines[str_detect(lines, "^SCORE:\\s+total_score")]
    if (length(header_line) == 0) return(NULL)
    
    # Extract column headers
    cols <- unlist(str_split(str_trim(header_line), "\\s+"))
    cols <- cols[cols != "SCORE:"] # Remove the "SCORE:" label prefix
    
    # Filter data rows containing numerical scores
    data_lines <- lines[str_detect(lines, "^SCORE:") & !str_detect(lines, "total_score")]
    
    # Convert data rows to a clean data frame
    df <- read_table(
      paste(data_lines, collapse = "\n"), 
      col_names = FALSE, 
      show_col_types = FALSE
    )
    
    # Drop the leading 'SCORE:' column from values
    df <- df %>% select(-1)
    colnames(df) <- cols
    
    # Ensure numerical types for total_score
    df <- df %>% mutate(total_score = as.numeric(total_score))
    
    # Extract variant ID (e.g., "MEFV_protein_variant_1" or "MEFV_protein_path_variant_1")
    # Matches description string like "TP53_protein_path_variant_1_Repair_0006"
    df <- df %>%
      mutate(Variant = str_extract(description, paste0(gene, "_protein_(path_)?variant_\\d+")))
    
    return(df)
  }) %>%
  # Group by variant and select the structure with minimum energy
  group_by(Variant) %>%
  slice_min(order_by = total_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# select only total_score and Variant
rosetta_data <- rosetta_data %>% select(total_score, Variant)

# merge Rosetta minimum score back into samples
updated_samples <- samples %>% left_join(rosetta_data, by = "Variant")

updated_samples <- updated_samples %>% rename("Total_Energy" = "total_score")

# save the output
write.table(
  updated_samples, 
  paste0("results/", gene, "/", gene, "_protein_haplotypes_quant_rosetta.tsv"), 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)


