
gene <- "MEFV"

setwd("~/git/MORPHA/")
samples <- read.delim(paste0("results/", gene, "/", gene, "_protein_haplotypes_matched_to_igsr.tsv"), sep = "\t", header = TRUE, check.names = FALSE)

# FOR NOW, ONLY SELECT WT
samples <- subset(samples, Origin=="WT")

head(samples)


# /Users/joeri/git/MORPHA/results/MEFV/folded/MEFV_protein_variant_12_Repair_0_ST.fxout


library(dplyr)
library(readr)
library(stringr)
library(purrr)

# 1. Define the correct FoldX column names based on your output mapping
foldx_cols <- c(
  "Pdb_file", "Total_Energy", "Backbone_Hbond", "Sidechain_Hbond", 
  "Van_der_Waals", "Electrostatics", "Solvation_Polar", "Solvation_Hydrophobic", 
  "Van_der_Waals_clashes", "Entropy_Side_Chain", "Entropy_Main_Chain", 
  "Sloop_Entropy", "Mloop_Entropy", "Cis_Bond", "Torsional_Clash", 
  "Backbone_Clash", "Helix_Dipole", "Water_Bridge", "Disulfide", 
  "Electrostatic_Kon", "Partial_Covalent_Bonds", "Energy_Ionisation", 
  "Entropy_Complex", "Residue_Count"
)

# 2. List all the FoldX output files in your current working directory
# (Adjust path inside list.files() if they are in a different folder)
# 2. List all the FoldX output files with FULL paths
# Adding full.names = TRUE ensures read_tsv can actually find them!
target_path <- paste0("results/", gene, "/folded/")
fxout_files <- list.files(pattern = "MEFV_protein_variant_\\d+_Repair_0_ST\\.fxout",path = target_path, full.names = TRUE)

# 3. Read and combine all FoldX files into a single master data frame
foldx_data <- fxout_files %>%
  map_dfr(function(file_path) {
    # Read the tab-delimited file using its full path
    df <- read_tsv(file_path, col_names = foldx_cols, show_col_types = FALSE)
    
    # Extract just the variant identifier from the path/filename string
    # e.g., extracts "MEFV_protein_variant_1" out of "results/MEFV/folded/MEFV_protein_variant_1_Repair_0_ST.fxout"
    df <- df %>% 
      mutate(Variant = str_extract(file_path, "MEFV_protein_variant_\\d+"))
    
    return(df)
  })

# 4. Merge the FoldX data back into your original samples data frame
# This matches on the "Variant" column present in both datasets
updated_samples <- samples %>%
  left_join(foldx_data, by = "Variant")

updated_samples$Pdb_file <- NULL

# 5. View your newly updated data frame
head(updated_samples)

write.table(updated_samples, paste0("results/", gene, "/", gene, "_protein_haplotypes_quant.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

