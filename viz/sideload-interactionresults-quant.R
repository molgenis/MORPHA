### sideload "quant" in world-plot-1KG-subpop.R script with this
# quick way to produce heatmap of results
load_interaction_energies <- function(file_a_path, file_b_path) {
  # Read tab-separated files
  df_a <- read.delim(file_a_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  df_b <- read.delim(file_b_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Extract clean Variant ID from File A's Pdb file paths
  # e.g., "./MEFV_protein_path_variant_845_Repair.pdb" -> "MEFV_protein_path_variant_845"
  df_a$Variant <- sub("_Repair\\.pdb$", "", basename(df_a$Pdb))
  
  # Merge File A's Interaction Energy into File B based on Variant
  merged_df <- merge(
    df_b[, c("Sample_name", "Haplotype", "Variant", "Origin", "VariantID")],
    df_a[, c("Variant", "Interaction Energy")],
    by = "Variant"
  )
  
  # Reorder columns to match File B's exact structure
  final_df <- merged_df[, c("Sample_name", "Haplotype", "Variant", "Origin", "VariantID", "Interaction Energy")]
  
  return(final_df)
}

# Usage:
# result <- load_interaction_energies("path/to/fileA.tsv", "path/to/fileB.tsv")
quant <- load_interaction_energies(
  #"results/MEFV/complex_14-3-3_MEFV_Foldx.tsv",
  "results/MEFV/complex_ASC_MEFV_Foldx.tsv",
  paste0("results/", gene, "/", gene, "_protein_haplotypes_quant.tsv")
)
names(quant)[names(quant) == "Interaction Energy"] <- "Total_Energy"
