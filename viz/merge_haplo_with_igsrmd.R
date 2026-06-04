
gene <- "MEFV"

setwd("~/git/MORPHA/")

# read sample names and and expand to haplotypes
samples_df <- read.delim("viz/igsr_samples.tsv", sep = "\t", header = TRUE, check.names = FALSE)
expanded_list <- lapply(1:nrow(samples_df), function(i) {
  row <- samples_df[i, ]
  data.frame(
    Sample_name = row$`Sample name`,
    Haplotype = c(paste0(row$`Sample name`, ".h1"), paste0(row$`Sample name`, ".h2")),
    stringsAsFactors = FALSE
  )
})
igsr_samples_expanded <- do.call(rbind, expanded_list)
dim(igsr_samples_expanded)
head(igsr_samples_expanded)

# --- helper Function to Parse FASTA files ---
parse_fasta <- function(file_path, variant_type) {
  fasta_lines <- readLines(file_path)
  header_lines <- fasta_lines[str_detect(fasta_lines, "^>")]
  
  parsed_df <- map_dfr(header_lines, function(header) {
    # Extract variant ID (e.g., MEFV_protein_path_variant_243)
    variant_name <- str_match(header, "^>([^|]+)")[, 2]
    
    # Extract sample string
    samples_part <- str_match(header, "samples=([^|]+)")[, 2]
    haps <- unlist(str_split(samples_part, ","))
    
    # NEW: Extract the ClinVarID trailing numeric value (e.g., ..._T.917511 -> 917511)
    # If it's a WT variant without a ClinVar dot ID, it assigns NA safely.
    clinvar_ids <- str_match(haps, "\\.(\\d+)$")[, 2]
    
    # Clean the haplotype names (e.g., 'NA18974.h1.chr16...' -> 'NA18974.h1')
    haps_clean <- str_match(haps, "^([^.]+\\.h[12])")[, 2]
    
    data.frame(
      Haplotype = haps_clean,
      Variant = variant_name,
      Origin = variant_type,
      ClinVarID = clinvar_ids,
      stringsAsFactors = FALSE
    )
  })
  
  return(parsed_df)
}

# Parse both FASTA files and combine
wt_variants <- parse_fasta(paste0("results/", gene, "/", gene, "_protein_haplotypes.fa"), "WT")
pg_variants <- parse_fasta(paste0("results/", gene, "/", gene, "_protein_path_haplotypes.fa"), "PG")
all_variants <- bind_rows(wt_variants, pg_variants)
dim(all_variants)
head(all_variants)


# ==============================================================================
# BEFORE WE CONTINUE - DO INTEGRITY CHECKS
# ==============================================================================

# Check A: Do we have haplotypes in the FASTA that DO NOT exist in the IGSR sample data?
unmapped_fasta <- all_variants$Haplotype[!(all_variants$Haplotype %in% igsr_samples_expanded$Haplotype)]
if (length(unmapped_fasta) > 0) {
  cat("\n[ERROR] FASTA haplotypes found that do not match IGSR metadata:\n")
  print(unique(unmapped_fasta))
  stop("Pipeline halted: FASTA variants contain orphan samples missing from your IGSR sheet.")
} else {
  cat("SUCCESS: All FASTA haplotypes map to IGSR metadata.\n")
}
# Check B: Do we have haplotypes in the expanded IGSR sample list that DO NOT exist in the FASTA haplotypes?
unmapped_igsr <- igsr_samples_expanded$Haplotype[!(igsr_samples_expanded$Haplotype %in% all_variants$Haplotype)]
if (length(unmapped_igsr) > 0) {
  cat("\n[ERROR] IGSR samples found that do not match FASTA haplotypes:\n")
  print(unique(unmapped_igsr))
  stop("Script halted")
} else {
  cat("SUCCESS: All IGSR samples map to the FASTA haplotypes.\n")
}

# now merge with the master haplotype sheet
final_df <- igsr_samples_expanded %>%
  left_join(all_variants, by = "Haplotype") %>%
  select(Sample_name, Haplotype, Variant, Origin, ClinVarID)

# reorder columns
final_df <- final_df[, c("Sample_name", "Haplotype", "Variant", "Origin", "ClinVarID")]

# save the output back in 'results' folder
write.table(final_df, paste0("results/", gene, "/", gene, "_protein_haplotypes_matched_to_igsr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            
            