
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
  header_lines <- fasta_lines[grep("^>", fasta_lines)]
  
  parsed_list <- list()
  for (header in header_lines) {
    # extract variant ID (e.g., MEFV_protein_path_variant_243)
    variant_name <- sub("^>([^|]+)\\|.*", "\\1", header)
    
    # extract sample string
    samples_part <- sub(".*samples=([^|]+).*", "\\1", header)
    haps <- unlist(strsplit(samples_part, ","))
    
    # clean the haplotype names (e.g., 'NA18974.h1.chr16...' -> 'NA18974.h1')
    haps_clean <- sub("^([^.]+\\.h[12]).*", "\\1", haps)
    
    parsed_list[[length(parsed_list) + 1]] <- data.frame(
      Haplotype = haps_clean,
      Variant = variant_name,
      Origin = variant_type,
      stringsAsFactors = FALSE
    )
  }
  return(do.call(rbind, parsed_list))
}

# parse both FASTA files for wild-type (WT) and pathogenic (PG) and combine
wt_variants <- parse_fasta(paste0("results/", gene, "/", gene, "_protein_haplotypes.fa"), "WT")
pg_variants   <- parse_fasta(paste0("results/", gene, "/", gene, "_protein_path_haplotypes.fa"), "PG")
all_variants <- rbind(wt_variants, pg_variants)
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
final_df <- merge(igsr_samples_expanded, all_variants, by = "Haplotype")

# reorder columns
final_df <- final_df[, c("Sample_name", "Haplotype", "Variant", "Origin")]

# save the output back in 'results' folder
write.table(final_df, paste0("results/", gene, "/", gene, "_protein_haplotypes_matched_to_igsr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
            
            