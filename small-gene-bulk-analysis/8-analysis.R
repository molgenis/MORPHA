library(R.utils)
library(data.table)

setwd("~/git/MORPHA/")

###################################
# Read in genes for this analysis #
###################################
genes <- read.csv("small-gene-bulk-analysis/uniprot_signal_peptides.csv", stringsAsFactors = FALSE)
cat("Total genes:", dim(genes)[1], "\n")

########################
# Read in variant VCFs #
########################
genes_and_variants <- list()
for (i in seq_len(nrow(genes))) {
  gene <- genes$gene[i]
  status <- genes$status[i]
  vcf_path <- file.path("small-gene-bulk-analysis", "genes", gene, "variants.vcf.gz")
  if (file.exists(vcf_path)) {
    cat("Reading VCF for:", gene, "\n")
    vcf <- fread(vcf_path, skip = "#CHROM")
    setnames(vcf, "#CHROM", "CHROM")
    cols_to_drop <- intersect(c("QUAL", "FILTER", "INFO"), names(vcf))
    if (length(cols_to_drop) > 0) {
      vcf[, (cols_to_drop) := NULL]
    }
    genes_and_variants[[gene]] <- vcf
  } else {
    cat("Warning: File does not exist ->", vcf_path, "\n")
  }
}
vcf<-NULL
cat("Total genes with variants:", length(genes_and_variants), "\n")
cat("Total number of variants:", sum(sapply(genes_and_variants, nrow)), "\n")

#################################################
# Read in haplotype energy minimization results #
#################################################
em_list <- list()
for (gene_name in names(genes_and_variants)) {
  cat("Reading haplo energy minimization for gene", gene_name,"\n")
  em_path <- file.path("results", gene_name, paste0(gene_name,"_minimum_per_haplotype.csv"))
  if (file.exists(em_path)) {
    df <- read.csv(em_path, header = FALSE, col.names = c("haplotype", "REU"), stringsAsFactors = FALSE)
    df$gene <- gene_name
    em_list[[gene_name]] <- df
  } else {
    cat("Warning: File missing for gene", gene_name, "->", em_path, "\n")
  }
  # TODO also load in full results (all totals, for stat. testing)
}
energy_min_per_haplo <- do.call(rbind, em_list)
df<-NULL
em_list<-NULL
rownames(energy_min_per_haplo) <- NULL
cat("Total haplotype x variant combinations:", dim(energy_min_per_haplo)[1], "\n")

################################################
# Read in haplo mapping AA seq and mapping to individuals/alleles
################################################
df_list <- list()
for (gene_name in names(genes_and_variants)) {
  cat("Reading haplo energy minimization for gene", gene_name,"\n")
  
}



# 1. Read FASTA file lines
gene_name<-"UFM1"
haplo_path <- file.path("results", gene_name, paste0(gene_name,"_protein_haplotypes.fa"))
fasta_lines <- readLines(haplo_path)

# 2. Identify header positions
header_idx <- grep("^>", fasta_lines)

# 3. Initialize vectors for columns
gene_vec     <- character(length(header_idx))
variant_vec  <- character(length(header_idx))
n_vec        <- integer(length(header_idx))
samples_vec  <- character(length(header_idx))
seq_vec      <- character(length(header_idx))

# 4. Parse headers and extract multiline sequences
for (i in seq_along(header_idx)) {
  h_line <- fasta_lines[header_idx[i]]
  
  # Remove leading '>'
  h_clean <- sub("^>", "", h_line)
  
  # Split header components by '|'
  parts <- strsplit(h_clean, "\\|")[[1]]
  
  # Extract variant name and gene name
  variant_id <- parts[1]
  gene_name  <- strsplit(variant_id, "_")[[1]][1]
  
  # Extract 'n' count
  n_str <- parts[grep("^n=", parts)]
  n_val <- as.integer(sub("n=", "", n_str))
  
  # Extract samples string
  samples_str <- parts[grep("^samples=", parts)]
  samples_val <- sub("samples=", "", samples_str)
  
  # Determine start and end lines for the sequence
  start_seq <- header_idx[i] + 1
  end_seq   <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(fasta_lines)
  
  # Concatenate sequence lines
  seq_str <- paste(fasta_lines[start_seq:end_seq], collapse = "")
  
  # Store in vectors
  gene_vec[i]    <- gene_name
  variant_vec[i] <- variant_id
  n_vec[i]       <- n_val
  samples_vec[i] <- samples_val
  seq_vec[i]     <- seq_str
}

# 5. Build final wide-format data frame
df_wide <- data.frame(
  gene = gene_vec,
  variant = variant_vec,
  n = n_vec,
  samples = samples_vec,
  sequence = seq_vec,
  stringsAsFactors = FALSE
)

# Inspect the result
head(df_wide)





# Combine

