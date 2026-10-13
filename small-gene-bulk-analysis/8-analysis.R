library(R.utils)
library(data.table)

setwd("~/git/MORPHA/")

#############################################
# Read in the target genes of this analysis #
#############################################
genes <- read.csv("small-gene-bulk-analysis/uniprot_signal_peptides.csv", stringsAsFactors = FALSE)
cat("Total gene targets:", dim(genes)[1], "\n")
head(genes)

########################
# Read in variant VCFs #
########################
var_list <- list()
for (i in seq_len(nrow(genes))) {
  gene <- genes$gene[i]
  status <- genes$status[i]
  vcf_path <- file.path("small-gene-bulk-analysis", "genes", gene, "variants.vcf.gz")
  if (file.exists(vcf_path)) {
    #cat("Reading VCF for:", gene, "\n")
    vcf <- fread(vcf_path, skip = "#CHROM")
    setnames(vcf, "#CHROM", "CHROM")
    cols_to_drop <- intersect(c("QUAL", "FILTER", "INFO"), names(vcf))
    if (length(cols_to_drop) > 0) {
      vcf[, (cols_to_drop) := NULL]
    }
    vcf$gene <- gene
    var_list[[gene]] <- vcf
  } else {
    cat("Warning: File does not exist ->", vcf_path, "\n")
  }
}
vcf<-NULL
genes_and_variants <- do.call(rbind, var_list)
var_list <- NULL
cat("Total genes with variants:", length(unique(genes_and_variants$gene)), "\n")
cat("Total number of variants:", dim(genes_and_variants)[1], "\n")
head(genes_and_variants)

#################################################
# Read in haplotype energy minimization results #
#################################################
em_list <- list()
em_all_list <- list()
for (gene_name in unique(genes_and_variants$gene)) {
  #cat("Reading haplo energy minimization for gene", gene_name,"\n")
  em_path <- file.path("results", gene_name, paste0(gene_name,"_minimum_per_haplotype.csv"))
  em_all_path <- file.path("results", gene_name, paste0(gene_name,"_all_dg_per_haplotype.csv"))
  # redundant since minimums are also in "all_dg", but minimums ARE the key result
  if (file.exists(em_path)) {
    df <- read.csv(em_path, header = FALSE, col.names = c("haplotype", "REU"), stringsAsFactors = FALSE)
    df$gene <- gene_name
    em_list[[gene_name]] <- df
  } else {
    cat("Warning: minimum_per_haplotype missing for gene", gene_name, "->", em_path, "\n")
  }
  if (file.exists(em_all_path)) {
    df <- read.csv(em_all_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c("haplotype", paste0("REU", 1:(ncol(df) - 1)))
    df$gene <- gene_name
    em_all_list[[gene_name]] <- df
  } else {
    cat("Warning: all_dg_per_haplotype missing for gene", gene_name, "->", em_path, "\n")
  }
}
energy_min_per_haplo <- do.call(rbind, em_list)
energy_all_per_haplo <- do.call(rbind, em_all_list)
df<-NULL
em_list<-NULL
em_all_list<-NULL
rownames(energy_min_per_haplo) <- NULL
rownames(energy_all_per_haplo) <- NULL
cat("Total haplotype x variant combinations in minimums:", dim(energy_min_per_haplo)[1], "\n")
cat("Total haplotype x variant combinations in all DGs:", dim(energy_all_per_haplo)[1], "\n")
head(energy_min_per_haplo)
energy_all_per_haplo[1:6,1:6]

######################################################
# Read in the protein haplotypes for 1kG individuals #
######################################################
haplo_list <- list()

for (gene_name in unique(genes_and_variants$gene)) {
  
  # Paths to both FASTA files
  fasta_paths <- c(file.path("results", gene_name, paste0(gene_name, "_protein_haplotypes.fa")), file.path("results", gene_name,paste0(gene_name, "_protein_path_haplotypes.fa")))
  gene_haplotypes <- list()
  
  for (haplo_path in fasta_paths) {
    
    # Skip missing FASTA files
    if (!file.exists(haplo_path)) {
      cat("Warning: FASTA file does not exist ->", haplo_path, "\n")
      next
    }
    
    fasta_lines <- readLines(haplo_path, warn = FALSE)
    header_idx <- grep("^>", fasta_lines)
    
    if (length(header_idx) == 0) {
      cat("Warning: No FASTA headers found in", haplo_path, "\n")
      next
    }
    
    # Initialize output vectors
    n_records <- length(header_idx)
    
    gene_vec      <- character(n_records)
    variant_vec   <- character(n_records)
    n_vec         <- rep(NA_integer_, n_records)
    samples_vec   <- rep(NA_character_, n_records)
    sequence_vec  <- character(n_records)
    variantID_vec <- rep(NA_character_, n_records)
    
    # Parse each FASTA record
    for (i in seq_along(header_idx)) {
      
      h_line <- sub("^>", "", fasta_lines[header_idx[i]])
      parts <- strsplit(h_line, "\\|")[[1]]
      
      # Extract haplotype identifier
      variant_id <- parts[1]
      
      gene_vec[i] <- gene_name
      variant_vec[i] <- variant_id
      
      # Extract sample count
      n_str <- parts[grep("^n=", parts)]
      
      if (length(n_str) > 0) {
        n_vec[i] <- as.integer(sub("^n=", "", n_str[1]))
      }
      
      # Extract sample names
      samples_str <- parts[grep("^samples=", parts)]
      
      if (length(samples_str) > 0) {
        
        samples_val <- sub("^samples=", "", samples_str[1])
        
        # Extract variant IDs only for pathogenic haplotypes
        if (grepl("_protein_path_variant_", variant_id)) {
          
          # Split the comma-separated sample names
          sample_names <- strsplit(
            samples_val,
            ",",
            fixed = TRUE
          )[[1]]
          
          # Identify numeric variant IDs at the end of sample names
          has_variantID <- grepl("\\.[0-9]+$", sample_names)
          
          if (any(has_variantID)) {
            
            # Extract variant IDs from sample-name suffixes
            sample_variantIDs <- sub(
              "^.*\\.([0-9]+)$",
              "\\1",
              sample_names[has_variantID]
            )
            
            # Store the unique variant IDs for this haplotype
            variantID_vec[i] <- paste(
              unique(sample_variantIDs),
              collapse = ","
            )
            
            # Remove variant IDs from sample names
            sample_names[has_variantID] <- sub(
              "\\.[0-9]+$",
              "",
              sample_names[has_variantID]
            )
          }
          
          # Recombine the cleaned sample names
          samples_vec[i] <- paste(sample_names, collapse = ",")
        } else {
          # Ordinary protein haplotypes: no variant ID
          samples_vec[i] <- samples_val
        }
      }
      
      # Determine sequence boundaries
      start_seq <- header_idx[i] + 1
      
      end_seq <- if (i < length(header_idx)) {
        header_idx[i + 1] - 1
      } else {
        length(fasta_lines)
      }
      
      # Concatenate sequence lines
      if (start_seq <= end_seq) {
        sequence_vec[i] <- paste(fasta_lines[start_seq:end_seq], collapse = "")
      } else {
        sequence_vec[i] <- ""
      }
    }
    
    # Build data frame for this FASTA file
    df_wide <- data.frame(gene = gene_vec, variant = variant_vec, n = n_vec, variantID = as.character(variantID_vec), samples = samples_vec, sequence = sequence_vec, stringsAsFactors = FALSE)
    gene_haplotypes[[length(gene_haplotypes) + 1]] <- df_wide
  }
  
  # Combine both FASTA files for this gene
  if (length(gene_haplotypes) > 0) {
    haplo_list[[gene_name]] <- do.call(
      rbind,
      gene_haplotypes
    )
  }
}

haplo_per_indiv <- do.call(rbind, haplo_list)
rownames(haplo_per_indiv) <- NULL
haplo_list <- NULL
df_wide <- NULL
cat("Total nr of haplotypes:", dim(haplo_per_indiv)[1], "\n")
haplo_per_indiv[1:5,]
haplo_per_indiv[1:15,1:4]

