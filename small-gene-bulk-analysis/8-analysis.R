library(R.utils)
library(data.table)
library(this.path)
library(ggplot2)

script_dir <- this.dir()
parent_dir <- dirname(script_dir)
setwd(dirname(script_dir))

#############################################
# Read in the target genes of this analysis #
#############################################
genes <- read.csv("small-gene-bulk-analysis/uniprot_signal_peptides.csv", stringsAsFactors = FALSE)
cat("Total gene targets:", dim(genes)[1], "\n")
head(genes)

############################################
# Read in their corresponding variant VCFs #
############################################
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
cat("Total 1kGP haplotypes in minimums:", sum(grepl("_protein_variant", energy_min_per_haplo$haplotype)), "\n")
cat("Total 1kGP haplotypes in all DGs:", sum(grepl("_protein_variant", energy_all_per_haplo$haplotype)), "\n")
cat("Total 1kGP hap x variant combinations in minimums:", sum(grepl("_protein_path_variant", energy_min_per_haplo$haplotype)), "\n")
cat("Total 1kGP hap x variant combinations in all DGs:", sum(grepl("_protein_path_variant", energy_all_per_haplo$haplotype)), "\n")
cat("Total nr of delta Gibbs REU estimate across all haplotype:", dim(energy_all_per_haplo)[2]-2, "\n")
head(energy_min_per_haplo)
energy_all_per_haplo[1:10,1:6]

#######################################################
# Read in the protein haplotypes for 1kGP individuals #
#######################################################
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


# since e have already separated 'variantID' from the sample/haplotype names, so remove duplicate samples
haplo_per_indiv$samples <- sapply(haplo_per_indiv$samples, function(x) {
  if (is.na(x) || x == "") return(x)
  items <- unlist(strsplit(x, ",", fixed = TRUE))
  unique_items <- unique(items)
  paste(unique_items, collapse = ",")
})

cat("Total nr of haplotypes:", dim(haplo_per_indiv)[1], "\n")
genes_without_natural_variation <- haplo_per_indiv[grepl("_protein_variant_[0-9]+$", haplo_per_indiv[,"variant"]) & haplo_per_indiv[,"n"] == 6404, ]
genes_without_natural_variation[,1:4]
genes_to_drop <- genes_without_natural_variation$gene
cat("Dropping",length(genes_to_drop),"genes without natural variation\n")
haplo_per_indiv <- haplo_per_indiv[!haplo_per_indiv$gene %in% genes_to_drop, ]
cat("Total nr of haplotypes now:", dim(haplo_per_indiv)[1], "\n")
cat(" ..of which natural:", dim(haplo_per_indiv[grepl("_protein_variant_[0-9]+$", haplo_per_indiv[,"variant"]), ]), "\n")
cat(" ..of which pathogenic:", dim(haplo_per_indiv[grepl("_protein_path_variant_[0-9]+$", haplo_per_indiv[,"variant"]), ]), "\n")
#haplo_per_indiv[1:5,]
#haplo_per_indiv[1:15,1:4]

# cleanup
rownames(haplo_per_indiv) <- haplo_list <- df_wide <- genes_without_natural_variation <- genes_to_drop <- NULL

#multiple_varIDs <- subset(haplo_per_indiv, grepl(",", variantID, fixed = TRUE))
#dim(multiple_varIDs)
#multiple_varIDs[1:5,1:4]

#all_natural_variation <- haplo_per_indiv[grepl("_protein_variant_[0-9]+$", haplo_per_indiv[,"variant"]), ]
#all_induced_path_variation <- haplo_per_indiv[grepl("_protein_path_variant_[0-9]+$", haplo_per_indiv[,"variant"]), ]

#dim(all_natural_variation)
#colnames(all_natural_variation)
#all_natural_variation[1:5,1:4]

#dim(all_induced_path_variation)
#all_induced_path_variation[1:5,1:4]

############################################
# Connect natural to pathogenic haplotypes #
############################################

dt <- copy(as.data.table(haplo_per_indiv))
dt[, row_id := .I]

# expand the comma-separated sample+allele identifiers and add gene, variant and variantID
sample_long <- dt[, .(sample = trimws(unlist(strsplit(samples, ",", fixed = TRUE)))), by = row_id]
sample_long <- merge(sample_long, dt[, .(row_id, gene, variant, variantID)], by = "row_id", all.x = TRUE)

# split into natural and pathogenic protein haplotypes
natural_long <- unique(sample_long[!grepl("_protein_path_variant_", variant, fixed = TRUE), .(gene, sample, natural_haplotype = variant)])
pathogenic_long <- unique(sample_long[grepl("_protein_path_variant_", variant, fixed = TRUE), .(gene, sample, pathogenic_haplotype = variant, variantID)])

# validate that each gene + sample+allele combination has exactly one natural background
natural_check <- natural_long[, .(n_natural = uniqueN(natural_haplotype)), by = .(gene, sample)][n_natural > 1]
if (nrow(natural_check) > 0) { stop("Some gene/sample combinations map to multiple natural haplotypes.") }

# grab and attach the natural background for the pathogenic haplotypes, since we know they are unique by gene x sample combination
matched <- merge(pathogenic_long, natural_long, by = c("gene", "sample"), allow.cartesian = TRUE)

# collapse into one row per unique combination of gene, natural haplotype, pathogenic haplotype, and variantID
haplotype_relationships <- matched[, .(n_shared_sample_alleles = uniqueN(sample)), by = .(gene, natural_haplotype, pathogenic_haplotype, variantID)]

# expanded mapping of variantID to corresponding pathogenic allele
pathogenic_variant_map <- unique(dt[grepl("_protein_path_variant_", variant, fixed = TRUE), .(gene, pathogenic_haplotype = variant, variantID)])[, .(variantID = trimws(unlist(strsplit(variantID, ",", fixed = TRUE)))), by = .(gene, pathogenic_haplotype)]

# cleanup
dt <- sample_long <- natural_long <- pathogenic_long <- natural_check <- matched <- NULL

cat("Nr of haplotype relationships:",dim(haplotype_relationships)[1],"\n")
cat("Nr of unique pathogenic haplotypes represented:",length(unique(haplotype_relationships$pathogenic_haplotype)),"\n")
cat("Nr of variants affecting alleles:",dim(pathogenic_variant_map)[1],"\n")

#####################################################
# Add DG and DDG for pathogenic and base haplotypes #
#####################################################

# convert to data table
setDT(energy_all_per_haplo)

# the 100 energy columns
reu_cols <- paste0("REU", 1:100)

# check that each gene + haplotype has exactly one energy row
if (energy_all_per_haplo[, .N, by = .(gene, haplotype)][N > 1, .N] > 0) { stop("energy_all_per_haplo has duplicate gene + haplotype keys.") }

# prepare natural and pathogenic haplotype energies
natural_energy <- copy(energy_all_per_haplo[, c("gene", "haplotype", reu_cols), with = FALSE])
setnames(natural_energy, old = c("haplotype", reu_cols), new = c("natural_haplotype", paste0(reu_cols, "_natural")))
pathogenic_energy <- copy(energy_all_per_haplo[, c("gene", "haplotype", reu_cols), with = FALSE])
setnames(pathogenic_energy, old = c("haplotype", reu_cols), new = c("pathogenic_haplotype", paste0(reu_cols, "_pathogenic")))

# add natural and pathogenic haplotype energies
haplotype_relationships <- merge(haplotype_relationships, natural_energy, by = c("gene", "natural_haplotype"), all.x = TRUE, sort = FALSE)
haplotype_relationships <- merge(haplotype_relationships, pathogenic_energy, by = c("gene", "pathogenic_haplotype"), all.x = TRUE, sort = FALSE)

# calculate DDG for each corresponding REU pair, i.e. pathogenic energy - natural energy
for (i in 1:100) {
  set( haplotype_relationships, j = paste0("DDG", i), value = haplotype_relationships[[paste0("REU", i, "_pathogenic")]] - haplotype_relationships[[paste0("REU", i, "_natural")]] )
}

# cleanup
natural_energy <- pathogenic_energy <- reu_cols <- NULL

########################
# Test for significance
###################



