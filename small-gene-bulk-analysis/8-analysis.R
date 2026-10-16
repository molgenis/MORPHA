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

if(sum(grepl("_protein_variant", energy_min_per_haplo$haplotype)) != sum(grepl("_protein_variant", energy_all_per_haplo$haplotype))) { stop("Min energies ne all energies")}
cat("Total 1kGP haplotypes in energies:", sum(grepl("_protein_variant", energy_all_per_haplo$haplotype)), "\n")
cat("Total 1kGP hap x variant combinations energies:", sum(grepl("_protein_path_variant", energy_all_per_haplo$haplotype)), "\n")
cat("Total energy estimate per haplotype:", dim(energy_all_per_haplo)[2]-2, "\n")


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
      
      if (length(n_str) > 0) { n_vec[i] <- as.integer(sub("^n=", "", n_str[1])) }
      
      # Extract sample names
      samples_str <- parts[grep("^samples=", parts)]
      
      if (length(samples_str) > 0) {
        
        samples_val <- sub("^samples=", "", samples_str[1])
        
        # Extract variant IDs only for pathogenic haplotypes
        if (grepl("_protein_path_variant_", variant_id)) {
          
          # Split the comma-separated sample names
          sample_names <- strsplit(samples_val, ",", fixed = TRUE)[[1]]
          
          # Identify numeric variant IDs at the end of sample names
          has_variantID <- grepl("\\.[0-9]+$", sample_names)
          
          if (any(has_variantID)) {
            
            # Extract variant IDs from sample-name suffixes
            sample_variantIDs <- sub("^.*\\.([0-9]+)$", "\\1", sample_names[has_variantID])
            
            # Store the unique variant IDs for this haplotype
            variantID_vec[i] <- paste(unique(sample_variantIDs), collapse = ",")
            
            # Remove variant IDs from sample names
            #sample_names[has_variantID] <- sub("\\.[0-9]+$","",sample_names[has_variantID])
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
    haplo_list[[gene_name]] <- do.call(rbind, gene_haplotypes)
  }
}
all_haplotypes <- do.call(rbind, haplo_list)
rownames(all_haplotypes) <- haplo_list <- df_wide <- gene_haplotypes <- NULL
cat("Total nr of haplotypes:", dim(all_haplotypes)[1], "\n")


##########################################
# Remove genes without natural variation #
##########################################
genes_without_natural_variation <- all_haplotypes[grepl("_protein_variant_[0-9]+$", all_haplotypes[,"variant"]) & all_haplotypes[,"n"] == 6404, ]
genes_without_natural_variation[,1:4]
genes_to_drop <- genes_without_natural_variation$gene
cat("Dropping",length(genes_to_drop),"genes without natural variation\n")
all_haplotypes <- all_haplotypes[!all_haplotypes$gene %in% genes_to_drop, ]
cat("Total nr of haplotypes now:", dim(all_haplotypes)[1], "\n")
cat(" ..of which natural:", dim(all_haplotypes[grepl("_protein_variant_[0-9]+$", all_haplotypes[,"variant"]), ]), "\n")
cat(" ..of which pathogenic:", dim(all_haplotypes[grepl("_protein_path_variant_[0-9]+$", all_haplotypes[,"variant"]), ]), "\n")

energy_all_per_haplo <- energy_all_per_haplo[!energy_all_per_haplo$gene %in% genes_to_drop, ]
energy_min_per_haplo <- energy_min_per_haplo[!energy_min_per_haplo$gene %in% genes_to_drop, ]
if(sum(grepl("_protein_variant", energy_min_per_haplo$haplotype)) != sum(grepl("_protein_variant", energy_all_per_haplo$haplotype))) { stop("Min energies ne all energies")}
cat("Dropping from energies as well, now representing", dim(energy_all_per_haplo)[1], "proteins\n")
genes_without_natural_variation <- genes_to_drop <- NULL


############################################
# Connect natural to pathogenic haplotypes #
############################################
dt <- copy(as.data.table(all_haplotypes))
dt[, row_id := .I]
# expand the comma-separated sample+allele identifiers and add gene, variant
sample_long <- dt[, .(sample = trimws(unlist(strsplit(samples, ",", fixed = TRUE)))), by = row_id]
sample_long <- merge(sample_long, dt[, .(row_id, gene, variant)], by = "row_id", all.x = TRUE)
# split into natural and pathogenic protein haplotypes
natural_long <- unique(sample_long[!grepl("_protein_path_variant_", variant, fixed = TRUE), .(gene, sample, natural_haplotype = variant)])
pathogenic_long <- unique(sample_long[grepl("_protein_path_variant_", variant, fixed = TRUE), .(gene, sample, pathogenic_haplotype = variant)])
# split off variantID from sample in extra column
pathogenic_long[, variantID := sub(".*\\.", "", sample)]
pathogenic_long[, sample := sub("\\.[^.]+$", "", sample)]
# validate that each gene + sample+allele combination has exactly one natural background
natural_check <- natural_long[, .(n_natural = uniqueN(natural_haplotype)), by = .(gene, sample)][n_natural > 1]
if (nrow(natural_check) > 0) { stop("Some gene/sample combinations map to multiple natural haplotypes.") }
# merge natural and pathogenic long tables on gene and sample
merged_dt <- merge(pathogenic_long[, .(gene, sample, pathogenic_haplotype, variantID)],natural_long[, .(gene, sample, natural_haplotype)],by = c("gene", "sample"))
# group into an 'expanded' list of all natural haplotypes, introduced variant, and resulting pathogenic haplotype
haplo_link_expanded <- merged_dt[, .(n_samples = uniqueN(sample)), by = .(natural_haplo = natural_haplotype, pathogenic_haplo = pathogenic_haplotype, variantID)]
haplo_link_expanded[, gene := sub("_.*", "", natural_haplo)]
setcolorder(haplo_link_expanded, c("gene", "natural_haplo", "variantID", "pathogenic_haplo", "n_samples"))
#subset(haplo_link_expanded, gene=="CSTB" & variantID=="1791620") # this one was problematic before
dt <- sample_long <- natural_long <- pathogenic_long <- natural_check <- merged_dt <- NULL

# sanity check: collapse n_samples per unique variantID and check for variants where total_shared_alleles is NOT equal to 6404
variant_allele_collapse <- haplo_link_expanded[, .(
  total_shared_alleles = sum(n_samples, na.rm = TRUE),
  max_shared_alleles   = max(n_samples, na.rm = TRUE),
  n_occurrences        = .N
), by = variantID]
anomalous_variants <- variant_allele_collapse[total_shared_alleles != 6404]
print(anomalous_variants)
anomalous_variants <- NULL

# collapse into unique natural-patho haplo pairs
haplo_link_collapse_unique_pairs <- haplo_link_expanded[, .(
  variantID = paste(sort(unique(variantID)), collapse = ","),
  n_samples = paste(sort(unique(n_samples)), collapse = ","),
  total_variants_count = uniqueN(variantID)
), by = .(gene, natural_haplo, pathogenic_haplo)]

cat("Nr of haplotype relationships:",dim(haplo_link_collapse_unique_pairs)[1],"\n")
cat("Nr of unique pathogenic haplotypes represented:",length(unique(haplo_link_collapse_unique_pairs$pathogenic_haplo)),"\n")
cat("Nr of variants affecting alleles:",dim(variant_allele_collapse)[1],"\n")

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
setnames(natural_energy, old = c("haplotype", reu_cols), new = c("natural_haplo", paste0(reu_cols, "_natural")))
pathogenic_energy <- copy(energy_all_per_haplo[, c("gene", "haplotype", reu_cols), with = FALSE])
setnames(pathogenic_energy, old = c("haplotype", reu_cols), new = c("pathogenic_haplo", paste0(reu_cols, "_pathogenic")))

# add natural and pathogenic haplotype energies
haplo_link_expanded <- merge(haplo_link_expanded, natural_energy, by = c("gene", "natural_haplo"), all.x = TRUE, sort = FALSE)
haplo_link_expanded <- merge(haplo_link_expanded, pathogenic_energy, by = c("gene", "pathogenic_haplo"), all.x = TRUE, sort = FALSE)

# calculate DDG for each corresponding REU pair, i.e. pathogenic energy - natural energy
for (i in 1:100) {
  set( haplo_link_expanded, j = paste0("DDG", i), value = haplo_link_expanded[[paste0("REU", i, "_pathogenic")]] - haplo_link_expanded[[paste0("REU", i, "_natural")]] )
}

# cleanup
natural_energy <- pathogenic_energy <- reu_cols <- NULL

###############################
# Test for significance IN PROGRESS...
##########################

# remove DG columns, no longer needed
dim(haplo_link_expanded)
reu_natural_cols <- paste0("REU", 1:100, "_natural")
haplo_link_expanded[, (reu_natural_cols) := NULL]
reu_path_cols <- paste0("REU", 1:100, "_pathogenic")
haplo_link_expanded[, (reu_path_cols) := NULL]
dim(haplo_link_expanded)

# grab all energy columns (DDG1 to DDG100)
ddg_cols <- grep("^DDG[0-9]+$", names(haplo_link_expanded), value = TRUE)


# 3. Melt the wide DDG columns into long format, keeping variantID, gene, and the groups (haplo_pair)
# We create 'haplo_pair' to group distinct natural-pathogenic trajectories per variant
haplo_link_expanded[, haplo_pair := paste(natural_haplo, pathogenic_haplo, sep = "__vs__")]

long_dt <- melt(
  haplo_link_expanded,
  id.vars = c("gene", "variantID", "haplo_pair"),
  measure.vars = ddg_cols,
  variable.name = "trajectory",
  value.name = "DDG"
)

# Filter out missing values to prevent test errors
long_dt <- long_dt[!is.na(variantID) & !is.na(DDG)]

# 4. Run Kruskal-Wallis test and calculate Epsilon-squared effect size per variantID
variant_results <- long_dt[, {
  k <- uniqueN(haplo_pair)
  n <- .N
  
  if (n > 2 && k > 1) {
    kt <- kruskal.test(DDG ~ haplo_pair)
    
    h_stat <- kt$statistic
    p_val  <- kt$p.value
    
    # Epsilon-squared effect size: (H - k + 1) / (n - k)
    eps_sq <- (h_stat - k + 1) / (n - k)
    eps_sq <- max(0, min(1, eps_sq, na.rm = TRUE), na.rm = TRUE)
    
    list(
      pvalue      = as.numeric(p_val),
      effectsize  = as.numeric(eps_sq),
      n_groups    = as.integer(k)
    )
  } else {
    list(
      pvalue      = NA_real_,
      effectsize  = NA_real_,
      n_groups    = as.integer(k)
    )
  }
}, by = .(gene, variantID)]

# apply Bonferroni nultiple testing correction and 1% significance
variant_results[!is.na(pvalue), adj_p_value := p.adjust(pvalue, method = "bonferroni")]
variant_results[, significant_bonf_1pct := fifelse(is.na(adj_p_value), NA_character_, fifelse(adj_p_value < 0.01, "YES", "NO"))]
pct_summary <- variant_results[, .N, by = significant_bonf_1pct][, percentage := (N / sum(N)) * 100]
print(pct_summary)

# visualization of the results (Volcano-style plot, highlighting Bonferroni hits)
variant_results[, neg_log_p := -log10(pvalue)]

ggplot(variant_results[is.finite(neg_log_p)], aes(x = effectsize, y = neg_log_p, color = significant_bonf_1pct)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("NO" = "steelblue", "YES" = "firebrick")) +
  labs(
    title = "Kruskal-Wallis Results with Bonferroni Correction",
    subtitle = "Significance threshold: adjusted p < 0.01 (Red)",
    x = "Effect Size (Epsilon-Squared)",
    y = "-log10(P-value)",
    color = "Significant (1% Bonf)"
  ) +
  theme_minimal()




# Pick a specific variantID from your results to inspect
target_variant <- variant_results[order(pvalue)][1, variantID] # selects the most significant one
target_variant <- "987020"

ggplot(long_dt[variantID == target_variant], aes(x = haplo_pair, y = DDG, fill = haplo_pair)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
  labs(
    title = paste("Trajectory Stability Distribution for Variant:", target_variant),
    x = "Pathogenic Haplotype Pathway",
    y = "DDG (Pathogenic - Natural)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")


