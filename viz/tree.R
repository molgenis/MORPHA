#!/usr/bin/env Rscript

library(seqinr)
library(tidyverse)
library(phangorn)
library(ggtree)
library(ape)

gene <- "MEFV"
setwd("~/git/MORPHA/")

fasta_file <- paste0("results/", gene, "/", gene, "_protein_haplotypes.fa")
metadata_file <- paste0("viz/igsr_samples.tsv")
exclude_samples <- "HG01783"

count_string <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  paste(paste0(names(tab), "=", as.integer(tab)), collapse = ", ")
}

dominant <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return("unknown")
  names(sort(table(x), decreasing = TRUE))[1]
}

# Read FASTA: headers look like variant|n=...|samples=sample.h1,sample.h2,...
fasta <- read.fasta(fasta_file, seqtype = "AA", as.string = TRUE,
                    whole.header = TRUE, strip.desc = FALSE)
headers <- names(fasta)
seqs <- vapply(fasta, paste0, character(1))
variant <- sub("\\|.*$", "", headers)
names(seqs) <- variant

header_tbl <- tibble(
  label = variant,
  samples = str_match(headers, "samples=([^|]+)$")[, 2],
  sequence_length = nchar(seqs)
)

# Summarize subpopulation and superpopulation composition per protein variant
metadata <- read.delim(metadata_file, sep = "\t", check.names = FALSE)

## CUSTOM LABELS
pop_count_string <- function(x, max_pops = Inf) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return("no mapped samples")
  tab <- sort(table(x), decreasing = TRUE)
  # Optional: collapse smaller populations into "rest"
  if (is.finite(max_pops) && length(tab) > max_pops) {
    shown <- tab[seq_len(max_pops)]
    rest <- sum(tab[-seq_len(max_pops)])
    tab <- c(shown, rest = rest)
  }
  paste0(names(tab), "=", as.integer(tab), collapse = ", ")
}

ann <- header_tbl %>%
  select(label, samples) %>%
  separate_rows(samples, sep = ",") %>%
  rename(sample_haplotype = samples) %>%
  mutate(sample_id = sub("\\.h[0-9]+$", "", sample_haplotype)) %>%
  filter(!sample_id %in% exclude_samples) %>%
  left_join(metadata, by = c("sample_id" = "Sample name")) %>%
  group_by(label) %>%
  summarise(
    n_haplotypes_after_exclusion = n(),
    subpopulation_counts = pop_count_string(`Population code`, max_pops = 3),
    superpopulation = dominant(`Superpopulation code`),
    superpopulation_label = pop_count_string(`Superpopulation code`, max_pops = Inf),
    .groups = "drop"
  ) %>%
  right_join(header_tbl %>% select(label, sequence_length), by = "label") %>%
  mutate(
    n_haplotypes_after_exclusion = replace_na(n_haplotypes_after_exclusion, 0L),
    subpopulation_counts = replace_na(subpopulation_counts, "no mapped samples"),
    superpopulation = replace_na(superpopulation, "unknown"),
    superpopulation_label = replace_na(superpopulation_label, "unknown"),
    
    tip_label = paste0(
      label,
      " (n=", n_haplotypes_after_exclusion,
      ", ", subpopulation_counts,
      ")"
    )
  )

write.table(ann, "MEFV_protein_tree_tip_annotations.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Align with MAFFT if lengths differ; phangorn expects aligned AA characters
if (length(unique(nchar(seqs))) > 1) {
  mafft <- Sys.which("mafft")
  if (mafft == "") {
    stop("Sequences have unequal lengths; install MAFFT, e.g. conda install -c bioconda mafft")
  }
  in_fa <- tempfile(fileext = ".fa")
  out_fa <- tempfile(fileext = ".fa")
  write.fasta(as.list(seqs), names(seqs), file.out = in_fa, as.string = TRUE)
  system2(mafft, c("--auto", "--quiet", in_fa), stdout = out_fa)
  aln <- read.fasta(out_fa, seqtype = "AA", as.string = TRUE)
  seqs <- vapply(aln, paste0, character(1))
}

stopifnot(length(unique(nchar(seqs))) == 1)
aa_matrix <- do.call(rbind, strsplit(seqs, ""))
rownames(aa_matrix) <- names(seqs)

# Protein distance tree: JTT AA model + neighbor joining
# CAN GIVE NEG DISTANCE:
aa <- phyDat(aa_matrix, type = "AA")
tree <- NJ(dist.ml(aa, model = "JTT", exclude = "pairwise")) %>% ladderize()
tree$edge.length[tree$edge.length < 0] <- 0 # set neg to 0
# ALTERNATIVE HC:
#aa <- phyDat(aa_matrix, type = "AA")
#d <- dist.ml(aa, model = "JTT", exclude = "pairwise")
#tree <- upgma(d) %>% ladderize()
#write.tree(tree, "MEFV_protein_tree.nwk")

superpop_cols <- c(AFR = "#E41A1C", AMR = "#377EB8", EAS = "#4DAF4A",
                   EUR = "#984EA3", SAS = "#FF7F00", unknown = "grey50")




p <- ggtree(tree) %<+% ann +
  geom_tippoint(aes(color = superpopulation), size = 2.2) +
  geom_tiplab(aes(label = tip_label, color = superpopulation),
              size = 4, align = TRUE, linetype = "dotted", linesize = 0.2, show.legend = FALSE) +
  scale_color_manual(values = superpop_cols, na.value = "grey50") +
  xlim(0,0.007)+
  theme_tree2() +
  theme(
    plot.title = element_text(size = 24, color = "black"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.key.size = unit(1.2, "cm")
    ) +
  guides(
    colour = guide_legend(override.aes = list(size = 5)),
    fill   = guide_legend(override.aes = list(size = 5))
  ) +
  labs(title = "MEFV protein haplotype tree across global ancestries",
       subtitle = "Based on maximum-likelihood amino-acid evolutionary distance under the JTT substitution model on all 57 unique MEFV protein haplotypes present in 6,402 haplotypes from 3,201 individuals in the 1000 Genomes Project",
       x = "Branch length (JTT amino-acid substitutions per site)",
       color = "Dominant\nsuperpopulation")
p
#ggsave("MEFV_protein_tree_dendrogram.pdf", p, width = 16, height = 12, limitsize = FALSE)
ggsave("MEFV_protein_tree_dendrogram.png", p, width = 16, height = 12, dpi = 300, limitsize = FALSE)

