# requires:
# results/<gene>/<gene>_protein_haplotypes_quant.tsv -> produced by viz/getQuant-FoldX.R
# viz/igsr_samples.tsv -> static resource
# results/<gene>/variant_info.tsv -> manually created

library(dplyr)
library(ggplot2)
library(maps)
library(tidyr)
library(stringr)
library(ggrepel)
library(rstatix)
library(tidytext)

genes <- c("TTR", "GBA1", "MLH1", "TP53", "MEFV")
genes <- c("TTR")

setwd("~/git/MORPHA/")

# read IGSR sample metadata and clean up
samples <- read.delim("viz/igsr_samples.tsv", sep = "\t", check.names = FALSE)
samples$`Biosample ID` <- NULL
samples$`Data collections` <- NULL
samples$`Superpopulation code` <- NULL
samples$`Superpopulation name` <- NULL
samples$`Population elastic ID` <- NULL
samples <- samples %>%
  rename(population_code = `Population code`) %>%
  rename(population_name = `Population name`) %>%
  filter(`Sample name` != "HG01783") # drop 1 mixed ancestry individual

head(samples)
dim(samples) # 3201 instead of 3202
unique(samples$population_code)

# Initialize empty lists to store results per gene
samples_quant_all <- data.frame()

# Loop through each gene individually
for (gene in genes) {
  
  # debug
  #gene = "TTR"
  
  # read quantifications
  quant <- read.delim( paste0("results/", gene, "/", gene, "_protein_haplotypes_quant.tsv"),  sep = "\t", check.names = FALSE)
  quant <- quant %>% filter(Sample_name != "HG01783")
  head(quant)
  
  # Check: Do we have sample names in IGSR list that DO NOT exist in the quantifications?
  unmapped_samplesA <- samples$`Sample name`[!(samples$`Sample name` %in% quant$Sample_name)]
  unmapped_samplesA
  # Check: Do we have sample names in quantification that DO NOT exist in the IGSR list?
  unmapped_samplesB <- quant$Sample_name[!(quant$Sample_name %in% samples$`Sample name`)]
  unmapped_samplesB
  
  # Merge the data, now in long format
  samples_quant_long <- quant %>%
    left_join(samples, by = c("Sample_name" = "Sample name"))
  
  # Add variant info
  variant_info <- read.delim( paste0("results/", gene, "/", "variant_info.tsv"),  sep = "\t", check.names = FALSE)
  samples_quant_long <- samples_quant_long %>%
    left_join(variant_info, by = c("VariantID" = "id"))
  
  ##########################################################################################
  # Pivot the data so that WT/pathogenic get their own columns and each row is a haplotype #
  ##########################################################################################
  
  # Isolate the haplotype-specific WT energies as baseline
  wt_baselines <- samples_quant_long %>%
    filter(Origin == "WT") %>%
    select(Haplotype, WT_Energy = Total_Energy)
  
  # Count occurence of energy levels and get top one as reference value (assuming this is the WT/reference genome)
  wt_energy_counts <- wt_baselines %>% dplyr::count(WT_Energy, name = "Count", sort = TRUE)
  head(wt_energy_counts)
  
  # !!! TP53 hg38 is not the most common allele - here, select index 2 if gene is "TP53", otherwise select index 1
  # for TTR, GBA1, and MLH1, the top one is also hg38.
  top_ref_energy <- ifelse(gene == "TP53", wt_energy_counts[[1]][2], wt_energy_counts[[1]][1])
  ref_energy <- top_ref_energy
  
  # Process and reshape the data frame based strictly on VariantID
  samples_quant <- samples_quant_long %>%
    # Filter down to just the pathogenic rows (.. including Benign ..)
    filter(Origin == "PG") %>%
    
    # Attach the correct WT background energy for this specific haplotype
    left_join(wt_baselines, by = "Haplotype") %>%
    
    # Calculate DDG using the individual's specific WT background
    mutate(DDG = Total_Energy - WT_Energy)
  
  # View the flattened haplotype-level data frame
  head(samples_quant)
  
  samples_quant <- samples_quant %>% select(Haplotype, population_code, label, DDG) %>% mutate(gene = gene)
  
  # combine all per-gene data frames into one master data frame
  samples_quant_all <- bind_rows(samples_quant_all, samples_quant)
}

# View summary of combined data frame
head(samples_quant_all)

# prepare data: combine gene and variant label for the X-axis categories
df_plot <- samples_quant_all %>%
  filter(!is.na(DDG)) %>%
  mutate(gene_variant = paste0(gene, "-", label))




##################
#### PLOTTING ####
##################


## WITH GLM / LM TRENDLINE THROUGH POPULATION MEANS


# ---------------------------------------------------------
# 1. Calculate one mean DDG per population within each gene
# ---------------------------------------------------------

pop_means <- df_plot %>%
  group_by(gene_variant, population_code) %>%
  summarise(
    pop_mean = mean(DDG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gene_variant) %>%
  arrange(pop_mean, .by_group = TRUE) %>%
  mutate(
    pop_index = row_number()
  ) %>%
  ungroup()

# ---------------------------------------------------------
# 2. Create the ordered factor used by ggplot
# ---------------------------------------------------------

pop_means <- pop_means %>%
  mutate(
    pop_ordered = reorder_within(
      population_code,
      -pop_mean,
      gene_variant
    )
  )

# ---------------------------------------------------------
# 3. Fit a linear model to the population means
# ---------------------------------------------------------

model_slopes <- pop_means %>%
  group_by(gene_variant) %>%
  summarise(
    model = list(lm(pop_mean ~ pop_index)),
    slope = coef(lm(pop_mean ~ pop_index))[2],
    .groups = "drop"
  ) %>%
  mutate(
    facet_label = paste0(
      gene_variant,
      ", AMECC: ",
      round(slope, 3)
    )
  )

# ---------------------------------------------------------
# 4. Generate predicted values for the trendline
# ---------------------------------------------------------

trendline <- pop_means %>%
  left_join(
    model_slopes %>% select(gene_variant, facet_label),
    by = "gene_variant"
  ) %>%
  group_by(gene_variant) %>%
  mutate(
    predicted = predict(
      lm(pop_mean ~ pop_index),
      newdata = data.frame(pop_index = pop_index)
    )
  ) %>%
  ungroup()

# ---------------------------------------------------------
# 5. Add facet labels to the raw data
# ---------------------------------------------------------

df_faceted <- df_plot %>%
  left_join(
    pop_means %>%
      select(
        gene_variant,
        population_code,
        pop_mean,
        pop_index,
        pop_ordered
      ),
    by = c("gene_variant", "population_code")
  ) %>%
  left_join(
    model_slopes %>%
      select(gene_variant, facet_label),
    by = "gene_variant"
  )

# ---------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------

ggplot(
  df_faceted,
  aes(
    x = pop_ordered,
    y = DDG,
    color = population_code
  )
) +
  
  # Raw individual haplotype points
  geom_point(
    position = position_jitter(width = 0.2),
    alpha = 0.4,
    size = 1.1
  ) +
  
  # Mean  SE
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    color = "black",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  
  # Population mean
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    color = "black",
    show.legend = FALSE
  ) +
  
  # -------------------------------------------------------
# Linear trend fitted to the POPULATION MEANS
# -------------------------------------------------------
geom_line(
  data = trendline,
  aes(
    x = pop_ordered,
    y = predicted,
    group = 1
  ),
  color = "blue",
  linewidth = 1,
  inherit.aes = FALSE
) +
  
  # Zero reference
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  
  # Facets
  facet_wrap(
    ~ facet_label,
    scales = "free",
    ncol = 2
  ) +
  
  scale_x_reordered() +
  
  labs(
    title = "DDG distributions with ancestry-dependent mutation energy change coefficient (AMECC)",
    subtitle = "Populations sorted high to low by mean DDG (Blue line = linear fit)",
    x = "Population Code",
    y = expression(paste(Delta, Delta, "G (energy change in REU)")),
    color = "Population"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),
    legend.position = "none"
  )

