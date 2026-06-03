library(dplyr)
library(ggplot2)
library(maps)

gene <- "MEFV"

setwd("~/git/MORPHA/")

samples <- read.delim("viz/igsr_samples.tsv", sep = "\t", check.names = FALSE)
samples$`Biosample ID` <- NULL
samples$`Data collections` <- NULL
samples$`Superpopulation code` <- NULL
samples$`Superpopulation name` <- NULL
samples$`Population elastic ID` <- NULL
samples <- samples %>%
  rename(population_code = `Population code`) %>%
  filter(`Sample name` != "HG01783") # drop 1 mixed ancestry individual

head(samples)
unique(samples$population_code)

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

# MERGE THE DATA
samples_quant <- quant %>%
  left_join(samples, by = c("Sample_name" = "Sample name"))

# aggregate metrics per population
pop_quant <- samples_quant %>%
  group_by(population_code) %>%
  summarise(
    sample_size = n(),                         # Total number of haplotypes
    mean_energy = mean(Total_Energy, na.rm = TRUE),
    median_energy = median(Total_Energy, na.rm = TRUE),
    sd_energy = sd(Total_Energy, na.rm = TRUE),
    .groups = "drop"                           # Keeps downstream data manipulation clean
  )

head(pop_quant)



# aggregrate data to get total sample sizes per group
#pop_quant <- samples_clean %>% count(population_code, name = "sample_size")


# 3. Build a precise Coordinate and Ethnic Name Lookup Table
# This breaks groups out by their specific sampling or historical locations
pop_coordinates <- tibble(
  population_code = c(
    "CHB", "JPT", "CHS", "CDX", "KHV",               # East Asian Ancestry (EAS)
    "CEU", "TSI", "GBR", "FIN", "IBS",               # European Ancestry (EUR)
    "YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB", # African Ancestry (AFR)
    "MXL", "PUR", "CLM", "PEL",                      # American Ancestry (AMR)
    "GIH", "PJL", "BEB", "ITU", "STU"                # South Asian Ancestry (SAS)
  ),
  population_name = c(
    "Han Chinese in Beijing, China", "Japanese in Tokyo, Japan", "Han Chinese South", "Chinese Dai in Xishuangbanna, China", "Kinh in Ho Chi Minh City, Vietnam", 
    "Utah residents (CEPH) with Northern and Western European ancestry", "Toscani in Italy", "British in England and Scotland", "Finnish in Finland", "Iberian populations in Spain",
    "Yoruba in Ibadan, Nigeria", "Luhya in Webuye, Kenya", "Gambian in Western Division, The Gambia - Mandinka", "Mende in Sierra Leone", "Esan in Nigeria", "African Ancestry in Southwest US", "African Caribbean in Barbados",
    "Mexican Ancestry in Los Angeles, California", "Puerto Rican in Puerto Rico", "Colombian in Medellin, Colombia", "Peruvian in Lima, Peru",
    "Gujarati Indians in Houston, TX", "Punjabi in Lahore, Pakistan", "Bengali in Bangladesh", "Indian Telugu in the UK", "Sri Lankan Tamil in the UK"
  ),
  longitude = c(
    116.40, 139.70, 113.30, 100.80, 106.66,          # EAS
    2.00, 11.25, -1.50, 24.90, -3.70,                # EUR (CEU moved to Northwest Europe/France)
    3.90, 37.90, -16.60, -11.80, 6.80, -2.00, 1.00,  # AFR (ASW and ACB moved to West Africa near Nigeria/Ghana)
    -99.13, -66.60, -75.60, -77.04,                  # AMR (MXL moved to Central Mexico)
    72.57, 74.35, 90.40, 79.20, 80.70                # SAS (GIH, ITU, STU all in South Asia)
  ),
  latitude = c(
    39.90, 35.68, 23.20, 22.00, 10.76,               # EAS
    47.00, 43.77, 52.50, 60.17, 40.41,               # EUR
    7.38, -0.10, 13.45, 8.48, 6.20, 5.00, 6.00,      # AFR
    19.43, 18.22, 6.25, -12.04,                      # AMR
    23.02, 31.52, 23.81, 14.50, 7.50                 # SAS
  )
)

# 4. Merge sample counts with geographic points
map_data_points <- pop_quant %>%
  left_join(pop_coordinates, by = "population_code")

# 5. Plot using ggplot2
world_map <- map_data("world")
ggplot() +
  # Base Layer: A subtle, neutral background of world borders
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "#f2f2f2", color = "#d9d9d9", linewidth = 0.2) +
  
  # Diversity Layer: Points representing distinct ethnic cohorts
  # FIX 1: Added position_jitter so overlapping points nudge away from each other slightly
  geom_point(data = map_data_points, 
             aes(x = longitude, y = latitude, color = population_name, size = mean_energy),
             alpha = 0.9,
             position = position_jitter(width = 2.5, height = 2.5, seed = 42)) +
  
  # Aesthetics and scaling
  scale_size_continuous(range = c(3, 8), name = "Mean energy") +
  
  # FIX 2: Use the "Turbo" palette which maxes out distinct hues across 26 categories
  scale_color_viridis_d(option = "turbo") +
  
  labs(
    title = "1000 Genomes Project Sample Populations",
    subtitle = "Visualizing sub-national regional and ethnic distribution",
    x = "Longitude", y = "Latitude",
    color = "Specific Population / Ethnicity"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#edf4f9", color = NA),
    legend.position = "right",
    legend.text = element_text(size = 8)
  )



# ======= DENSITY PLOT ======= #



pop_lookup <- tibble(
  population_code = c(
    "CHB", "JPT", "CHS", "CDX", "KHV",               # EAS
    "CEU", "TSI", "GBR", "FIN", "IBS",               # EUR
    "YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB", # AFR
    "MXL", "PUR", "CLM", "PEL",                      # AMR
    "GIH", "PJL", "BEB", "ITU", "STU"                # SAS
  ),
  superpopulation_code = c(
    rep("EAS", 5),
    rep("EUR", 5),
    rep("AFR", 7),
    rep("AMR", 4),
    rep("SAS", 5)
  ),
  superpopulation_name = c(
    rep("East Asian Ancestry", 5),
    rep("European Ancestry", 5),
    rep("African Ancestry", 7),
    rep("American Ancestry", 4),
    rep("South Asian Ancestry", 5)
  )
)

# ==============================================================================
# STEP 2: Map the Lookup Table to your Raw, Unaggregated Data
# ==============================================================================
# Assuming 'merged_df' is your haplotype-level data frame containing 'Total_Energy' 
# and 'population_code' columns.
plot_data <- samples_quant %>%
  left_join(pop_lookup, by = "population_code")

# ==============================================================================
# STEP 3: Generate the Comparative Density Plot
# ==============================================================================
ggplot(plot_data, aes(x = Total_Energy, fill = superpopulation_name, color = superpopulation_name)) +
  # Add geom_density with transparency (alpha) so overlapping distributions are visible
  geom_density(alpha = 0.35, linewidth = 0.7) +
  
  # Standard 1000 Genomes-style clean color palette
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  
  # Labeling
  labs(
    title = "Distribution of MEFV Structural Total Energy",
    subtitle = "Global comparison across 1000 Genomes continental superpopulations",
    x = "Total Energy (kcal/mol)",
    y = "Density",
    fill = "Superpopulation",
    color = "Superpopulation"
  ) +
  
  # Professional theme configurations
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 11),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )