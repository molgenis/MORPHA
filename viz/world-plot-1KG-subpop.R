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

head(samples)
unique(samples$`Population code`)
unique(samples$`Population name`)

quant <- read.delim( paste0("results/", gene, "/", gene, "_protein_haplotypes_quant.tsv"),  sep = "\t", check.names = FALSE)
head(quant)



# clean data (dropping a single mixed ancestry individual)
samples_clean <- samples %>%
  rename(population_code = `Population code`) %>%
  filter(`Sample name` != "HG01783")


merged_df <- quant %>%
  left_join(samples_clean, by = c("Sample_name" = "Sample name"))

# aggregate metrics per population
pop_quant <- merged_df %>%
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
    "CHB", "JPT", "CHS", "CDX", "KHV", 
    "CEU", "TSI", "GBR", "FIN", "IBS", 
    "YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB", 
    "MXL", "PUR", "CLM", "PEL", 
    "GIH", "PJL", "BEB", "ITU", "STU"
  ),
  population_name = c(
    "Han Chinese (Beijing, North China)", "Japanese (Tokyo)", "Southern Han Chinese", "Chinese Dai (Xishuangbanna, Yunnan)", "Kinh (Vietnam)",
    "Utah Residents (NW European Ancestry)", "Toscani (Italy)", "British (UK)", "Finnish", "Iberian (Spain)",
    "Yoruba (Ibadan, Nigeria)", "Luhya (Webuye, Kenya)", "Gambian (Western Division)", "Mende (Sierra Leone)", "Esan (Nigeria)", "African Ancestry (SW USA)", "African Caribbean (Barbados)",
    "Mexican Ancestry (LA, USA)", "Puerto Rican", "Colombian (Medellin)", "Peruvian (Lima)",
    "Gujarati Indian (Houston, USA)", "Punjabi (Lahore, Pakistan)", "Bengali (Bangladesh)", "Indian Telugu (UK)", "Sri Lankan Tamil (UK)"
  ),
  # Approximate center-point Coordinates for the cohorts
  longitude = c(
    116.40, 139.70, 113.00, 100.80, 105.80,   # EAS
    -111.90, 11.25, -1.50, 24.90, -3.70,      # EUR
    3.90, 37.90, -16.60, -11.80, 6.00, -118.25, -59.50, # AFR
    -118.24, -66.60, -75.60, -77.04,          # AMR
    71.70, 74.35, 90.40, 80.60, 80.00         # SAS
  ),
  latitude = c(
    39.90, 35.68, 25.00, 22.00, 21.03,        # EAS
    40.76, 43.77, 52.50, 60.17, 40.41,        # EUR
    7.38, -0.10, 13.45, 8.48, 6.50, 34.05, 13.19,       # AFR
    34.05, 18.22, 6.25, -12.04,               # AMR
    22.30, 31.52, 23.81, 16.50, 9.66          # SAS
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
  geom_point(data = map_data_points, 
             aes(x = longitude, y = latitude, color = population_name, size = mean_energy),
             alpha = 0.8) +
  
  # Aesthetics and scaling
  scale_size_continuous(range = c(3, 8), name = "Mean energy") +
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
