library(dplyr)
library(ggplot2)
library(maps)
library(tidyr)
library(stringr)
library(ggrepel)

gene <- "MEFV"
setwd("~/git/MORPHA/")


########################################################################
# Merge IGSR sample metadata and protein haplotype quantification data #
########################################################################

# read IGSR sample metadata and clean up
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

# Merge the data, now in long format
samples_quant_long <- quant %>%
  left_join(samples, by = c("Sample_name" = "Sample name"))


##########################################################################################
# Pivot the data so that WT/pathogenic get their own columns and each row is a haplotype #
##########################################################################################

# Isolate the haplotype-specific WT energies as baseline
wt_baselines <- samples_quant_long %>%
  filter(Origin == "WT") %>%
  select(Haplotype, WT_Energy = Total_Energy)

# Process and reshape the data frame based strictly on ClinVarID
samples_quant <- samples_quant_long %>%
  # Step A: Filter down to just the pathogenic rows
  filter(Origin == "PG") %>%
  
  # Step B: Attach the correct WT background energy for this specific haplotype
  left_join(wt_baselines, by = "Haplotype") %>%
  
  # Step C: Calculate DDG using the individual's specific WT background
  mutate(DDG = Total_Energy - WT_Energy) %>%
  
  # Step D: Prefix the ClinVarID column headers so they are clear (e.g., ClinVar_917511)
  mutate(Column_Name = paste0("ClinVar_", ClinVarID)) %>%
  
  # Step E: Pivot wider so each ClinVar ID gets its own column containing the DDG
  pivot_wider(
    id_cols = c(Sample_name, Haplotype, WT_Energy, Sex, population_code, `Population name`),
    names_from = Column_Name,
    values_from = DDG
  )

# View the flattened haplotype-level data frame
head(samples_quant)


#############################################################
# Aggregate WT/pathogenic quantifications per subpopulation #
#############################################################

pop_quant <- samples_quant %>%
  group_by(population_code) %>%
  summarise(
    sample_size = n(),  # Total number of haplotypes in the population
    
    # Apply summary statistics across WT_Energy AND all ClinVar columns dynamically
    across(
      .cols = c(WT_Energy, starts_with("ClinVar_")),
      .fns = list(
        mean   = ~ mean(.x, na.rm = TRUE)
        #median = ~ median(.x, na.rm = TRUE),
        #sd     = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.fn}_{.col}"  # Formats output column names (e.g., mean_ClinVar_917511)
    ),
    .groups = "drop"
  )

# Preview the summary table
head(pop_quant)


################################################
# Map subpopulations to geographic coordinates #
################################################

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
    116.40, 138, 113.30, 100.80, 106.66,          # EAS
    10, 11.25, -2, 27, -3.70,                # EUR (CEU moved to Northwest Europe/France)
    4, 37.90, -16.20, -11.80, 8, -110.5, -59.5,  # AFR (ASW and ACB moved to West Africa near Nigeria/Ghana)
    -102, -66.60, -75.60, -77.04,                  # AMR (MXL moved to Central Mexico)
    72.57, 74.35, 90.40, 79.20, 80.70                # SAS (GIH, ITU, STU all in South Asia)
  ),
  latitude = c(
    39.90, 36.2048, 23.20, 22.00, 10.76,               # EAS
    54, 43.77, 54, 65, 40.41,               # EUR
    8.5, -0.10, 13.38, 8.48, 6.20, 37.5, 13.17,      # AFR
    22, 18.22, 6.25, -12.04,                      # AMR
    23.02, 31.4, 23.81, 14.50, 7.50                 # SAS
  )
)


###############################################################
# Merge with subpopulation aggregates and load world polygons #
###############################################################

map_data_points <- pop_quant %>%
  left_join(pop_coordinates, by = "population_code")
world_map <- map_data("world")
map_data_points$population_name_code <- paste0(map_data_points$population_name, " (",map_data_points$population_code,")")

######################
# Plot using ggplot2 #
######################


make_world_plot <- function(main_title, sub_title, legend_title, variable_name, save_file_name)
{
  
  col_oceans = "lightblue" #"#edf4f9"
  col_countries = "gray" #f2f2f2"
  col_borders = "black" # "#d9d9d9"
  col_legend_bg = "#ffffffcc"
  col_legend_border = "#d9d9d9"
  
  # legend tucked into the lower right corner
  plot_mean_WT_Energy <- ggplot() +
    # Base Layer: A subtle, neutral background of world borders
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
                 fill = col_countries, color = col_borders, linewidth = 0.2) +
    
    # Diversity Layer: Points representing distinct ethnic cohorts
    geom_point(data = map_data_points, 
               aes(x = longitude, y = latitude, color = population_name_code, size = {{variable_name}}),
               alpha = 0.9) +
    #geom_label(data = map_data_points, 
    #          aes(x = longitude, y = latitude, label = population_code), 
    #          size = 2, fontface = "bold", color = "black") +
    geom_label_repel(
      data = map_data_points, 
      aes(x = longitude, y = latitude, label = population_code), 
      size = 2.2, 
      fontface = "bold", 
      color = "black",
      fill = "#ffffffb3",          # Solid-ish white background box for readability
      
      # --- LINE FORCE PARAMETERS ---
      min.segment.length = 0,      # FIX: Forces lines to ALWAYS draw even if very short
      segment.size = 0.3,          # FIX: Replaced segment.linewidth to support your version
      segment.color = "gray20",    # Dark gray lines look cleaner than stark black lines
      
      # --- PADDING AND DISTANCE ---
      box.padding = 0.6,           # Pushes the labels a bit further out to give the lines room to breathe
      point.padding = 0.4,         # Keeps the line from touching the very center of the circle
      force = 8,                   # Increases the push away strength to display clear lines
      max.overlaps = Inf
    )+
    
    # Aesthetics and scaling
    scale_size_continuous(range = c(3, 8), name = legend_title) +
    scale_color_viridis_d(option = "turbo") +
    
    labs(
      title = main_title,
      subtitle = sub_title,
      x = "Longitude", y = "Latitude",
      color = "Specific Population / Ethnicity"
    ) +
    
    # Set explicit layout order for the legends ---
    # Mean energy gets order = 1 (Left), Population gets order = 2 (Right)
    guides(
      size = guide_legend(order = 1),
      color = guide_legend(ncol = 2, order = 2, override.aes = list(size = 3))
    ) +
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = col_oceans, color = "black"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.25, "cm"),
      legend.key.spacing.y = unit(0.02, "cm"),
      
      # Positioning the main container inside the lower right
      legend.position = "inside",              
      legend.position.inside = c(0.99, 0.01),   
      legend.justification = c(1, 0),          
      
      # Arrange distinct legend boxes horizontally side-by-side
      legend.box = "horizontal",               # Places size and color scales next to each other
      legend.box.just = "bottom",              # Aligns their bottom baselines cleanly
      
      # Styling background box
      legend.background = element_rect(fill = col_legend_bg, color =col_legend_border, linewidth = 0.2),
      legend.margin = margin(6, 6, 6, 6)       
    )
  
  ggsave(filename=save_file_name, plot = plot_mean_WT_Energy, bg = "white", width = 10, height = 6.5)
  
}

make_world_plot("1000 Genomes Project Sample Populations", "Visualizing sub-national regional and ethnic distribution", "Folding energy ΔG (kJ/mol)", mean_WT_Energy, "plot_mean_WT_Energy.png")
make_world_plot("bla", "bla", "Folding energy ΔΔG (kJ/mol)", mean_ClinVar_917511, "mean_ClinVar_917511.png")



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
ggplot(plot_data, aes(x = ClinVar_917511, fill = superpopulation_name, color = superpopulation_name)) +
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




