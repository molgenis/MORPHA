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
wt_energy_counts <- wt_baselines %>% count(WT_Energy, name = "Count", sort = TRUE)
head(wt_energy_counts)
top_ref_energy <- wt_energy_counts[[1]][1]
# ALTERNTIVE: Overall mean WT delta G energy across all haplotypes
mean_WT_energy <- mean(wt_baselines$WT_Energy, na.rm = TRUE)
# assign it here
ref_energy <- top_ref_energy

# Process and reshape the data frame based strictly on VariantID
samples_quant <- samples_quant_long %>%
  # Filter down to just the pathogenic rows (.. including Benign ..)
  filter(Origin == "PG") %>%
  
  # Attach the correct WT background energy for this specific haplotype
  left_join(wt_baselines, by = "Haplotype") %>%
  
  # Calculate DDG using the individual's specific WT background
  mutate(DDG = Total_Energy - WT_Energy) %>%
  
  # WT deviation from the global mean WT delta G
  mutate(WT_diff = WT_Energy - ref_energy) %>%
  
  # Prefix the ClinVarID column headers so they are clear (e.g., ClinVar_917511)
  mutate(Column_Name = paste0("Variant_", VariantID)) %>%
  
  # Pivot wider so each ClinVar ID gets its own column containing the DDG
  pivot_wider(
    id_cols = c(Sample_name, Haplotype, WT_Energy, WT_diff, Sex, population_code, `Population name`),
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
      .cols = c(WT_diff, starts_with("Variant_")), # used to be 'WT_Energy', now 'WT_diff' to be comparible to DDG values
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
      aes(x = longitude, y = latitude, label = paste0(population_code, " ", round({{variable_name}}))), 
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

# plot on world map!
#make_world_plot("1000 Genomes Project Sample Populations", "Visualizing sub-national regional and ethnic distribution", "Folding energy ΔG (kJ/mol)", mean_WT_Energy, "plot_mean_WT_Energy.png")


# ======= SUPERPOPULATION MAPPING ======= #

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



# ======= HEATMAP ======= #

# 1. Isolate the population code and the specific 'mean' columns you want
library(pheatmap)
variant_columns <- paste0("mean_Variant_", variant_info$id)
heatmap_data <- pop_quant %>%
  select(
    population_code,
    mean_WT_diff, # used to be 'mean_WT_Energy', now 'mean_WT_diff' to be comparible to DDG values
    all_of(variant_columns)
  )

# Convert the data frame into a matrix with population_code as row names
heatmap_matrix <- heatmap_data %>%
  tibble::column_to_rownames("population_code") %>%
  as.matrix()

# Match matrix columns to variant metadata
column_metadata <- tibble(
  matrix_column = colnames(heatmap_matrix)
) %>%
  mutate(
    VariantID = str_remove(matrix_column, "^mean_Variant_")
  ) %>%
  left_join(
    variant_info %>%
      select(id, label, classification),
    by = c("VariantID" = "id")
  ) %>%
  mutate(
    # Special handling for the WT column
    label = if_else(
      matrix_column == "mean_WT_diff",
      "Δhg38",
      coalesce(label, VariantID)
    ),
    classification = if_else(
      matrix_column == "mean_WT_diff",
      "Δhg38",
      classification
    )
  )

# Column annotation must have matrix column names as row names
annotation_col <- column_metadata %>%
  select(matrix_column, Classification = classification) %>%
  tibble::column_to_rownames("matrix_column") %>%
  as.data.frame()


# PREPARE THE METADATA ANNOTATION DATA FRAME
# Filter lookup down to only the populations currently present in your matrix rows
annotation_df <- pop_lookup %>%
  filter(population_code %in% rownames(heatmap_matrix)) %>%
  select(population_code, Superpopulation = superpopulation_code) %>%
  # pheatmap strictly requires metadata to be a data frame with matching row names
  tibble::column_to_rownames("population_code") %>%
  as.data.frame()

# DEFINE EXPLICIT COLORS FOR THE SUPERPOPULATIONS and Pathogenic/Benign
annotation_colors <- list(
  Superpopulation = c(
    EAS = "#7FC97F",
    EUR = "#BEAED4",
    AFR = "#FDC086",
    AMR = "#FFFF99",
    SAS = "#386CB0"
  ),
  
  Classification = c(
    `Δhg38` = "grey60",
    Benign = "#2CA25F",
    Pathogenic = "#DE2D26"
  )
)

# compute breaks
max_abs <- max(
  abs(heatmap_matrix),
  na.rm = TRUE
)
heatmap_breaks <- seq(
  -max_abs,
  max_abs,
  length.out = 51
)

# GENERATE ANNOTATED HEATMAP
pheatmap(
  mat = heatmap_matrix,
  scale = "none", # used to be column for within-variant compare, but can now compare all
  
  # Clustering logic
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  # Row/column Annotations, colors, labels
  annotation_row = annotation_df,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  labels_col = column_metadata$label,
  
  # Cell numbers display logic (Original raw values)
  display_numbers = round(heatmap_matrix, 2),
  number_color = "black",
  fontsize_number = 8,
  
  # Aesthetics & Background configuration
  color = colorRampPalette(c("#313695", "#f7f7f7", "#a50026"))(50),
  breaks = heatmap_breaks,
  
  border_color = "white",       
  background_color = "white",   
  
  # Labels and titles
  main = "Population vs Mutation Effects\n(Heatmap: column-wise Z-score | Text: Raw DDG)",
  # "3201 individuals aggregated at haplotype level 6402 per sub-population"
  angle_col = 45,
  
  # filename and dimensions here
  filename = "superpop_mutation_heatmap.png",
  width = 10,
  height = 6.5
)









# =========================
# INDIVIDUAL-LEVEL HEATMAP
# =========================
# Link the raw haplotype dataset to the ancestry lookup table
indiv_data <- samples_quant %>%
  left_join(pop_lookup, by = "population_code") %>%
  filter(!is.na(superpopulation_code)) %>%
  # Create a unique row identifier combining Sample and Haplotype designation
  mutate(Haplotype_ID = paste0(Sample_name, "_", str_extract(Haplotype, "h[12]$")))

# Convert to a pure numeric matrix for the heatmap grid
indiv_matrix <- indiv_data %>%
  select(
    Haplotype_ID,
    WT_Energy,
    starts_with("ClinVar_")
  ) %>%
  tibble::column_to_rownames("Haplotype_ID") %>%
  as.matrix()
# Clean up column names for display
colnames(indiv_matrix) <- gsub("ClinVar_", "ClinVar ", colnames(indiv_matrix))
# This creates a tracking bar mapping every individual haplotype to its ancestry
indiv_annotation <- indiv_data %>%
  select(Haplotype_ID, Superpopulation = superpopulation_code) %>%
  tibble::column_to_rownames("Haplotype_ID") %>%
  as.data.frame()
# Explicit colors for the ancestry tracker (matching your previous plots)
annotation_colors <- list(
  Superpopulation = c(
    EAS = "#7FC97F", # Soft Green
    EUR = "#BEAED4", # Soft Purple
    AFR = "#FDC086", # Soft Orange
    AMR = "#FFFF99", # Soft Yellow
    SAS = "#386CB0"  # Soft Blue
  )
)
pheatmap(
  mat = indiv_matrix,
  scale = "column", 
  
  # Clustering configuration
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  # Ancestry metadata tracking
  annotation_row = indiv_annotation,
  annotation_colors = annotation_colors,
  
  # Hide unreadable labels
  show_rownames = FALSE,         
  show_colnames = TRUE,
  cluster_rows = TRUE,          
  display_numbers = FALSE,
  
  # 1. Give NA cells a definitive color so they act as a clean canvas
  na_col = "#EAEAEA", 
  
  # 2. Cap the gradient scale dynamically to force vibrant contrast
  # This stops extreme outliers (-20 or +20) from washing out the entire plot
  breaks = seq(-3, 3, length.out = 51), 
  
  # 3. Clean color ramp running from deep variant-suppressed to highly destabilized
  color = colorRampPalette(c("#4575b4", "#f7f7f7", "#d73027"))(50), 
  
  # 4. Remove cell borders entirely so the data smoothly blends together
  border_color = NA,            
  background_color = "white",   
  
  # Labels and saving parameters
  main = "Individual Haplotype Variant Profiles\n(Background Grid: Non-Carriers | Colored Blocks: Column Z-Score)",
  angle_col = 45,
  
  filename = "individual_haplotype_heatmap_fixed.png",
  width = 8,
  height = 10                  
)




# ======= HEATMAP: Haplotype Counts by WT Energy Level for subpopulations ======= #

# 1. Calculate haplotype counts per subpopulation at each WT Energy level
heatmap_data <- samples_quant %>%
  group_by(population_code, WT_Energy) %>%
  summarise(haplotype_count = n(), .groups = "drop") %>%
  # Pivot wider so WT energy levels become columns (X-axis)
  pivot_wider(
    names_from = WT_Energy, 
    values_from = haplotype_count, 
    values_fill = 0
  )
# 2. Convert the data frame into a matrix with population_code as row names
heatmap_matrix <- heatmap_data %>%
  tibble::column_to_rownames("population_code") %>%
  as.matrix()
# 3. Sort column names numerically so the X-axis scales cleanly if needed
heatmap_matrix <- heatmap_matrix[, order(as.numeric(colnames(heatmap_matrix)))]
# PREPARE THE METADATA ANNOTATION DATA FRAME
annotation_df <- pop_lookup %>%
  filter(population_code %in% rownames(heatmap_matrix)) %>%
  select(population_code, Superpopulation = superpopulation_code) %>%
  tibble::column_to_rownames("population_code") %>%
  as.data.frame()
# DEFINE EXPLICIT COLORS FOR THE SUPERPOPULATIONS
annotation_colors <- list(
  Superpopulation = c(
    EAS = "#7FC97F", # Soft Green
    EUR = "#BEAED4", # Soft Purple
    AFR = "#FDC086", # Soft Orange
    AMR = "#FFFF99", # Soft Yellow
    SAS = "#386CB0"  # Soft Blue
  )
)
# 1. Find the maximum raw value in your matrix to anchor the top of the scale
max_val <- max(heatmap_matrix, na.rm = TRUE)
palette_length <- 100
# 2. Build a smooth, expanded warm color palette
heatmap_colors <- colorRampPalette(c("#ffffff", "#ffe5cc", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#7f0000"))(palette_length)
# 3. Create non-linear breaks using a square-root curve.
# This forces multiple distinct color changes down between 0 and 10, 
# while letting the huge 300+ counts share the deep red colors at the top.
heatmap_breaks <- seq(0, max_val^(1/2), length.out = palette_length + 1)^2
# GENERATE CLEANED COUNT HEATMAP WITH RAW VALUES AND COMPRESSED COLOR SCALE
pheatmap(
  mat = heatmap_matrix,   # Using the raw, untouched matrix values
  scale = "none",         # No mathematical transformation of the data
  
  # Clustering logic (retained per your preference)
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  # Row Annotations & Colors
  annotation_row = annotation_df,
  annotation_colors = annotation_colors,
  
  # Inject the custom color-compression arrays
  color = heatmap_colors,
  breaks = heatmap_breaks,
  
  # Structural formatting
  display_numbers = FALSE,
  border_color = NA,  
  
  # Labels and titles
  main = "Haplotype Count per Subpopulation by WT Energy Level\n(Color Scale Compressed to Highlight Lower Ranges 1-10)",
  
  angle_col = 90,
  fontsize_col = 6,   
  fontsize_row = 9,   
  
  # Output file settings
  filename = "wt_energy_counts_heatmap_compressed.png",
  width = 12,
  height = 7
)



# ======= HEATMAP: Haplotype Counts by WT Energy Level for superpopulations ======= #

# 1. Map subpopulations to superpopulations, then count haplotypes
heatmap_data <- samples_quant %>%
  left_join(pop_lookup, by = "population_code") %>%
  group_by(superpopulation_name, WT_Energy) %>% # Using full names for cleaner row labels
  summarise(haplotype_count = n(), .groups = "drop") %>%
  # Pivot wider so WT energy levels become columns (X-axis)
  pivot_wider(
    names_from = WT_Energy, 
    values_from = haplotype_count, 
    values_fill = 0
  )

# 2. Convert into a matrix with superpopulation_name as row names
heatmap_matrix <- heatmap_data %>%
  tibble::column_to_rownames("superpopulation_name") %>%
  as.matrix()

# 3. Sort column names numerically (X-axis)
heatmap_matrix <- heatmap_matrix[, order(as.numeric(colnames(heatmap_matrix)))]

# Find the new maximum value (binned superpopulation numbers will be larger)
max_val <- max(heatmap_matrix, na.rm = TRUE)
palette_length <- 100

# Expanded warm color palette
heatmap_colors <- colorRampPalette(c("#ffffff", "#ffe5cc", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#7f0000"))(palette_length)

# Dynamic square-root breaks adjusted to the new superpopulation data volume
heatmap_breaks <- seq(0, max_val^(1/2), length.out = palette_length + 1)^2

# GENERATE SUPERPOPULATION COUNT HEATMAP
pheatmap(
  mat = heatmap_matrix,   # Matrix now has 5 rows instead of 26
  scale = "none",         
  
  # Clustering logic
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  # Removed annotation_row & annotation_colors because rows are now the superpopulations!
  
  # Inject the custom color-compression arrays
  color = heatmap_colors,
  breaks = heatmap_breaks,
  
  # Structural formatting
  display_numbers = FALSE,
  border_color = NA,  
  
  # Labels and titles
  main = "Haplotype Count per Superpopulation by WT Energy Level\n(Color Scale Compressed to Highlight Lower Ranges)",
  
  angle_col = 90,
  fontsize_col = 6,   
  fontsize_row = 10,   # Boosted row font size slightly since there is more room now
  
  # Output file settings (Height reduced to 4.5 because 5 rows need less vertical space)
  filename = "superpop_energy_counts_heatmap.png",
  width = 12,
  height = 4.5
)
