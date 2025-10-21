# Define library path
lib_path <- "/rds/projects/h/hallly-microbiome/software/R_tools/R_LIB/lxm697_lib/"

# Load required libraries
library(ggplot2, lib.loc = lib_path)
library(sf, lib.loc = lib_path)
library(rnaturalearth, lib.loc = lib_path)
library(rnaturalearthdata, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(grid, lib.loc = lib_path)
library(gridExtra, lib.loc = lib_path)

# Read the CSV file
sample_coords <- read.csv("Malawi_metadata_test_import_gps_coordinates.csv", stringsAsFactors = FALSE)

# Filter out rows with missing or invalid coordinates
sample_coords <- sample_coords %>% 
  filter(!is.na(gps_longitude) & !is.na(gps_latitude) & gps_longitude != "" & gps_latitude != "")

# Convert sample coordinates to sf object
samples_sf <- st_as_sf(sample_coords, coords = c("gps_longitude", "gps_latitude"), crs = 4326)

# Get map data for Malawi
malawi <- ne_countries(country = "Malawi", scale = "medium", returnclass = "sf")

# Get map data for Africa (for inset map)
africa <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")

# Create main map of Malawi with sample points
main_map <- ggplot() +
  geom_sf(data = malawi, fill = "lightgray", color = "black") +
  geom_sf(data = samples_sf, color = "red", size = 0.5, shape = 21, fill = "red") +
  theme_minimal() +
  labs(title = "Sample Locations in Malawi", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(32, 36, by = 2))  # Reduce frequency of longitude labels

# Create inset map showing Malawi in Africa
inset_map <- ggplot() +
  geom_sf(data = africa, fill = "lightgray", color = "black") +
  geom_sf(data = malawi, fill = "red", color = "black") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  coord_sf(xlim = c(-20, 55), ylim = c(-35, 37))  # Focus on Africa

# Combine maps using annotation_custom (move inset to top-right corner)
inset_grob <- ggplotGrob(inset_map)

final_plot <- main_map +
  annotation_custom(
    grob = inset_grob,
    xmin = 36,  # Longitude just beyond Malawi's eastern boundary
    xmax = 40,  # Adjust to fit inset width
    ymin = -15,  # Latitude just above Malawi's northern boundary
    ymax = -11   # Adjust to fit inset height
  )

# Display the plot
print(final_plot)

# Save the plot
ggsave("malawi_samples_map.png", final_plot, width = 8, height = 6, dpi = 300)