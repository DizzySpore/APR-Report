# Load required libraries
library(WDI)
library(dplyr)
library(countrycode)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(scales)
library(viridis)

# Create output directory
if (!dir.exists("feeding_maps_output")) {
  dir.create("feeding_maps_output")
}

# ============================================================================
# SETUP: Load world map from Natural Earth
# ============================================================================
cat("Loading world map from Natural Earth...\n")

world <- ne_countries(scale = "medium", returnclass = "sf")
cat(sprintf("Loaded map with %d countries/territories\n", nrow(world)))

# Initialize year variables
exclusive_bf_year <- NULL
early_bf_year <- NULL
continued_bf_1yr_year <- NULL
continued_bf_2yr_year <- NULL
bottle_feeding_year <- NULL

# Aggregate codes to exclude
aggregate_codes <- c("1A", "1W", "4E", "7E", "8S", "B8", "EU", "F1", "JG", 
                     "OE", "S1", "S2", "S3", "S4", "T2", "T3", "T4", "T5", 
                     "T6", "T7", "V1", "V2", "V3", "V4", "XC", "XD", "XE", 
                     "XF", "XG", "XH", "XI", "XJ", "XK", "XL", "XM", "XN", 
                     "XO", "XP", "XQ", "XT", "XU", "XY", "Z4", "Z7", "ZF", 
                     "ZG", "ZH", "ZI", "ZJ", "ZQ", "ZT")

# ============================================================================
# STEP 1: Download Data
# ============================================================================
cat("\nDownloading exclusive breastfeeding data...\n")

exclusive_bf_df <- tryCatch({
  bf_data <- WDI(
    country = "all", 
    indicator = "SH.STA.BFED.ZS",  # Exclusive breastfeeding (% of children under 6 months)
    start = 2000, 
    end = 2025,
    extra = FALSE
  )
  
  cat(sprintf("  Downloaded: %d observations\n", nrow(bf_data)))
  
  df <- bf_data %>%
    filter(!iso2c %in% aggregate_codes) %>%
    group_by(iso2c) %>%
    filter(!is.na(SH.STA.BFED.ZS)) %>%
    filter(year == max(year)) %>%
    ungroup() %>%
    mutate(iso_a3 = countrycode(iso2c, "iso2c", "iso3c", warn = FALSE)) %>%
    filter(!is.na(iso_a3), !is.na(SH.STA.BFED.ZS)) %>%
    select(iso_a3, exclusive_bf_rate = SH.STA.BFED.ZS, year)
  
  cat(sprintf("  Processed: %d countries\n", nrow(df)))
  cat(sprintf("  Matched to map: %d countries\n", sum(df$iso_a3 %in% world$iso_a3_eh)))
  
  # Show which regions are missing data
  total_countries <- world %>% filter(!is.na(iso_a3_eh)) %>% nrow()
  missing_count <- total_countries - sum(df$iso_a3 %in% world$iso_a3_eh)
  cat(sprintf("  Countries without data: %d out of %d mappable countries\n", missing_count, total_countries))
  
  # Show example missing countries
  missing_countries <- world %>% 
    filter(!iso_a3_eh %in% df$iso_a3) %>%
    filter(!is.na(name)) %>%
    select(name, iso_a3_eh) %>%
    head(10)
  
  if (nrow(missing_countries) > 0) {
    cat("  Example countries without data:\n")
    for (i in 1:min(5, nrow(missing_countries))) {
      cat(sprintf("    - %s (%s)\n", missing_countries$name[i], missing_countries$iso_a3_eh[i]))
    }
  }
  
  if (nrow(df) > 0) {
    cat(sprintf("  Range: %.1f%% to %.1f%%\n", min(df$exclusive_bf_rate), max(df$exclusive_bf_rate)))
    cat(sprintf("  Latest year in data: %s\n", max(df$year)))
  }
  
  exclusive_bf_year <- max(df$year)
  df %>% select(iso_a3, exclusive_bf_rate)
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  data.frame()
})

# ============================================================================
# CREATE MAPS
# ============================================================================

# MAP 1: Exclusive Breastfeeding Rates
if (nrow(exclusive_bf_df) > 0) {
  cat("\n--- Creating Exclusive Breastfeeding Map ---\n")
  tryCatch({
    map_data <- world %>%
      left_join(exclusive_bf_df, by = c("iso_a3_eh" = "iso_a3"))
    
    p1 <- ggplot(data = map_data) +
      geom_sf(aes(fill = exclusive_bf_rate), color = "white", size = 0.1) +
      scale_fill_viridis_c(
        option = "viridis",
        name = "Exclusive BF Rate (%)",
        na.value = "grey90",
        direction = -1,
        labels = function(x) sprintf("%.0f", x)
      ) +
      labs(
        title = "Exclusive Breastfeeding Rates by Country",
        subtitle = sprintf("Percentage of infants 0-5 months exclusively breastfed (Data up to %s)", 
                          ifelse(!is.null(exclusive_bf_year), exclusive_bf_year, "latest available year"))
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
    
    print(p1)
    ggsave("feeding_maps_output/01_exclusive_breastfeeding.png", p1, width = 14, height = 8, dpi = 300)
    cat("✓ Saved to: feeding_maps_output/01_exclusive_breastfeeding.png\n")
    cat(sprintf("  Mapped: %d countries with data\n", sum(!is.na(map_data$exclusive_bf_rate))))
  }, error = function(e) {
    cat("✗ Error creating map:", e$message, "\n")
  })
}

cat("\n======================================\n")
cat("SCRIPT COMPLETED\n")
cat("======================================\n")
cat("Maps saved to: feeding_maps_output/\n")