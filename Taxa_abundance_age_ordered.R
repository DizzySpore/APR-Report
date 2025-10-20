# Bracken Data Processing and Visualization Script
# This script processes Bracken files, creates relative abundance visualizations
# with original sample names ordered by age from Malawi_Samples_age.csv, and uses a color
# palette for species.
# Load required libraries explicitly
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(patchwork)
# Set up for headless environment (no X11 display)
options(bitmapType = 'cairo')
if (!interactive()) {
  # Force cairo device for PNG output on headless systems
  if (capabilities("cairo")) {
    options(device = function(...) cairo_pdf(...))
  }
}
# ============================================================================
# DATA PROCESSING FUNCTIONS
# ============================================================================
# Function to process a single Bracken file
process_bracken_file <- function(file_path) {
  tryCatch({
    # Extract sample name from file path
    sample_name <- gsub("_prokaryote_bracken_abundance\\.txt$", "", basename(file_path))
    sample_name <- gsub("_prok-bracken_abundance\\.txt$", "", sample_name)
   
    # Read the data with more robust parameters
    data <- read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "", # Disable quote parsing
      fill = TRUE, # Handle uneven columns
      check.names = FALSE # Preserve column names
    )
   
    # Verify expected columns are present
    required_cols <- c("name", "taxonomy_lvl", "fraction_total_reads")
    if (!all(required_cols %in% colnames(data))) {
      warning(sprintf("Missing required columns in file: %s", basename(file_path)))
      return(NULL)
    }
   
    # Filter for species-level data and remove NA values
    data <- data %>%
      filter(
        taxonomy_lvl == "S",
        !is.na(fraction_total_reads)
      )
   
    # Add sample column
    data$sample <- sample_name
   
    # Reorder columns to have sample first
    data <- data[, c("sample", names(data)[names(data) != "sample"])]
   
    return(data)
   
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", basename(file_path), e$message))
    return(NULL)
  })
}
# ============================================================================
# COLOR PALETTE FUNCTION
# ============================================================================
# Function to create unified color palette with multiple options
create_unified_color_palette <- function(data, top_n = 20, palette_option = "qualitative") {
  # Get top species by total abundance
  top_species <- data %>%
    group_by(name_clean) %>%
    summarise(total_abundance = sum(fraction_total_reads)) %>%
    arrange(desc(total_abundance)) %>%
    head(top_n) %>%
    pull(name_clean)
  # Number of colors needed
  n_colors <- length(top_species)
  cat(sprintf("Selected %d top species for palette\n", n_colors))
  # Print the top species for debugging
  cat("Top species:\n")
  print(top_species)
  # Define different palette options
  palette <- switch(palette_option,
                    "default" = {
                      base_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "vibrant" = {
                      base_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "pastel" = {
                      base_colors <- c(brewer.pal(8, "Pastel1"), brewer.pal(8, "Pastel2"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "rainbow" = {
                      colorRampPalette(rainbow(10))(n_colors)
                    },
                    "qualitative" = {
                      base_colors <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set3"))
                      colorRampPalette(base_colors)(n_colors)
                    },
                    "spectral" = {
                      colorRampPalette(rev(brewer.pal(11, "Spectral")))(n_colors)
                    },
                    {
                      base_colors <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Accent"))
                      colorRampPalette(base_colors)(n_colors)
                    }
  )
  # Add gray for "Other" category
  palette <- c(palette, "gray70")
  names(palette) <- c(top_species, "Other")
  return(list(palette = palette, top_species = top_species))
}
# Function to create normalized abundance plots
create_abundance_plot <- function(processed_data, top_species, color_palette, title_prefix, interactive = FALSE, age_data) {
  top_n <- min(20, length(top_species))
  top_categories <- top_species[1:top_n]
  all_samples <- unique(processed_data$sample)
  
  # Get actual abundances for top species (sum if multiple entries, but unlikely)
  actual_tops <- processed_data %>%
    filter(name_clean %in% top_categories) %>%
    group_by(sample, name_clean) %>%
    summarise(abundance = sum(fraction_total_reads), .groups = "drop") %>%
    rename(species_category = name_clean)
  
  # Create full grid for top categories, filling missing with 0
  top_grid <- expand_grid(sample = all_samples, species_category = top_categories) %>%
    left_join(actual_tops, by = c("sample", "species_category")) %>%
    mutate(abundance = ifelse(is.na(abundance), 0, abundance))
  
  # Calculate sum of top abundances per sample
  sum_tops <- top_grid %>%
    group_by(sample) %>%
    summarise(sum_tops = sum(abundance), .groups = "drop")
  
  # Create Other data
  other_data <- expand_grid(sample = all_samples) %>%
    left_join(sum_tops, by = "sample") %>%
    mutate(
      species_category = "Other",
      abundance = pmax(1.0 - sum_tops, 0)
    )
  
  # Combine all data
  plot_data <- bind_rows(top_grid, other_data) %>%
    left_join(age_data, by = "sample") %>%
    mutate(
      sample = factor(sample, levels = age_data$sample[order(age_data$age)]),
      species_category = factor(species_category, levels = c(top_categories, "Other"))
    )
  
  # Debug print
  cat(sprintf("For %s plot: %d unique species categories (including Other)\n", title_prefix, nlevels(plot_data$species_category)))
  cat(sprintf("  Top species count: %d\n", top_n))
  
  # Create the plot with normalized values
  p <- ggplot(plot_data, aes(x = sample, y = abundance, fill = species_category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = color_palette, name = "Species", drop = FALSE) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = paste(title_prefix, ": Relative Abundance of Top", top_n, "Species Across Samples"),
      subtitle = "Remaining species are grouped as 'Other', samples ordered by age (months)",
      x = "Sample",
      y = "Relative Abundance (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))
  if (interactive) {
    # Ensure plotly is loaded before calling ggplotly
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package is not installed. Skipping interactive plot.")
      return(p)
    }
    p <- tryCatch({
      ggplotly(p, tooltip = c("x", "y", "fill"))
    }, error = function(e) {
      warning("Error in ggplotly: ", e$message, ". Returning static plot.")
      return(p)
    })
  }
  return(p)
}
# ============================================================================
# MAIN PROCESSING PIPELINE
# ============================================================================
cat("Starting Bracken data processing...\n")
# --- 1. PROCESS BRACKEN DATA ---
cat("\n=== 1. PROCESSING BRACKEN DATA ===\n")
bracken_files <- list.files(path = "kraken2_bracken_out", pattern = "*_prokaryote_bracken_abundance\\.txt$|*_prok-bracken_abundance\\.txt$", full.names = TRUE)
bracken_processed <- NULL
if (length(bracken_files) > 0) {
  cat("Found", length(bracken_files), "Bracken files.\n")
  all_bracken_data <- map_dfr(bracken_files, process_bracken_file)
  bracken_processed <- all_bracken_data %>%
    mutate(name_clean = gsub("_", " ", gsub("s__", "", sub(".*\\|", "", name)))) %>%
    filter(!is.na(fraction_total_reads) & name_clean != "") %>%
    group_by(sample) %>%
    mutate(fraction_total_reads = fraction_total_reads / sum(fraction_total_reads, na.rm = TRUE)) %>%
    ungroup()
  cat("- Bracken data processed and normalized to 100% per sample.\n")
  cat(sprintf("- Unique species in Bracken: %d\n", length(unique(bracken_processed$name_clean))))
  
  # Debug: Print first 10 Bracken sample names
  cat("\nFirst 10 Bracken sample names (sorted):\n")
  print(head(sort(unique(bracken_processed$sample)), 10))
} else {
  cat("No Bracken files found.\n")
}
# --- 2. LOAD AGE DATA AND PREPARE COLOR MAPPING ---
cat("\n=== 2. LOADING AGE DATA AND CREATING COLOR MAPPING ===\n")
# Safely combine sample names and species data
all_sample_names <- unique(c(
  if (!is.null(bracken_processed)) bracken_processed$sample else NULL
))
cat(sprintf("Total unique Bracken samples: %d\n", length(all_sample_names)))
# Load age data
age_data <- read.csv("Malawi_Samples_age.csv", stringsAsFactors = FALSE)
colnames(age_data) <- c("sample", "age")

# Debug: Print first 10 age file samples
cat("\nFirst 10 age file samples:\n")
print(head(age_data$sample, 10))

cat("Number of age samples before filter:", nrow(age_data), "\n")

age_data <- age_data %>%
  filter(sample %in% all_sample_names) %>%
  arrange(age)

cat("Number of age samples after filter:", nrow(age_data), "\n")

if (nrow(age_data) > 0) {
  cat("- Age data loaded for", nrow(age_data), "samples from 'Malawi_Samples_age.csv'.\n")
} else {
  cat("- No valid age data found.\n")
  stop("Age data is required for ordering samples.")
}
# Create a single color palette for all species
data_to_combine <- list()
if (!is.null(bracken_processed)) {
  data_to_combine$bracken <- select(bracken_processed, name_clean, fraction_total_reads)
}
all_species_data <- bind_rows(data_to_combine)
if (nrow(all_species_data) > 0) {
  # Filter to only species with positive total abundance
  all_species_data <- all_species_data %>%
    group_by(name_clean) %>%
    filter(sum(fraction_total_reads) > 0) %>%
    ungroup()
  cat(sprintf("- Total unique species with >0 abundance across all data: %d\n", length(unique(all_species_data$name_clean))))
  color_info <- create_unified_color_palette(all_species_data, top_n = 20, palette_option = "qualitative")
  unified_color_palette <- color_info$palette
  top_species_list <- color_info$top_species
} else {
  cat("- No species data to generate a color palette.\n")
  stop("No species data available for visualization.")
}
# ============================================================================
# CREATE AND SAVE VISUALIZATIONS
# ============================================================================
cat("\n=== 3. creating VISUALIZATIONS ===\n")
if (!is.null(bracken_processed) && nrow(bracken_processed) > 0) {
  cat("Creating abundance plots...\n")
  # Create static plot
  bracken_plot <- create_abundance_plot(
    processed_data = bracken_processed,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = "Bracken",
    age_data = age_data
  )
 
  # Save static PNG
  png("bracken_relative_abundance_stacked_bar.png",
      width = 16, height = 9,
      units = "in", res = 300, type = "cairo")
  print(bracken_plot)
  dev.off()
 
  # Create and save interactive plot
  bracken_interactive <- create_abundance_plot(
    processed_data = bracken_processed,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = "Bracken",
    interactive = TRUE,
    age_data = age_data
  )
  tryCatch({
    saveWidget(bracken_interactive, "bracken_abundance_interactive.html",
               selfcontained = TRUE)
    cat("- Interactive HTML saved successfully.\n")
  }, error = function(e) {
    warning("Failed to save interactive HTML (pandoc may be missing): ", e$message, "\nStatic PNG is available.\n")
  })
 
  cat("- Bracken plots saved: .png and (if possible) .html files\n")
} else {
  cat("- Insufficient data to create plots\n")
}
cat("\nScript execution completed successfully!\n")