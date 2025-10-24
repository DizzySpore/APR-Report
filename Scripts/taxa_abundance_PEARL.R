# Bracken Data Processing and Visualization Script
# This script processes Bracken files for Prokaryote, Archaea, and Eukaryote, creates relative abundance visualizations
# with original sample names ordered alphabetically, and uses a color palette for species.
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
# Function to process a single Prokaryote Bracken file
process_prok_bracken_file <- function(file_path) {
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
# Function to process a single Archaea Bracken file
process_arc_bracken_file <- function(file_path) {
  tryCatch({
    # Extract sample name from file path
    sample_name <- gsub("_arc-test_bracken_abundance\\.txt$", "", basename(file_path))
   
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
# Function to process a single Eukaryote Bracken file
process_euk_bracken_file <- function(file_path) {
  tryCatch({
    # Extract sample name from file path
    sample_name <- gsub("_euk-fin_bracken_abundance\\.txt$", "", basename(file_path))
   
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
# Function to create normalized abundance plots (stacked bar) without age-based separators
create_abundance_plot <- function(processed_data, top_species, color_palette, title_prefix, interactive = FALSE) {
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
    mutate(
      sample = factor(sample, levels = sort(all_samples)),
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
      subtitle = paste("Remaining species grouped as 'Other'; samples ordered alphabetically"),
      x = "Sample",
      y = "Relative Abundance (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
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
cat("Starting Bracken data processing for Prokaryote, Archaea, and Eukaryote...\n")

# --- PROCESS PROKARYOTE DATA ---
cat("\n=== 1a. PROCESSING PROKARYOTE BRACKEN DATA ===\n")
prok_files <- list.files(path = "kraken2_bracken_out", pattern = "*_prokaryote_bracken_abundance\\.txt$|*_prok-bracken_abundance\\.txt$", full.names = TRUE)
prok_processed <- NULL
if (length(prok_files) > 0) {
  cat("Found", length(prok_files), "Prokaryote Bracken files.\n")
  all_prok_data <- map_dfr(prok_files, process_prok_bracken_file)
  prok_processed <- all_prok_data %>%
    mutate(name_clean = gsub("_", " ", gsub("s__", "", sub(".*\\|", "", name)))) %>%
    filter(!is.na(fraction_total_reads) & name_clean != "") %>%
    group_by(sample) %>%
    mutate(fraction_total_reads = fraction_total_reads / sum(fraction_total_reads, na.rm = TRUE)) %>%
    ungroup()
  cat("- Prokaryote Bracken data processed and normalized to 100% per sample.\n")
  cat(sprintf("- Unique species in Prokaryote Bracken: %d\n", length(unique(prok_processed$name_clean))))
} else {
  cat("No Prokaryote Bracken files found.\n")
}

# --- PROCESS ARCHAEA DATA ---
cat("\n=== 1b. PROCESSING ARCHAEA BRACKEN DATA ===\n")
arc_files <- list.files(path = "kraken2_bracken_out", pattern = "*_arc-test_bracken_abundance\\.txt$", full.names = TRUE)
arc_processed <- NULL
if (length(arc_files) > 0) {
  cat("Found", length(arc_files), "Archaea Bracken files.\n")
  all_arc_data <- map_dfr(arc_files, process_arc_bracken_file)
  arc_processed <- all_arc_data %>%
    mutate(name_clean = gsub("_", " ", gsub("s__", "", sub(".*\\|", "", name)))) %>%
    filter(!is.na(fraction_total_reads) & name_clean != "") %>%
    group_by(sample) %>%
    mutate(fraction_total_reads = fraction_total_reads / sum(fraction_total_reads, na.rm = TRUE)) %>%
    ungroup()
  cat("- Archaea Bracken data processed and normalized to 100% per sample.\n")
  cat(sprintf("- Unique species in Archaea Bracken: %d\n", length(unique(arc_processed$name_clean))))
} else {
  cat("No Archaea Bracken files found.\n")
}

# --- PROCESS EUKARYOTE DATA ---
cat("\n=== 1c. PROCESSING EUKARYOTE BRACKEN DATA ===\n")
euk_files <- list.files(path = "kraken2_bracken_out", pattern = "*_euk-fin_bracken_abundance\\.txt$", full.names = TRUE)
euk_processed <- NULL
if (length(euk_files) > 0) {
  cat("Found", length(euk_files), "Eukaryote Bracken files.\n")
  all_euk_data <- map_dfr(euk_files, process_euk_bracken_file)
  euk_processed <- all_euk_data %>%
    mutate(name_clean = gsub("_", " ", gsub("s__", "", sub(".*\\|", "", name)))) %>%
    filter(!is.na(fraction_total_reads) & name_clean != "") %>%
    group_by(sample) %>%
    mutate(fraction_total_reads = fraction_total_reads / sum(fraction_total_reads, na.rm = TRUE)) %>%
    ungroup()
  cat("- Eukaryote Bracken data processed and normalized to 100% per sample.\n")
  cat(sprintf("- Unique species in Eukaryote Bracken: %d\n", length(unique(euk_processed$name_clean))))
} else {
  cat("No Eukaryote Bracken files found.\n")
}

# --- 2. PREPARE DATASETS ---
cat("\n=== 2. PREPARING DATASETS ===\n")

# Function to prepare a dataset
prepare_dataset <- function(processed_data, dataset_name) {
  if (is.null(processed_data) || nrow(processed_data) == 0) {
    cat(sprintf("- No %s data to prepare.\n", dataset_name))
    return(NULL)
  }
  
  cat(sprintf("Number of %s samples: %d\n", dataset_name, length(unique(processed_data$sample))))
  
  return(list(data = processed_data))
}

# Process each dataset
prok_with_data <- prepare_dataset(prok_processed, "Prokaryote")
arc_with_data <- prepare_dataset(arc_processed, "Archaea")
euk_with_data <- prepare_dataset(euk_processed, "Eukaryote")

# --- 3. CREATE VISUALIZATIONS FOR EACH DATASET ---
cat("\n=== 3. CREATING VISUALIZATIONS FOR EACH DATASET ===\n")

# Function to create and save plots for a dataset
create_and_save_plots <- function(data_with_prep, dataset_name, suffix) {
  if (is.null(data_with_prep)) {
    cat(sprintf("- Skipping plots for %s (no data).\n", dataset_name))
    return()
  }
  
  processed_data <- data_with_prep$data
  
  # Create color palette for this dataset
  all_species_data <- select(processed_data, name_clean, fraction_total_reads)
  all_species_data <- all_species_data %>%
    group_by(name_clean) %>%
    filter(sum(fraction_total_reads) > 0) %>%
    ungroup()
  cat(sprintf("- Total unique species with >0 abundance in %s: %d\n", dataset_name, length(unique(all_species_data$name_clean))))
  color_info <- create_unified_color_palette(all_species_data, top_n = 20, palette_option = "qualitative")
  unified_color_palette <- color_info$palette
  top_species_list <- color_info$top_species
  
  cat(sprintf("Creating %s abundance plots...\n", dataset_name))
  
  # Create static bar plot
  bar_plot <- create_abundance_plot(
    processed_data = processed_data,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = paste(dataset_name, "Bracken")
  )
 
  # Save static PNG for bar plot
  png_filename <- paste0("bracken_", suffix, "_relative_abundance_stacked_bar.png")
  png(png_filename,
      width = 16, height = 9,
      units = "in", res = 300, type = "cairo")
  print(bar_plot)
  dev.off()
 

 
  # Create and save interactive bar plot
  interactive_plot <- create_abundance_plot(
    processed_data = processed_data,
    top_species = top_species_list,
    color_palette = unified_color_palette,
    title_prefix = paste(dataset_name, "Bracken"),
    interactive = TRUE
  )
  interactive_filename <- paste0("bracken_", suffix, "_abundance_interactive.html")
  tryCatch({
    saveWidget(interactive_plot, interactive_filename,
               selfcontained = TRUE)
    cat("- Interactive HTML saved successfully.\n")
  }, error = function(e) {
    warning("Failed to save interactive HTML (pandoc may be missing): ", e$message, "\nStatic PNG is available.\n")
  })
 
  cat(sprintf("- %s plots saved: .png (bar and heatmap) and (if possible) .html files\n", dataset_name))
}

# Generate plots for each
create_and_save_plots(prok_with_data, "Prokaryote", "prok")
create_and_save_plots(arc_with_data, "Archaea", "arc")
create_and_save_plots(euk_with_data, "Eukaryote", "euk")

cat("\nScript execution completed successfully!\n")