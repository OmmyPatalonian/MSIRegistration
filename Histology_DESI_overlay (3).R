# Histology DESI Overlay Application
# A Shiny application for overlaying histology images with DESI mass spectrometry data

# Check and install required packages
required_packages <- c("shiny", "png", "jpeg", "ggplot2", "grid", "abind", "tools", "DT")
missing_packages <- c()

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
  
  # Load packages after installation
  for (pkg in missing_packages) {
    library(pkg, character.only = TRUE)
  }
}

# Check Cardinal (Bioconductor package)
if (!require("Cardinal", quietly = TRUE)) {
  cat("Installing Cardinal from Bioconductor...\n")
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    library(BiocManager)
  }
  BiocManager::install("Cardinal")
  library(Cardinal)
}

# Load all required libraries
library(shiny)
library(png)
library(jpeg)
library(ggplot2)
library(grid)
library(Cardinal)
library(abind)
library(tools)
library(DT)

# Source the plotting module
source("plot_Card_overlay_NEW (3).R")

# Function to extract MSI metadata and spatial information
extract_msi_metadata <- function(msi_data) {
  metadata_info <- list()
  
  tryCatch({
    # Basic information
    metadata_info$class <- class(msi_data)[1]
    metadata_info$dimensions <- paste(dim(msi_data), collapse = " x ")
    metadata_info$n_features <- nrow(msi_data)
    metadata_info$n_pixels <- ncol(msi_data)
    
    # Try to get pixel coordinates and calculate resolution
    if (exists("coord", where = msi_data) || "coord" %in% slotNames(msi_data)) {
      coords <- tryCatch(coord(msi_data), error = function(e) NULL)
      if (!is.null(coords)) {
        # Calculate pixel spacing in x and y directions
        unique_x <- sort(unique(coords$x))
        unique_y <- sort(unique(coords$y))
        
        if (length(unique_x) > 1) {
          x_spacing <- median(diff(unique_x), na.rm = TRUE)
          metadata_info$x_spacing <- x_spacing
        }
        if (length(unique_y) > 1) {
          y_spacing <- median(diff(unique_y), na.rm = TRUE)
          metadata_info$y_spacing <- y_spacing
        }
        
        metadata_info$x_range <- paste(range(coords$x), collapse = " to ")
        metadata_info$y_range <- paste(range(coords$y), collapse = " to ")
        metadata_info$spatial_dims <- paste(length(unique_x), "x", length(unique_y))
      }
    }
    
    # Try to get m/z information
    if (exists("mz", where = msi_data) || "mz" %in% slotNames(msi_data)) {
      mz_values <- tryCatch(mz(msi_data), error = function(e) NULL)
      if (!is.null(mz_values)) {
        metadata_info$mz_range <- paste(round(range(mz_values), 4), collapse = " to ")
        metadata_info$n_mz_values <- length(mz_values)
      }
    }
    
    # Try to extract pixel size from metadata if available
    if (exists("metadata", where = msi_data) || "metadata" %in% slotNames(msi_data)) {
      meta <- tryCatch(metadata(msi_data), error = function(e) NULL)
      if (!is.null(meta)) {
        # Look for common spatial resolution fields
        pixel_fields <- grep("pixel|spacing|resolution|size", names(meta), ignore.case = TRUE, value = TRUE)
        if (length(pixel_fields) > 0) {
          metadata_info$pixel_metadata <- meta[pixel_fields]
        }
      }
    }
    
    # Try to get run information
    if (exists("run", where = msi_data) || "run" %in% slotNames(msi_data)) {
      runs <- tryCatch(run(msi_data), error = function(e) NULL)
      if (!is.null(runs)) {
        metadata_info$n_runs <- length(unique(runs))
        metadata_info$runs <- paste(unique(runs), collapse = ", ")
      }
    }
    
  }, error = function(e) {
    metadata_info$error <- paste("Error extracting metadata:", e$message)
  })
  
  return(metadata_info)
}

# Function to estimate MSI pixel resolution from Cardinal data
estimate_msi_resolution <- function(msi_data) {
  resolution_estimate <- NULL
  
  tryCatch({
    # Try to get pixel coordinates
    coords <- coord(msi_data)
    
    if (!is.null(coords) && nrow(coords) > 1) {
      # Calculate median spacing between adjacent pixels
      unique_x <- sort(unique(coords$x))
      unique_y <- sort(unique(coords$y))
      
      x_spacing <- NA
      y_spacing <- NA
      
      if (length(unique_x) > 1) {
        x_spacing <- median(diff(unique_x), na.rm = TRUE)
      }
      if (length(unique_y) > 1) {
        y_spacing <- median(diff(unique_y), na.rm = TRUE)
      }
      
      # Use average spacing as pixel resolution estimate
      if (!is.na(x_spacing) && !is.na(y_spacing)) {
        resolution_estimate <- mean(c(x_spacing, y_spacing), na.rm = TRUE)
      } else if (!is.na(x_spacing)) {
        resolution_estimate <- x_spacing
      } else if (!is.na(y_spacing)) {
        resolution_estimate <- y_spacing
      }
      
      # For typical DESI data, convert to micrometers if needed
      # (This is an assumption - may need adjustment based on actual data units)
      if (!is.null(resolution_estimate) && resolution_estimate < 1) {
        # Likely in mm, convert to micrometers
        resolution_estimate <- resolution_estimate * 1000
      } else if (!is.null(resolution_estimate) && resolution_estimate > 1000) {
        # Likely already in micrometers, keep as is
        resolution_estimate <- resolution_estimate
      }
    }
  }, error = function(e) {
    cat("Could not estimate MSI resolution:", e$message, "\n")
  })
  
  return(resolution_estimate)
}

# Function to generate sample data for testing
generate_sample_data <- function() {
  cat("Generating sample data for testing...\n")
  
  # Create sample MSI data using Cardinal
  set.seed(123)
  sample_msi <- Cardinal::simulateImage(
    preset = 1,
    npeaks = 20,
    dim = c(20, 20),
    peakheight = c(3, 1),
    peakwidth = 500,
    sdnoise = 0.1,
    resolution = 1000,
    baseline = 1
  )
  
  # Save the sample MSI data
  saveRDS(sample_msi, "sample_msi_data.rds")
  
  # Create a simple sample histology image
  img_size <- 200
  sample_image <- array(0, dim = c(img_size, img_size, 3))
  
  # Create a gradient background with circular features
  for (i in 1:img_size) {
    for (j in 1:img_size) {
      sample_image[i, j, 1] <- 0.8 * (i / img_size)
      sample_image[i, j, 2] <- 0.6 * (j / img_size)
      sample_image[i, j, 3] <- 0.4
    }
  }
  
  # Add circular tissue-like structures
  centers <- list(c(50, 50), c(150, 50), c(100, 150), c(50, 150), c(150, 150))
  for (center in centers) {
    cx <- center[1]; cy <- center[2]; radius <- 25
    for (i in max(1, cx - radius):min(img_size, cx + radius)) {
      for (j in max(1, cy - radius):min(img_size, cy + radius)) {
        distance <- sqrt((i - cx)^2 + (j - cy)^2)
        if (distance <= radius) {
          factor <- 1 - (distance / radius) * 0.5
          sample_image[i, j, 1] <- sample_image[i, j, 1] * factor
          sample_image[i, j, 2] <- sample_image[i, j, 2] * factor
          sample_image[i, j, 3] <- sample_image[i, j, 3] * factor + 0.3 * (1 - factor)
        }
      }
    }
  }
  
  # Save the sample histology image
  png::writePNG(sample_image, "sample_histology.png")
  
  cat("Sample data generated successfully!\n")
  cat("Files created: sample_msi_data.rds, sample_histology.png\n")
}

# Uncomment the line below to generate sample data when sourcing this file
# generate_sample_data()

# Function to get MSI dimensions in pixels
get_msi_dimensions <- function(msi_data) {
  dims <- list(width = NULL, height = NULL)
  
  tryCatch({
    coords <- coord(msi_data)
    if (!is.null(coords) && nrow(coords) > 0) {
      unique_x <- unique(coords$x)
      unique_y <- unique(coords$y)
      dims$width <- length(unique_x)
      dims$height <- length(unique_y)
    }
  }, error = function(e) {
    cat("Error getting MSI dimensions:", e$message, "\n")
  })
  
  return(dims)
}

  # Set maximum file size for uploads (500MB)
options(shiny.maxRequestSize = 500*1024^2)

ui <- fluidPage(
  titlePanel("Histology DESI Overlay Application"),
  
  # Help text
  fluidRow(
    column(12,
      wellPanel(
        h4("Instructions:"),
        p("1. Upload a histology image (PNG/JPEG) and MSI dataset (.rds file)"),
        p("2. Use the transformation controls to align the histology image with MSI data"),
        p("3. Adjust transparency to see the overlay effect"),
        p("4. Save/load settings to preserve your transformations"),
        p(tags$strong("Note:"), "MSI data should be saved as .rds files containing Cardinal-compatible objects"),
        tags$hr(),
        p(tags$strong("Sample data available:"), 
          "sample_histology.png and sample_msi_data.rds in the current directory")
      )
    )
  ),
  
  # Use a sidebar layout with a panel for inputs and a main panel for the image output
  sidebarLayout(
    sidebarPanel(
      # File uploads
      fileInput("histology_upload", "Upload Histology Image (PNG/JPEG)", 
                accept = c("image/png", "image/jpeg", ".png", ".jpg", ".jpeg")),
      
      fileInput("msi_upload", "Upload MSI dataset (imzML + ibd or .rds)", 
                accept = c('.rds', '.RDS', '.imzML', '.ibd', 'application/octet-stream'),
                multiple = TRUE),
      
      # Image scaling
      sliderInput("scale_box", "Multiplier for images", 
                  min = 0.1, max = 10, value = 1, step = 0.1),
      
      # Save and load settings
      downloadButton("save_button", "Save Settings"),
      fileInput("load_button", "Load Settings", accept = c(".rds")),
      actionButton("restore_settings", "Restore Settings"),
      
      # MSI Metadata Display
      h4("MSI Data Information"),
      wellPanel(
        h5("MSI Dataset Metadata"),
        div(id = "msi_metadata", style = "font-size: 12px; background-color: #f8f9fa; padding: 10px; border-radius: 5px;",
            verbatimTextOutput("msi_metadata_display")
        ),
        br(),
        fluidRow(
          column(6, textOutput("msi_pixel_info")),
          column(6, textOutput("msi_spatial_info"))
        )
      ),
      
      # Resolution and scaling controls
      h4("Resolution & Scaling"),
      wellPanel(
        h5("Resolution Settings"),
        fluidRow(
          column(6, numericInput("histology_pixel_width", "Histology Width (pixels)", 
                               value = 1024, min = 100, max = 10000, step = 1)),
          column(6, numericInput("histology_pixel_height", "Histology Height (pixels)", 
                               value = 1024, min = 100, max = 10000, step = 1))
        ),
        fluidRow(
          column(6, numericInput("msi_pixel_width", "MSI Width (pixels)", 
                              value = 20, min = 1, max = 1000, step = 1)),
          column(6, numericInput("msi_pixel_height", "MSI Height (pixels)", 
                              value = 20, min = 1, max = 1000, step = 1))
        ),
        fluidRow(
          column(6, numericInput("histology_microns_per_pixel", "Histology μm/pixel", 
                               value = 0.5, min = 0.01, max = 100, step = 0.01)),
          column(6, numericInput("msi_microns_per_pixel", "MSI μm/pixel", 
                               value = 20, min = 0.1, max = 1000, step = 0.1))
        ),
        fluidRow(
          column(6, 
            checkboxInput("auto_scale_resolution", "Auto-scale based on resolution", value = TRUE)
          ),
          column(6,
            selectInput("scale_target", "Scale to target", 
                        choices = c("MSI Grid" = "msi", "Histology Grid" = "histology", "Best Fit" = "best_fit"),
                        selected = "best_fit")
          )
        ),
        fluidRow(
          column(12,
            div(id = "resolution_info", style = "font-size: 12px; color: #666;",
                textOutput("resolution_display")
            )
          )
        ),
        fluidRow(
          column(12,
            div(id = "physical_size_info", style = "font-size: 12px; color: #666;",
                textOutput("physical_size_display")
            )
          )
        )
      ),
      
      # Image transformation controls
      h4("Manual Transformations"),
      
      # Scaling controls
      fluidRow(
        column(6, sliderInput("scalex", "Scale X", 
                            min = 0.1, max = 5, value = 1, step = 0.001)),
        column(6, sliderInput("scaley", "Scale Y", 
                            min = 0.1, max = 5, value = 1, step = 0.001))
      ),
      
      # Rotation controls with buttons
      h5("Rotation"),
      fluidRow(
        column(3, actionButton("rotate_left", "↺ -5°", class = "btn-primary btn-sm")),
        column(3, actionButton("rotate_right", "↻ +5°", class = "btn-primary btn-sm")),
        column(6, numericInput("rotate", "Degrees:", value = 0, min = -180, max = 180, step = 0.1))
      ),
      fluidRow(
        column(6, actionButton("rotate_left_fine", "↺ -1°", class = "btn-outline-primary btn-sm")),
        column(6, actionButton("rotate_right_fine", "↻ +1°", class = "btn-outline-primary btn-sm"))
      ),
      
      # Translation controls with arrow buttons
      h5("Translation"),
      fluidRow(
        column(12, 
          div(style = "text-align: center; margin: 10px 0;",
            actionButton("move_up", "▲", class = "btn-success btn-sm", style = "display: block; margin: 0 auto; width: 40px;")
          )
        )
      ),
      fluidRow(
        column(4, actionButton("move_left", "◀", class = "btn-success btn-sm", style = "width: 40px;")),
        column(4, div(style = "text-align: center;", actionButton("center_image", "⌂", class = "btn-warning btn-sm", style = "width: 40px;"))),
        column(4, actionButton("move_right", "▶", class = "btn-success btn-sm", style = "width: 40px;"))
      ),
      fluidRow(
        column(12, 
          div(style = "text-align: center; margin: 10px 0;",
            actionButton("move_down", "▼", class = "btn-success btn-sm", style = "display: block; margin: 0 auto; width: 40px;")
          )
        )
      ),
      
      # Fine control sliders (for precise adjustments)
      h5("Fine Controls"),
      sliderInput("translate_x", "Fine X Position", 
                  min = -100, max = 100, value = 0, step = 0.1),
      sliderInput("translate_y", "Fine Y Position", 
                  min = -100, max = 100, value = 0, step = 0.1),
      
      # Display options
      checkboxInput("spatial", "Spatial plot only", value = TRUE),
      checkboxInput("debug", "Debug plot?", value = FALSE),
      
      # Advanced rendering options
      h4("Advanced Options"),
      wellPanel(
        selectInput("render_mode", "MSI Rendering Mode", 
                    choices = c("Standard" = "standard", "Clean Pixels" = "clean"),
                    selected = "standard"),
        helpText("Clean mode may improve overlay quality but take longer to render")
      )
    ),
    
    mainPanel(
      # Output the overlay plot
      plot_card_UI("hist_plot_card"),
      sliderInput("alpha", "Adjust Image Transparency", 
                  min = 0, max = 1, value = 0.5)
    )
  )
)
server <- function(input, output, session) {
  
  # Create a reactiveValues object for all inputs
  allInputs <- reactiveValues(
    # Add a default max scale factor to prevent memory issues with extreme resolution differences
    max_scale_factor = 5.0,
    # Flag to check if using sample data
    is_sample_data = FALSE,
    # Default values for resolution settings
    histology_microns_per_pixel = 0.5,
    msi_microns_per_pixel = 20,
    scale_target = "best_fit"
  )
  
  # Reactive value to store MSI data for metadata extraction
  msi_data_reactive <- reactiveVal(NULL)
  
  # Reactive values for settings
  settings <- reactiveValues(
    scalex = 1,
    scaley = 1,
    rotate = 0,
    translate_x = 0,
    translate_y = 0,
    alpha = 0.5
  )
  
  # Update allInputs whenever any input changes
  observe({
    inputList <- reactiveValuesToList(input)
    for (name in names(inputList)) {
      allInputs[[name]] <- inputList[[name]]
    }
    
    # Ensure render_mode is available
    if (is.null(allInputs$render_mode)) {
      allInputs$render_mode <- "standard"
    }
  })
  
  # Make sure dimension changes get reflected immediately and trigger recalculation
  update_dimension_scaling <- function() {
    # Ensure all necessary values are available
    if (is.null(input$histology_pixel_width) || is.null(input$msi_pixel_width) ||
        is.null(input$histology_pixel_height) || is.null(input$msi_pixel_height) ||
        is.na(input$histology_pixel_width) || is.na(input$msi_pixel_width) || 
        is.na(input$histology_pixel_height) || is.na(input$msi_pixel_height) || 
        input$histology_pixel_width <= 0 || input$msi_pixel_width <= 0 ||
        input$histology_pixel_height <= 0 || input$msi_pixel_height <= 0) {
      return()
    }
    
    # Update allInputs
    allInputs$histology_pixel_width <- input$histology_pixel_width
    allInputs$histology_pixel_height <- input$histology_pixel_height
    allInputs$msi_pixel_width <- input$msi_pixel_width
    allInputs$msi_pixel_height <- input$msi_pixel_height
    
    # Force a plot refresh
    allInputs$dimension_change_trigger <- runif(1)
    
    # Log the change
    cat("\nDimension change detected:\n")
    cat(" - Histology dimensions:", input$histology_pixel_width, "×", input$histology_pixel_height, "pixels\n")
    cat(" - MSI dimensions:", input$msi_pixel_width, "×", input$msi_pixel_height, "pixels\n")
  }
  
  # Observe changes to dimensions
  observeEvent(input$histology_pixel_width, {
    allInputs$histology_pixel_width <- input$histology_pixel_width
    update_dimension_scaling()
  })
  
  observeEvent(input$histology_pixel_height, {
    allInputs$histology_pixel_height <- input$histology_pixel_height
    update_dimension_scaling()
  })
  
  observeEvent(input$msi_pixel_width, {
    allInputs$msi_pixel_width <- input$msi_pixel_width
    update_dimension_scaling()
  })
  
  observeEvent(input$msi_pixel_height, {
    allInputs$msi_pixel_height <- input$msi_pixel_height
    update_dimension_scaling()
  })
  
  # Make sure resolution values are reflected immediately
  observeEvent(input$histology_microns_per_pixel, {
    allInputs$histology_microns_per_pixel <- input$histology_microns_per_pixel
    # Force plot refresh
    allInputs$resolution_change_trigger <- runif(1)
  }, priority = 10)
  
  observeEvent(input$msi_microns_per_pixel, {
    allInputs$msi_microns_per_pixel <- input$msi_microns_per_pixel
    # Force plot refresh
    allInputs$resolution_change_trigger <- runif(1)
  }, priority = 10)
  
  # Make sure scale target changes are reflected immediately
  observeEvent(input$scale_target, {
    allInputs$scale_target <- input$scale_target
    # Force plot refresh
    allInputs$resolution_change_trigger <- runif(1)
  }, priority = 10)
  
  # Make sure render mode changes are reflected immediately
  observeEvent(input$render_mode, {
    allInputs$render_mode <- input$render_mode
    # Force refresh of plot when render mode changes
    allInputs$render_mode_change <- runif(1)
  })
  
  # Main observer for MSI data processing
  observe({
    req(input$histology_upload)
    req(input$msi_upload)
    
    cat("Files uploaded:\n")
    cat("- Histology:", input$histology_upload$name, "\n")
    cat("- MSI:", paste(input$msi_upload$name, collapse = ", "), "\n")
    
    # Determine MSI file type and load accordingly
    msi_files <- input$msi_upload
    msi_paths <- msi_files$datapath
    msi_names <- msi_files$name
    
    # If .rds, load with readRDS
    if (any(grepl("\\.rds$|\\.RDS$", msi_names, ignore.case = TRUE))) {
      rds_idx <- which(grepl("\\.rds$|\\.RDS$", msi_names, ignore.case = TRUE))[1]
      tryCatch({
        cat("Attempting to load MSI data from:", msi_paths[rds_idx], "\n")
        dat_in <- readRDS(msi_paths[rds_idx])
        cat("MSI data loaded successfully. Class:", class(dat_in), "\n")
        
        # Store MSI data for metadata extraction
        msi_data_reactive(dat_in)
        
        # Extract MSI width and height and update UI
        tryCatch({
          coords <- coord(dat_in)
          if (!is.null(coords)) {
            msi_width <- length(unique(coords$x))
            msi_height <- length(unique(coords$y))
            cat("MSI dimensions:", msi_width, "×", msi_height, "pixels\n")
            
            # Update UI with actual MSI dimensions
            updateNumericInput(session, "msi_pixel_width", value = msi_width)
            updateNumericInput(session, "msi_pixel_height", value = msi_height)
          }
        }, error = function(e) {
          cat("Error extracting MSI dimensions:", e$message, "\n")
        })
        
        # Try to estimate and update MSI resolution
        estimated_resolution <- estimate_msi_resolution(dat_in)
        if (!is.null(estimated_resolution) && !is.na(estimated_resolution)) {
          updateNumericInput(session, "msi_microns_per_pixel", value = round(estimated_resolution, 2))
          showNotification(paste("Auto-detected MSI resolution:", round(estimated_resolution, 2), "μm/pixel"), 
                          type = "message")
        }
        
        # Set sample data flag if needed
        if (any(msi_names == "sample_msi_data.rds")) {
          allInputs$is_sample_data <- TRUE
          cat("Sample MSI data detected\n")
        }
        
        plot_card_server("hist_plot_card", dat_in, 
                         spatialOnly = input$spatial, 
                         allInputs = allInputs)
      }, error = function(e) {
        cat("Error loading MSI data:", e$message, "\n")
        showNotification(paste("Error loading MSI data:", e$message), 
                         type = "error")
      })
    } else if (any(grepl("\\.imzML$", msi_names, ignore.case = TRUE))) {
      # Find imzML and ibd files
      imzml_idx <- which(grepl("\\.imzML$", msi_names, ignore.case = TRUE))[1]
      ibd_idx <- which(grepl("\\.ibd$", msi_names, ignore.case = TRUE))[1]
      
      if (is.na(ibd_idx)) {
        showNotification("Please upload both .imzML and .ibd files", type = "error")
        return()
      }
      
      # Create temporary directory with proper file names
      temp_dir <- tempdir()
      
      # Get original file names (without extension)
      imzml_basename <- tools::file_path_sans_ext(msi_names[imzml_idx])
      
      # Copy files to temp directory with proper names
      temp_imzml_path <- file.path(temp_dir, paste0(imzml_basename, ".imzML"))
      temp_ibd_path <- file.path(temp_dir, paste0(imzml_basename, ".ibd"))
      
      tryCatch({
        # Copy files to temporary location with matching names
        file.copy(msi_paths[imzml_idx], temp_imzml_path, overwrite = TRUE)
        file.copy(msi_paths[ibd_idx], temp_ibd_path, overwrite = TRUE)
        
        cat("Attempting to load MSI data from imzML:", temp_imzml_path, "\n")
        cat("Corresponding ibd file:", temp_ibd_path, "\n")
        
        # Cardinal should now find the .ibd file
        dat_in <- Cardinal::readMSIData(temp_imzml_path)
        cat("MSI data loaded successfully from imzML. Class:", class(dat_in), "\n")
        cat("Data dimensions:", dim(dat_in), "\n")
        cat("Number of features:", nrow(dat_in), "\n")
        cat("Number of pixels:", ncol(dat_in), "\n")
        
        # Get MSI image dimensions directly from coordinates
        tryCatch({
          lengthx <- length(unique(coord(dat_in)$x))
          lengthy <- length(unique(coord(dat_in)$y))
          cat("MSI dimensions (from coords):", lengthx, "×", lengthy, "pixels\n")
          
          # Update UI with actual MSI dimensions
          updateNumericInput(session, "msi_pixel_width", value = lengthx)
          updateNumericInput(session, "msi_pixel_height", value = lengthy)
        }, error = function(e) {
          cat("Could not determine MSI dimensions from coordinates:", e$message, "\n")
        })
        
        # Store MSI data for metadata extraction
        msi_data_reactive(dat_in)
        
        # Try to estimate and update MSI resolution
        estimated_resolution <- estimate_msi_resolution(dat_in)
        if (!is.null(estimated_resolution) && !is.na(estimated_resolution)) {
          updateNumericInput(session, "msi_microns_per_pixel", value = round(estimated_resolution, 2))
          showNotification(paste("Auto-detected MSI resolution:", round(estimated_resolution, 2), "μm/pixel"), 
                          type = "message")
        }
        
        # Set sample data flag if needed
        if (any(msi_names == "sample_msi_data.rds")) {
          allInputs$is_sample_data <- TRUE
          cat("Sample MSI data detected\n")
        }
        
        # Try to call the plotting server with error handling
        tryCatch({
          plot_card_server("hist_plot_card", dat_in, 
                           spatialOnly = input$spatial, 
                           allInputs = allInputs)
          cat("plot_card_server called successfully\n")
        }, error = function(plot_error) {
          cat("Error in plot_card_server:", plot_error$message, "\n")
          showNotification(paste("Error in plotting module:", plot_error$message), type = "error")
        })
      }, error = function(e) {
        cat("Error loading MSI data (imzML):", e$message, "\n")
        showNotification(paste("Error loading MSI data (imzML):", e$message), 
                         type = "error")
      })
    } else {
      showNotification("Please upload a valid MSI dataset (.rds or .imzML + .ibd)", type = "error")
    }
  })
  
  # Debug observers for file uploads
  observeEvent(input$histology_upload, {
    if (!is.null(input$histology_upload)) {
      cat("Histology file uploaded:", input$histology_upload$name, 
          "Size:", input$histology_upload$size, "bytes\n")
      showNotification(paste("Histology image uploaded:", input$histology_upload$name), 
                       type = "message")
      
      # Check if this is the sample data
      if (input$histology_upload$name == "sample_histology.png") {
        allInputs$is_sample_data <- TRUE
        cat("Sample histology data detected\n")
      } else {
        allInputs$is_sample_data <- FALSE
      }
    }
  })
  
  observeEvent(input$msi_upload, {
    if (!is.null(input$msi_upload)) {
      cat("MSI file uploaded:", paste(input$msi_upload$name, collapse = ", "), 
          "Size:", paste(input$msi_upload$size, collapse = ", "), "bytes\n")
      
      # Create a proper single character string for notification
      file_names <- paste(input$msi_upload$name, collapse = ", ")
      notification_msg <- paste("MSI dataset uploaded:", file_names)
      showNotification(notification_msg, type = "message")
    }
  })

  # MSI Metadata outputs
  output$msi_metadata_display <- renderText({
    msi_data <- msi_data_reactive()
    if (is.null(msi_data)) {
      return("No MSI data loaded yet. Please upload an MSI dataset.")
    }
    
    metadata <- extract_msi_metadata(msi_data)
    
    # Format metadata for display
    metadata_text <- c()
    
    if (!is.null(metadata$class)) {
      metadata_text <- c(metadata_text, paste("Data Type:", metadata$class))
    }
    if (!is.null(metadata$dimensions)) {
      metadata_text <- c(metadata_text, paste("Dimensions:", metadata$dimensions))
    }
    if (!is.null(metadata$n_features)) {
      metadata_text <- c(metadata_text, paste("Number of m/z features:", metadata$n_features))
    }
    if (!is.null(metadata$n_pixels)) {
      metadata_text <- c(metadata_text, paste("Number of pixels:", metadata$n_pixels))
    }
    if (!is.null(metadata$mz_range)) {
      metadata_text <- c(metadata_text, paste("m/z range:", metadata$mz_range))
    }
    if (!is.null(metadata$spatial_dims)) {
      metadata_text <- c(metadata_text, paste("Spatial dimensions:", metadata$spatial_dims))
    }
    if (!is.null(metadata$x_range)) {
      metadata_text <- c(metadata_text, paste("X coordinate range:", metadata$x_range))
    }
    if (!is.null(metadata$y_range)) {
      metadata_text <- c(metadata_text, paste("Y coordinate range:", metadata$y_range))
    }
    if (!is.null(metadata$x_spacing)) {
      metadata_text <- c(metadata_text, paste("X pixel spacing:", round(metadata$x_spacing, 4)))
    }
    if (!is.null(metadata$y_spacing)) {
      metadata_text <- c(metadata_text, paste("Y pixel spacing:", round(metadata$y_spacing, 4)))
    }
    if (!is.null(metadata$n_runs)) {
      metadata_text <- c(metadata_text, paste("Number of runs:", metadata$n_runs))
    }
    if (!is.null(metadata$pixel_metadata)) {
      metadata_text <- c(metadata_text, "--- Pixel Metadata ---")
      for (i in seq_along(metadata$pixel_metadata)) {
        metadata_text <- c(metadata_text, paste(names(metadata$pixel_metadata)[i], ":", metadata$pixel_metadata[[i]]))
      }
    }
    if (!is.null(metadata$error)) {
      metadata_text <- c(metadata_text, paste("ERROR:", metadata$error))
    }
    
    if (length(metadata_text) == 0) {
      return("Metadata extraction completed but no information found.")
    }
    
    return(paste(metadata_text, collapse = "\n"))
  })
  
  output$msi_pixel_info <- renderText({
    msi_data <- msi_data_reactive()
    if (is.null(msi_data)) {
      return("MSI Pixel Info: Not available")
    }
    
    metadata <- extract_msi_metadata(msi_data)
    info_parts <- c()
    
    if (!is.null(metadata$n_pixels)) {
      info_parts <- c(info_parts, paste("Pixels:", metadata$n_pixels))
    }
    if (!is.null(metadata$spatial_dims)) {
      info_parts <- c(info_parts, paste("Grid:", metadata$spatial_dims))
    }
    
    if (length(info_parts) > 0) {
      return(paste("MSI Info:", paste(info_parts, collapse = ", ")))
    } else {
      return("MSI Pixel Info: Not available")
    }
  })
  
  output$msi_spatial_info <- renderText({
    msi_data <- msi_data_reactive()
    if (is.null(msi_data)) {
      return("Spatial Resolution: Not available")
    }
    
    estimated_res <- estimate_msi_resolution(msi_data)
    if (!is.null(estimated_res) && !is.na(estimated_res)) {
      return(paste("Estimated Resolution:", round(estimated_res, 2), "μm/pixel"))
    } else {
      return("Spatial Resolution: Could not estimate")
    }
  })

  # Utility functions for micron/pixel conversion
  phys_size_um <- function(dims_px, scale_um_per_px) {
    dims_px * scale_um_per_px
  }
  
  resampled_dims <- function(dims_px, source_um_per_px, target_um_per_px) {
    round(dims_px * (source_um_per_px / target_um_per_px))
  }
  
  px_to_um <- function(xy_px, scale_um_per_px) {
    xy_px * scale_um_per_px
  }
  
  um_to_px <- function(xy_um, scale_um_per_px) {
    xy_um / scale_um_per_px
  }
  
  # Resolution calculation and display
  output$resolution_display <- renderText({
    req(input$histology_microns_per_pixel, input$msi_microns_per_pixel)
    
    resolution_ratio <- input$msi_microns_per_pixel / input$histology_microns_per_pixel
    histology_field_width <- input$histology_pixel_width * input$histology_microns_per_pixel
    histology_field_height <- input$histology_pixel_height * input$histology_microns_per_pixel
    
    paste0("Resolution ratio (MSI:Histology): ", round(resolution_ratio, 2), 
           " | Histology field size: ", round(histology_field_width), "×", 
           round(histology_field_height), " μm")
  })
  
  # Physical size display
  output$physical_size_display <- renderText({
    req(input$histology_microns_per_pixel, input$msi_microns_per_pixel)
    req(input$histology_pixel_width, input$histology_pixel_height)
    
    # Get MSI dimensions
    msi_data <- msi_data_reactive()
    msi_width <- NULL
    msi_height <- NULL
    
    if (!is.null(msi_data)) {
      tryCatch({
        coords <- coord(msi_data)
        if (!is.null(coords)) {
          msi_width <- length(unique(coords$x))
          msi_height <- length(unique(coords$y))
        }
      }, error = function(e) {
        cat("Error getting MSI dimensions:", e$message, "\n")
      })
    }
    
    # Calculate physical sizes
    if (!is.null(msi_width) && !is.null(msi_height)) {
      msi_phys_width <- msi_width * input$msi_microns_per_pixel
      msi_phys_height <- msi_height * input$msi_microns_per_pixel
      
      histo_phys_width <- input$histology_pixel_width * input$histology_microns_per_pixel
      histo_phys_height <- input$histology_pixel_height * input$histology_microns_per_pixel
      
      if (input$scale_target == "msi") {
        histo_on_msi <- resampled_dims(
          c(input$histology_pixel_width, input$histology_pixel_height), 
          input$histology_microns_per_pixel, 
          input$msi_microns_per_pixel)
        
        return(paste0("Physical sizes - MSI: ", round(msi_phys_width, 1), "×", round(msi_phys_height, 1), " μm | ",
                     "Histology: ", round(histo_phys_width, 1), "×", round(histo_phys_height, 1), " μm | ",
                     "Histology on MSI grid: ", histo_on_msi[1], "×", histo_on_msi[2], " pixels"))
      } else if (input$scale_target == "histology") {
        msi_on_histo <- resampled_dims(
          c(msi_width, msi_height), 
          input$msi_microns_per_pixel, 
          input$histology_microns_per_pixel)
        
        return(paste0("Physical sizes - MSI: ", round(msi_phys_width, 1), "×", round(msi_phys_height, 1), " μm | ",
                     "Histology: ", round(histo_phys_width, 1), "×", round(histo_phys_height, 1), " μm | ",
                     "MSI on histology grid: ", msi_on_histo[1], "×", msi_on_histo[2], " pixels"))
      } else {
        return(paste0("Physical sizes - MSI: ", round(msi_phys_width, 1), "×", round(msi_phys_height, 1), " μm | ",
                     "Histology: ", round(histo_phys_width, 1), "×", round(histo_phys_height, 1), " μm"))
      }
    } else {
      histo_phys_width <- input$histology_pixel_width * input$histology_microns_per_pixel
      histo_phys_height <- input$histology_pixel_height * input$histology_microns_per_pixel
      
      return(paste0("Physical size - Histology: ", round(histo_phys_width, 1), "×", round(histo_phys_height, 1), " μm | ",
                   "MSI dimensions not available"))
    }
  })
  
  # Auto-scale based on resolution and dimensions when enabled
  observeEvent(list(input$auto_scale_resolution, 
                   input$histology_microns_per_pixel, input$msi_microns_per_pixel,
                   input$histology_pixel_width, input$histology_pixel_height,
                   input$msi_pixel_width, input$msi_pixel_height,
                   input$scale_target), {
    
    req(input$histology_microns_per_pixel, input$msi_microns_per_pixel)
    req(input$scale_target)
    
    # Update allInputs immediately with current resolution and dimension values
    allInputs$histology_microns_per_pixel <- input$histology_microns_per_pixel
    allInputs$msi_microns_per_pixel <- input$msi_microns_per_pixel
    allInputs$histology_pixel_width <- input$histology_pixel_width
    allInputs$histology_pixel_height <- input$histology_pixel_height
    allInputs$msi_pixel_width <- input$msi_pixel_width
    allInputs$msi_pixel_height <- input$msi_pixel_height
    allInputs$scale_target <- input$scale_target
    
    # Force refresh plot when parameters change
    allInputs$resolution_change_trigger <- runif(1)
    
    # Log current values for debugging
    cat("\nAuto-scale observer triggered with values:\n")
    cat(" - Histology dimensions:", input$histology_pixel_width, "×", input$histology_pixel_height, "pixels\n")
    cat(" - MSI dimensions:", input$msi_pixel_width, "×", input$msi_pixel_height, "pixels\n")
    cat(" - Histology resolution:", input$histology_microns_per_pixel, "μm/pixel\n")
    cat(" - MSI resolution:", input$msi_microns_per_pixel, "μm/pixel\n")
    cat(" - Scale target:", input$scale_target, "\n")
    cat(" - Auto-scale enabled:", input$auto_scale_resolution, "\n")
    
    # Calculate pixel dimension scale factors
    pixel_scale_x <- 1
    pixel_scale_y <- 1
    
    # Only apply pixel dimension scaling if all values are valid
    if (!is.null(input$histology_pixel_width) && !is.null(input$msi_pixel_width) &&
        !is.null(input$histology_pixel_height) && !is.null(input$msi_pixel_height) &&
        !is.na(input$histology_pixel_width) && !is.na(input$msi_pixel_width) && 
        !is.na(input$histology_pixel_height) && !is.na(input$msi_pixel_height) && 
        input$histology_pixel_width > 0 && input$msi_pixel_width > 0 &&
        input$histology_pixel_height > 0 && input$msi_pixel_height > 0) {
        
        # Calculate pixel dimension ratios
        pixel_scale_x <- input$msi_pixel_width / input$histology_pixel_width  # Inverted ratio from render code
        pixel_scale_y <- input$msi_pixel_height / input$histology_pixel_height # Inverted ratio from render code
        
        cat("Pixel dimension scale factors (MSI:Histology) - X:", pixel_scale_x, "Y:", pixel_scale_y, "\n")
    } else {
        cat("Pixel dimension scaling not applied (missing or invalid dimensions)\n")
    }
    
    # Update scales when auto-scale is enabled
    if (input$auto_scale_resolution) {
      # Determine target microns per pixel based on scale target
      target_um_per_px <- switch(input$scale_target,
        "msi"       = input$msi_microns_per_pixel,
        "histology" = input$histology_microns_per_pixel,
        "best_fit"  = sqrt(input$msi_microns_per_pixel * input$histology_microns_per_pixel)
      )
      
      # Calculate resolution-based scale factors
      scale_histology <- input$histology_microns_per_pixel / target_um_per_px
      scale_msi <- input$msi_microns_per_pixel / target_um_per_px
      
      # For diagnostic purposes
      cat("Resolution-based scaling:\n")
      cat(" - Target µm/pixel:", target_um_per_px, "\n")
      cat(" - Histology scale factor:", scale_histology, "\n")
      cat(" - MSI scale factor:", scale_msi, "\n")
      
      # Combine resolution scaling with pixel dimension scaling
      # We're scaling the histology image in this application
      new_scale_x <- scale_histology * pixel_scale_x
      new_scale_y <- scale_histology * pixel_scale_y
      
      cat("Combined scale factors - X:", new_scale_x, "Y:", new_scale_y, "\n")
      
      # Adjust for aspect ratio preservation
      if (abs(new_scale_x - new_scale_y) < 0.1) {
        # If scales are similar, use average to maintain aspect ratio
        avg_scale <- (new_scale_x + new_scale_y) / 2
        new_scale_x <- avg_scale
        new_scale_y <- avg_scale
        cat("Preserving aspect ratio with average scale:", avg_scale, "\n")
      }
      
      # Update both scales with combined factors
      updateSliderInput(session, "scalex", value = new_scale_x)
      updateSliderInput(session, "scaley", value = new_scale_y)
    } else {
      # If auto-scale is disabled but we have valid pixel dimensions,
      # still apply basic pixel ratio scaling
      if (pixel_scale_x != 1 || pixel_scale_y != 1) {
        cat("Auto-scale disabled but applying basic pixel dimension scaling\n")
        updateSliderInput(session, "scalex", value = pixel_scale_x)
        updateSliderInput(session, "scaley", value = pixel_scale_y)
      }
    }
  }, ignoreInit = FALSE)
  
  # Update histology dimensions when image is uploaded
  observe({
    req(input$histology_upload)
    
    tryCatch({
      # Try to read image dimensions
      img_path <- input$histology_upload$datapath
      file_type <- tolower(tools::file_ext(img_path))
      
      if (file_type %in% c("png")) {
        img <- png::readPNG(img_path)
        img_dims <- dim(img)
        # Immediately update allInputs as well to ensure these values are used
        allInputs$histology_pixel_width <- img_dims[2]
        allInputs$histology_pixel_height <- img_dims[1]
        updateNumericInput(session, "histology_pixel_width", value = img_dims[2])
        updateNumericInput(session, "histology_pixel_height", value = img_dims[1])
        
        showNotification(paste("Auto-detected histology dimensions:", 
                              img_dims[2], "×", img_dims[1], "pixels"), 
                        type = "message")
        
      } else if (file_type %in% c("jpg", "jpeg")) {
        img <- jpeg::readJPEG(img_path)
        img_dims <- dim(img)
        # Immediately update allInputs as well
        allInputs$histology_pixel_width <- img_dims[2]
        allInputs$histology_pixel_height <- img_dims[1]
        updateNumericInput(session, "histology_pixel_width", value = img_dims[2])
        updateNumericInput(session, "histology_pixel_height", value = img_dims[1])
        
        showNotification(paste("Auto-detected histology dimensions:", 
                              img_dims[2], "×", img_dims[1], "pixels"), 
                        type = "message")
      }
    }, error = function(e) {
      cat("Could not auto-detect image dimensions:", e$message, "\n")
    })
  })

  # Save settings download handler
  output$save_button <- downloadHandler(
    filename = function() {
      paste("settings-", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      settings_to_save <- list(
        scalex = input$scalex,
        scaley = input$scaley,
        rotate = input$rotate,
        translate_x = input$translate_x,
        translate_y = input$translate_y,
        alpha = input$alpha,
        histology_pixel_width = input$histology_pixel_width,
        histology_pixel_height = input$histology_pixel_height,
        histology_microns_per_pixel = input$histology_microns_per_pixel,
        msi_microns_per_pixel = input$msi_microns_per_pixel,
        auto_scale_resolution = input$auto_scale_resolution,
        scale_target = input$scale_target
      )
      saveRDS(settings_to_save, file)
    }
  )
  
  # Load settings
  settings_to_load <- reactive({
    req(input$load_button)
    tryCatch({
      readRDS(input$load_button$datapath)
    }, error = function(e) {
      showNotification("Error loading settings file", type = "error")
      return(NULL)
    })
  })
  
  # Restore settings button
  observeEvent(input$restore_settings, {
    loaded_settings <- settings_to_load()
    req(loaded_settings)
    
    # Update transformation inputs with loaded values
    updateSliderInput(session, "scalex", value = loaded_settings$scalex)
    updateSliderInput(session, "scaley", value = loaded_settings$scaley)
    updateNumericInput(session, "rotate", value = loaded_settings$rotate)
    updateSliderInput(session, "translate_x", value = loaded_settings$translate_x)
    updateSliderInput(session, "translate_y", value = loaded_settings$translate_y)
    updateSliderInput(session, "alpha", value = loaded_settings$alpha)
    
    # Update resolution inputs if they exist in saved settings
    if (!is.null(loaded_settings$histology_pixel_width)) {
      updateNumericInput(session, "histology_pixel_width", value = loaded_settings$histology_pixel_width)
    }
    if (!is.null(loaded_settings$histology_pixel_height)) {
      updateNumericInput(session, "histology_pixel_height", value = loaded_settings$histology_pixel_height)
    }
    if (!is.null(loaded_settings$histology_microns_per_pixel)) {
      updateNumericInput(session, "histology_microns_per_pixel", value = loaded_settings$histology_microns_per_pixel)
    }
    if (!is.null(loaded_settings$msi_microns_per_pixel)) {
      updateNumericInput(session, "msi_microns_per_pixel", value = loaded_settings$msi_microns_per_pixel)
    }
    if (!is.null(loaded_settings$auto_scale_resolution)) {
      updateCheckboxInput(session, "auto_scale_resolution", value = loaded_settings$auto_scale_resolution)
    }
    if (!is.null(loaded_settings$scale_target)) {
      updateSelectInput(session, "scale_target", selected = loaded_settings$scale_target)
    }
    
    showNotification("Settings restored successfully!", type = "message")
  })
  

  
  # Movement step sizes
  move_step_large <- 5  # pixels for arrow buttons
  move_step_small <- 1  # pixels for small adjustments
  
  # Translation button event handlers
  observeEvent(input$move_up, {
    new_value <- input$translate_y + move_step_large
    new_value <- max(-100, min(100, new_value))  # Keep within bounds
    updateSliderInput(session, "translate_y", value = new_value)
  })
  
  observeEvent(input$move_down, {
    new_value <- input$translate_y - move_step_large
    new_value <- max(-100, min(100, new_value))  # Keep within bounds
    updateSliderInput(session, "translate_y", value = new_value)
  })
  
  observeEvent(input$move_left, {
    new_value <- input$translate_x - move_step_large
    new_value <- max(-100, min(100, new_value))  # Keep within bounds
    updateSliderInput(session, "translate_x", value = new_value)
  })
  
  observeEvent(input$move_right, {
    new_value <- input$translate_x + move_step_large
    new_value <- max(-100, min(100, new_value))  # Keep within bounds
    updateSliderInput(session, "translate_x", value = new_value)
  })
  
  # Center button - reset translation to zero
  observeEvent(input$center_image, {
    updateSliderInput(session, "translate_x", value = 0)
    updateSliderInput(session, "translate_y", value = 0)
    showNotification("Image centered", type = "message")
  })
  
  # Rotation button event handlers
  observeEvent(input$rotate_left, {
    new_value <- input$rotate - 5
    new_value <- max(-180, min(180, new_value))  # Keep within bounds
    updateNumericInput(session, "rotate", value = new_value)
  })
  
  observeEvent(input$rotate_right, {
    new_value <- input$rotate + 5
    new_value <- max(-180, min(180, new_value))  # Keep within bounds
    updateNumericInput(session, "rotate", value = new_value)
  })
  
  observeEvent(input$rotate_left_fine, {
    new_value <- input$rotate - 1
    new_value <- max(-180, min(180, new_value))  # Keep within bounds
    updateNumericInput(session, "rotate", value = new_value)
  })
  
  observeEvent(input$rotate_right_fine, {
    new_value <- input$rotate + 1
    new_value <- max(-180, min(180, new_value))  # Keep within bounds
    updateNumericInput(session, "rotate", value = new_value)
  })
}

# Run the application
cat("Starting Histology DESI Overlay Application...\n")
cat("Version: 1.0\n")
cat("For sample data, uncomment generate_sample_data() and re-source this file.\n")
cat("Loading application...\n")

shinyApp(ui, server)
