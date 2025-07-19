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
      
      # Image transformation controls
      h4("Image Transformations"),
      
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
      checkboxInput("debug", "Debug plot?", value = FALSE)
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
  allInputs <- reactiveValues()
  
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
        alpha = input$alpha
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
    
    # Update slider inputs with loaded values
    updateSliderInput(session, "scalex", value = loaded_settings$scalex)
    updateSliderInput(session, "scaley", value = loaded_settings$scaley)
    updateSliderInput(session, "rotate", value = loaded_settings$rotate)
    updateSliderInput(session, "translate_x", value = loaded_settings$translate_x)
    updateSliderInput(session, "translate_y", value = loaded_settings$translate_y)
    updateSliderInput(session, "alpha", value = loaded_settings$alpha)
    
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
