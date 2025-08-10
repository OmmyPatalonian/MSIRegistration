# Plot Card Overlay Module for Histology DESI Overlay Application
# 
# This module provides the UI and server functions for displaying and manipulating
# mass spectrometry imaging data with histology image overlays.
# 
# Features:
# - Multiple visualization modes for MSI data
# - Ion selection and mathematical operations
# - Image transformation and overlay capabilities
# - Color palette options and enhancement tools
# - Interactive data tables for ion selection

# Required libraries
if (!require("shiny", quietly = TRUE)) install.packages("shiny")
if (!require("Cardinal", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("Cardinal")
}
if (!require("DT", quietly = TRUE)) install.packages("DT")
if (!require("grid", quietly = TRUE)) install.packages("grid")
if (!require("abind", quietly = TRUE)) install.packages("abind")
if (!require("tools", quietly = TRUE)) install.packages("tools")
if (!require("png", quietly = TRUE)) install.packages("png")
if (!require("jpeg", quietly = TRUE)) install.packages("jpeg")
if (!requireNamespace("viridisLite", quietly=TRUE)) install.packages("viridisLite")
if (!requireNamespace("pals", quietly=TRUE)) install.packages("pals")

library(shiny)
library(Cardinal)
library(DT)
library(grid)
library(abind)
library(tools)
library(png)
library(jpeg)

# Color palette function - handles various color schemes
cpal <- function(name) {
  tryCatch({
    switch(name,
      "Spectral" = rainbow(256),
      "Viridis" = {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::viridis(256)
        } else {
          rainbow(256)
        }
      },
      "Plasma" = {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::plasma(256)
        } else {
          heat.colors(256)
        }
      },
      "Inferno" = {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::inferno(256)
        } else {
          heat.colors(256)
        }
      },
      "Cividis" = {
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          viridisLite::cividis(256)
        } else {
          rainbow(256)
        }
      },
      # Default case - try hcl.colors, fallback to rainbow
      tryCatch(hcl.colors(256, name), error = function(e) rainbow(256))
    )
  }, error = function(e) {
    rainbow(256)  # Ultimate fallback
  })
}

# features within ±ppm around each target
features_in_ppm <- function(obj, targets, ppm){
  mz_all <- mz(obj)
  idxs <- unlist(lapply(targets, function(t){
    tol <- t * ppm * 1e-6
    which(mz_all >= t - tol & mz_all <= t + tol)
  }))
  unique(idxs)
}

# palette (no extra deps required; uses viridisLite if available)
get_palette <- function(name){
  if (name == "gray") return(grDevices::gray.colors(256))
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    switch(name,
      viridis = viridisLite::viridis(256),
      plasma  = viridisLite::plasma(256),
      magma   = viridisLite::magma(256),
      inferno = viridisLite::inferno(256),
      viridisLite::viridis(256)
    )
  } else grDevices::terrain.colors(256)
}

# clip + optional log + raster draw (keeps your Y flip)
draw_clean_image_enhanced <- function(img, pal_name, clip_q = c(0.01, 0.99), log10_scale = FALSE){
  z <- img$mat
  if (log10_scale) z <- log10(pmax(z, 1))

  nz <- as.numeric(z[is.finite(z) & z > 0])
  if (length(nz)) {
    q <- stats::quantile(nz, probs = pmin(pmax(clip_q, 0), 1), na.rm = TRUE)
    z <- pmin(pmax(z, q[1]), q[2])
    rng <- max(q[2] - q[1], .Machine$double.eps)
    z <- (z - q[1]) / rng
  } else {
    z[] <- 0
  }

  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i", bg = "white")

  # Calculate aspect ratio based on grid step sizes
  dx <- if (length(img$xs) > 1) median(diff(img$xs)) else 1
  dy <- if (length(img$ys) > 1) median(diff(img$ys)) else 1
  
  graphics::image(
    x = img$xs, y = img$ys, z = z,  # no transpose needed now
    useRaster = TRUE, axes = FALSE, xlab = "", ylab = "",
    ylim = rev(range(img$ys, na.rm = TRUE)),
    col = get_palette(pal_name),
    asp = (dy * length(img$ys)) / (dx * length(img$xs))  # preserve aspect ratio
  )
}

# Plotting parameters management function
vizi_par <- function(...) {
  if (length(list(...)) == 0) {
    # Return current parameters
    return(par(no.readonly = TRUE))
  } else {
    # Set parameters and return old values
    return(par(...))
  }
}

# Convert intensity vector to spatial image 
vector_to_image <- function(values, coords) {
  # Extract unique x and y coordinates
  xs <- sort(unique(coords$x))
  ys <- sort(unique(coords$y))
  nx <- length(xs); ny <- length(ys)

  # Create a matrix for the image - rows = x, cols = y
  mat <- matrix(0, nrow = nx, ncol = ny)
  
  # Get indices for all points
  x_idx <- match(coords$x, xs)  # 1..nx
  y_idx <- match(coords$y, ys)  # 1..ny
  
  # Fill the matrix with intensity values - use max if multiple spectra map to same pixel
  for (k in seq_along(values)) {
    i <- x_idx[k]; j <- y_idx[k]
    if (!is.na(i) && !is.na(j)) {
      # Use max value when multiple spectra map to the same pixel
      if (is.na(mat[i, j]) || values[k] > mat[i, j]) mat[i, j] <- values[k]
    }
  }
  
  return(list(mat = mat, xs = xs, ys = ys))
}

# Intensity vector for different visualization modes
intensity_vector_for_mode <- function(obj, mode, mz_values = NULL, ppm = 10, norm = "None") {
  I <- intensity(obj)           # features x pixels
  coords <- coord(obj)
  
  # Try to get physical pixel size from metadata (for reference/debugging)
  um_per_px <- tryCatch(pixelSize(obj), error = function(e) NULL)
  if (!is.null(um_per_px)) {
    cat("Physical pixel size from imzML metadata:", um_per_px, "µm/px\n")
  }

  if (mode == "First Ion") {
    vals <- as.numeric(I[1, ])
  } else if (mode == "All Ions (TIC)") {
    vals <- colSums(I, na.rm = TRUE)
  } else { # Custom ions with ppm window(s)
    idxs <- features_in_ppm(obj, mz_values, ppm)
    if (length(idxs) == 0) {
      vals <- rep(0, ncol(I))
    } else {
      vals <- colSums(I[idxs, , drop = FALSE], na.rm = TRUE)
    }
  }

  # normalization
  if (norm == "TIC") {
    tic <- colSums(I, na.rm = TRUE); tic[tic == 0] <- 1
    vals <- vals / tic
  } else if (norm == "Max") {
    m <- max(vals, na.rm = TRUE); if (is.finite(m) && m > 0) vals <- vals / m
  }

  vals[!is.finite(vals)] <- 0
  list(vals = as.numeric(vals), coords = coords)
}

# Read MSI data with progress bar
read_with_progress <- function(path) {
  withProgress(message = "Reading MSI data…", value = 0, {
    incProgress(0.2)
    obj <- Cardinal::readMSIData(path)
    incProgress(0.8)
    obj
  })
}

plot_card_UI <- function(id) {
    
    ns <- NS(id)
    
    col_choices<-c(hcl.pals())
    initial_cols<-c("Spectral", "Cividis", "Viridis", "Inferno", "Plasma", 
                    "Zissou 1", "Purple-Green", "Berlin", "PiYG", "Grays", 
                    "Batlow", "turku", "YlOrRd", "Terrain", 
                    "PrGn", "Green-Brown", "Hawaii", "Cork", "Rocket", "RdYlBu")
    
    col_choices<-c(initial_cols, setdiff(col_choices, initial_cols))
  
    tagList(  
      
        fluidRow(
          tags$head(
            #tags$style("label{font-family: BentonSans Book;}")
            #tags$style("label{font-size: 11px;} ")
          ),
          
          column(3, offset = 0,
                 radioButtons(ns("ion_viz3"), "Image visualization ions",
                                 c("First ion" = "viz_first", 
                                   "All ions (TIC)"="viz_all", 
                                    "Custom single / multiple"="custom")),
                 #numericInput(ns("mz_viz3"), "mz value for visualization",255.2),
                 uiOutput(ns("mz_viz3a")),
                 
                 fluidRow(
                   column(6, selectInput(ns("contrast3"), "Contrast enhancement", c( "none", "histogram",  "adaptive"))),
                   column(6, selectInput(ns("color3"), "Colorscale", col_choices, selected="Spectral"))
                 ),
                 fluidRow(
                   column(6, selectInput(ns("smooth3"), "Smoothing options", c("none", "mean", "gaussian", "bilateral", "adaptive", "diffusion", "guided"))),
                   column(6, checkboxInput(ns("normalize3"), "Scale multiple images?", value = TRUE))
                   ),
                 
                 fluidRow(
                   column(6, checkboxInput(ns("colorkey3"), "Draw color key?", value = TRUE)),
                   column(6, checkboxInput(ns("dark_bg"), "Dark background?", value = FALSE))
                   ),
                   
                 
                 
                 fluidRow(
                   column(6, numericInput(ns("width_im"), "Image plot width (px)", value = 800, step = 50)),
                   column(6, numericInput(ns("height_im"), "Image plot height (px)", value=600, step = 50))
                 ),
                 fluidRow(
                   column(12, checkboxInput(ns("use_1to1_export"), "Use 1:1 pixel export", value = FALSE),
                          helpText("When enabled, exports at exact MSI pixel resolution"))
                 ),
                 fluidRow(
                   column(12, verbatimTextOutput(ns("res_summary")))
                 ),
                 fluidRow(
                   column(12, sliderInput(ns("intensity_threshold"), "Noise threshold percentile", 
                                       min = 0, max = 0.2, value = 0.05, step = 0.01))
                 ),
                 numericInput(ns("ppm"), "m/z tolerance (ppm):", 10, min = 1, step = 1),
                 selectInput(ns("norm"), "Normalize:", c("None","TIC","Max"), selected = "None"),
                 checkboxInput(ns("logScale"), "Log10 scale", value = FALSE),
                 selectInput(ns("pal"), "Palette:", c("gray","viridis","plasma","magma","inferno"), "viridis"),
                 sliderInput(ns("clip"), "Intensity clip (quantiles %):", 0, 100, value = c(1, 99)),
                 checkboxInput(ns("smooth"), "Median smooth (fast)", value = FALSE),
                 checkboxInput(ns("plot_pdata"), "Plot Phenotype data?", value=FALSE),
                 uiOutput(ns("plotpdata")),
                 checkboxInput(ns("expand_fonts"), "Extended font options?", value=FALSE),
                 uiOutput(ns("fonts")),                 
                 checkboxInput(ns("expand_runs"), "Select individual runs for plotting only?", value=FALSE),
                 uiOutput(ns("select_runs")),
                 p("___________________________________________________________________")
                 

                 
          ),
          column(9, uiOutput(ns("plot.window"))
          )  
        ),
        
        uiOutput(ns("spectrum"))
    )
          
  }
  

  plot_card_server <- function(id, overview_peaks_sel, spatialOnly=FALSE, allInputs=NULL) {

    moduleServer(id, function(input, output, session){
      
      #for dynamic UI
      ns = session$ns
      
      graphics.off()
      
      cat("plot_card_server started\n")
      cat("Data class:", class(overview_peaks_sel), "\n")
      cat("Data dimensions:", dim(overview_peaks_sel), "\n")

      #create new overview_peaks_sel object with mean values
      if(is.null(fData(overview_peaks_sel)$mean)) {
        cat("Computing feature summaries...\n")
        tryCatch({
          overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel, verbose=F)
          cat("Feature summaries computed successfully\n")
        }, error = function(e) {
          cat("Error in summarizeFeatures:", e$message, "\n")
          showNotification(paste("Error computing feature summaries:", e$message), type="error")
          return()
        })
        
        if(class(overview_peaks_sel)=="try-error") {
          showNotification("No data available, please check your parameters or dataset", type="error")
          return()
        }
      }
      # browser()
      # #create new overview_peaks_sel object with media values
      # if(is.null(fData(overview_peaks_sel)$non_zero)) {
      #   overview_peaks_sel$non_zero<-Cardinal::summarizeFeatures(overview_peaks_sel, "nnzero")
      #   if(class(overview_peaks_sel)=="try-error") {
      #     showNotification("No data available, please check your parameters or dataset", type="error")
      #     return()
      #   }
      # }
      # 
 
      observe({
        
        output$plot.window <- renderUI({
          
          if(input$ion_viz3=="custom") {
            fluidRow(
              column(12, imageOutput(ns("plot3_pk"))),
              fluidRow(
                column(12, br()),  # Adding space between the image and the table
                column(12, DT::dataTableOutput(ns('tbl')))
              )
            )
          } else {
            fluidRow(
              column(12, imageOutput(ns("plot3_pk")))
            )
          }
          
        })
        
      })
      
      

      
         output$mz_viz3a <- renderUI ({
           req(overview_peaks_sel)
           
           if(input$ion_viz3=="custom") {
             
            
             mz_list<- mz(overview_peaks_sel)
             tbl<-as.data.frame(fData(overview_peaks_sel))
             if(!is.null(tbl$freq)) {
               tbl$freq<-round(tbl$freq, 2)
             }
             if(!is.null(tbl$mean)) {
               tbl$mean<-round(tbl$mean, 1)
             }
             tbl$mz<-round(tbl$mz, 4)
             
             
             updateNumericInput(session, ("width_im"), value=800, step=50)
             updateNumericInput(session, ("height_im"), value=450, step=50)
             
        
             
             #add something here at some point to only select mz, freq, count, and mean columns
             
             
             output$tbl <-DT::renderDataTable({
               tbl
               
             }, selection = "multiple")
  
             list(
               checkboxInput(ns("superpose"), "Superpose images?", value = FALSE),
               selectInput(ns("display_mode"), "Ion math?", 
                          c("none", "sum", "ratio", "subtract", "min", "max", 
                            "mean", "sd", "var", "multiply")),
               numericInput(ns("plusminus_viz3"), "+/- m/z for visualization", 0.05)
             )
           }
         
          })
         
        
    
      
         mz_viz3 <- reactive({
           req(overview_peaks_sel)
           if (is.null(input$tbl_rows_selected) || length(input$tbl_rows_selected) == 0) {
             return(NULL)
           }
           tbl <- as.data.frame(fData(overview_peaks_sel))
           return(tbl[input$tbl_rows_selected, "mz"])
         })

      

      observe({
        
        output$select_runs <- renderUI ({
          req(overview_peaks_sel)
          
          if(input$expand_runs) {
            
            run_list<- unique(run((overview_peaks_sel)))
            
            list(
              selectInput(ns("select_runs"),
                          label= "run selection (plot only)",
                          multiple=TRUE,
                          choices = run_list,
                          selected = run_list)
            )
          }
        })
        
        
      })
      
      # Add MSI resolution and grid information summary
      output$res_summary <- renderText({
        req(overview_peaks_sel)
        c <- coord(overview_peaks_sel)
        xs <- sort(unique(c$x)); ys <- sort(unique(c$y))
        nx <- length(xs); ny <- length(ys)
        dx <- if (length(xs) > 1) median(diff(xs)) else NA
        dy <- if (length(ys) > 1) median(diff(ys)) else NA
        
        # Try to get physical pixel size from metadata
        um <- tryCatch(pixelSize(overview_peaks_sel), error = function(e) NULL)
        
        # If not in metadata, use the user-provided or default value
        if (is.null(um)) {
          um_str <- if (!is.null(allInputs) && !is.null(allInputs$msi_microns_per_pixel)) {
            paste0(allInputs$msi_microns_per_pixel, " µm/px (user input)")
          } else {
            "unknown (not in imzML)"
          }
        } else {
          um_str <- paste0(um, " µm/px (from imzML)")
        }
        
        paste0(
          "MSI grid (pixels): ", nx, " × ", ny, "\n",
          "Grid step (coord units): dx=", round(dx,3), ", dy=", round(dy,3), "\n",
          "Physical pixel size: ", um_str
        )
      })
      
      observe({
        output$fonts <- renderUI ({
          
          if(input$expand_fonts){
            list(
              numericInput(ns("axis_size"), label = "Axis font scaling (%)", min = 0, value=100),
              numericInput(ns("title_size"), label = "Title font scaling (%)", min = 0, value=100),
              numericInput(ns("label_size"), label = "Label font scaling (%)", min = 0, value=100),
              numericInput(ns("subtitle_size"), label = "Subtitle font scaling (%)", min = 0, value=100)
            )
          }
          
          
        })
      })
      
      observe({
        output$plotpdata <- renderUI ({
          
          if(input$plot_pdata){
            list(
              selectInput(ns("pdata_var_plot"), label="pData variable to plot", 
                          choices = colnames(as.data.frame(pData(overview_peaks_sel)))[-c(1:2)]
                          )
            )
          }
          
          
        })
      })

       #add spectrum or not
      output$spectrum<-renderUI({
        
        if(input$ion_viz3=="custom") {
          DT::dataTableOutput(ns('tbl'))
        }

        
        if(spatialOnly==FALSE){
          tagList(
            fluidRow(
              column(6,uiOutput(ns("plot_ranges"))),
              column(6,uiOutput(ns("plot_ranges2"))),
            ),
           
            fluidRow(
              
              column(2,  uiOutput(ns("x_target"))),
              column(2,  numericInput(ns("x_tol"), "+/- m/z window", value=5)),
              column(2,  uiOutput(ns("slider_y_max"))),
              column(3,numericInput(ns("width_ms"), "Plot width (px)", value = 800)),
              column(3,numericInput(ns("height_ms"), "Plot height (px)", value=300))
            ),
            fluidRow(
              column(10, imageOutput(ns("plot4")), style = "margin-bottom: 0px; padding-bottom: 0px;")
              ),
            fluidRow(
              column(3, checkboxInput(ns("calc_ppm"), label = "Show ppm error?", width = '100%')),
              column(3, checkboxInput(ns("show_int"), label = "Show intensity?", width = '100%')),
              column(3, numericInput(ns("show_mz"), label = "# mz values to show?", width = '100%', value=0))
            ),
            checkboxInput(ns("spectrum_expand_fonts"), "Extended font options?", value=FALSE),
            uiOutput(ns("spectrum_fonts")),
          )
        } 
      })
      
      observe({
        output$spectrum_fonts <- renderUI ({
          
          if(input$spectrum_expand_fonts) {
            fluidRow(
              column(
                3,
                numericInput(
                  ns("spectrum_axis_size"),
                  label = "Axis font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_title_size"),
              #     label = "Title font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("spectrum_label_size"),
                  label = "Label font scaling (%)",
                  min = 0,
                  value = 100
                )
              ),
              # column(
              #   2,
              #   numericInput(
              #     ns("spectrum_subtitle_size"),
              #     label = "Subtitle font scaling (%)",
              #     min = 0,
              #     value = 100
              #   )
              # ),
              column(
                3,
                numericInput(
                  ns("linewidth"),
                  label="Linewidth",
                  min=0,
                  value=100
                )
              )
            )
            
          }
          
          
        })
      })
          
       observe({
        output$plot3_pk <- renderImage( {  #plot image in overview after peakpicking / reading file
          
          #req(overview_peaks_sel)
          
          # Add dependency on important parameters to trigger re-rendering
          alpha_val <- if (!is.null(allInputs) && !is.null(allInputs$alpha)) allInputs$alpha else 0.5
          resolution_trigger <- if (!is.null(allInputs)) allInputs$resolution_change_trigger else NULL
          render_trigger <- if (!is.null(allInputs)) allInputs$render_mode_change else NULL
          
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          # Use 1:1 pixel mapping if selected
          if (isTRUE(input$use_1to1_export)) {
            # Get exact MSI grid dimensions for 1:1 pixel mapping
            coords <- coord(overview_peaks_sel)
            nx <- length(unique(coords$x))
            ny <- length(unique(coords$y))
            cat("Using 1:1 pixel mapping for export: ", nx, "×", ny, " pixels\n")
            png(outfile, width = nx, height = ny)
          } else {
            # Use user-specified dimensions
            png(outfile, width = input$width_im, height = input$height_im)
          }
          
          #ion <- switch(input$mode,
          #              "p"=786,
          #             "n"=255.2)
          
          
          if(!is.null(input$select_runs)) {
            overview_peaks_sel <- subsetPixels(overview_peaks_sel, run %in% input$select_runs)
          }
          
          vp_orig<-vizi_par()
          
          #set sizes
          if(input$expand_fonts) {
            req(input$axis_size)
            
            vizi_par(
              cex.axis=input$axis_size/100,
              cex.lab=input$label_size/100,
              cex.main=input$title_size/100,
              cex.sub=input$subtitle_size/100
            )
            #get margins
            # cur_mar<-par()$mar
            # 
            # new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3]+cex.mainp/2, cur_mar[4])
            # 
            # #if drawing colorkey, add a little on the left
            # if(input$colorkey3) {
            #   new_mar[4]=new_mar[4]+cex.axisp
            # }
            # 
            # cur_mgp<-par()$mgp
            # 
            # #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
            # new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
            # 
            
            
            
            
           } else {
              vizi_par(
                cex.axis=1,
                cex.lab=1,
                cex.main=1,
                cex.sub=1
              )

           }
           
          if(input$plot_pdata){
            req(input$pdata_var_plot)
            
            #create list of arguments for image
            arg_list<-list(overview_peaks_sel, 
                       input$pdata_var_plot,
                        key=(input$colorkey3),
                        col=if (requireNamespace("pals", quietly = TRUE)) pals::alphabet() else rainbow(26))
                        
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            plt_tmp<-do.call(Cardinal::image, arg_list)
            
            
            # Cardinal::image(overview_peaks_sel,
            #                         input$pdata_var_plot,
            #                         key=(input$colorkey3),
            #                         #superpose=input$superpose,
            #                         col=if (requireNamespace("pals", quietly = TRUE)) pals::alphabet() else rainbow(26))
            print(plt_tmp,
                                  #cex.axis=req(cex.axisp),
                                  #cex.lab=cex.labp,
                                  #cex.main=cex.mainp,
                                  #cex.sub=cex.subp,
                                  #mar=new_mar,
                                  #mgp=new_mgp
                  )
            
            vizi_par(vp_orig)
          } else if (input$ion_viz3=="viz_all") {
            
            
            mz_range=range(mz(overview_peaks_sel))
            
            #find closest mz value to middle of range
            test_value <- mean(mz_range)
            
            
            # Calculate the absolute differences
            differences <- abs(mz(overview_peaks_sel) - test_value)
            
            # Find the index of the minimum difference
            closest_index <- which.min(differences)
            
            mz_set=mz(overview_peaks_sel)[closest_index]
            tol=max(differences) + differences[closest_index]+1
            
            plusminus=tol
            
            #old way-- may be more memory efficient?
            # smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
            # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
            # 
            # image_command <- paste("Cardinal::image(overview_peaks_sel, 
            #                       mz=mz_set,
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                       enhance_option,
            #                       smoothing_option,",
            #                       scale=input$normalize3,
            #                       #superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            # 
            
            #print(eval(parse(text = image_command)))
            
            
            
            smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
            enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
            
            
            arg_list<-list(overview_peaks_sel, 
                           mz=mz_set,
                           tolerance=round(plusminus,3), 
                           units='mz',
                           col=cpal(input$color3),
                            enhance=enhance_option,
                           smooth=smoothing_option,
                           scale=input$normalize3,
                           #superpose=input$superpose,
                           key=(input$colorkey3))
            
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            # Check if clean rendering mode is enabled for TIC visualization
            if (!is.null(allInputs) && !is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
              # Create clean MSI image with direct data extraction for TIC
              tryCatch({
                cat("Using clean rendering mode for TIC visualization\n")
                
                # Apply median smoothing if requested
                msset_clean <- overview_peaks_sel
                if (isTRUE(input$smooth)) {
                  msset_clean <- tryCatch({
                    Cardinal::smoothSignal(msset_clean, method = "median")
                  }, error = function(e) {
                    cat("Smoothing failed, using original data:", e$message, "\n")
                    msset_clean
                  })
                } else if (input$smooth3 != "none") {
                  msset_clean <- tryCatch({
                    Cardinal::smoothSignal(msset_clean, method = input$smooth3)
                  }, error = function(e) {
                    cat("Smoothing failed, using original data:", e$message, "\n")
                    msset_clean
                  })
                }
                
                # Get intensity vector using our new function
                iv <- tryCatch({
                  intensity_vector_for_mode(
                    msset_clean, 
                    "All Ions (TIC)", 
                    mz_values = NULL, 
                    ppm = if(!is.null(input$ppm)) input$ppm else 10, 
                    norm = if(!is.null(input$norm)) input$norm else "None"
                  )
                }, error = function(e) {
                  cat("Error in intensity_vector_for_mode:", e$message, "\n")
                  # Return default values to avoid errors
                  coords <- coord(msset_clean)
                  list(vals = rep(0, ncol(intensity(msset_clean))), coords = coords)
                })
                
                # Create image matrix
                img_mat <- vector_to_image(iv$vals, iv$coords)
                
                # Use the enhanced drawing function
                q <- if(!is.null(input$clip)) {
                  pmin(pmax(input$clip, 0), 100) / 100
                } else {
                  c(0.01, 0.99)  # Default if not set
                }
                
                draw_clean_image_enhanced(
                  img_mat,
                  pal_name = if(!is.null(input$pal)) input$pal else "viridis",
                  clip_q = c(min(q), max(q)),
                  log10_scale = isTRUE(input$logScale)
                )
                
                # Add title
                title(main = "Total Ion Current (TIC)")
                
                if (input$colorkey3) {
                  # Add color key
                  pal_colors <- get_palette(if(!is.null(input$pal)) input$pal else "viridis")
                  n_colors <- min(10, length(pal_colors))
                  legend_colors <- pal_colors[seq(1, length(pal_colors), length.out = n_colors)]
                  legend("right", 
                         legend = rep("", length(legend_colors)), 
                         fill = legend_colors, 
                         border = NA,
                         bty = "n")
                }
                
                cat("Clean MSI rendering completed for TIC\n")
              }, error = function(e) {
                cat("Error in clean MSI rendering for TIC:", e$message, "\n")
                # Fall back to standard rendering
                print(do.call(Cardinal::image, arg_list))
              })
            } else {
              # Standard Cardinal image rendering
              print(do.call(Cardinal::image, arg_list))
            }
            
            
            
            

            vizi_par(vp_orig)
            
            
          } else if (input$ion_viz3=="custom"){
            
            if(is.null(mz_viz3())){
              ion=mz(overview_peaks_sel[1,])
            } else {
              
              # observeEvent(input$mz_viz3,{
                ion=as.numeric(mz_viz3())
              # })
            }
            
            
            if(!is.null(input$display_mode) && input$display_mode!="none"){
              if(input$display_mode%in%c("min", "max", "min", "mean", "sum", "sd", "var")) { 
          
                
                
                select_vec<-as.character(mz(overview_peaks_sel)) %in% as.character(ion)
                #test to make sure there are 2 or more elements
                if(sum(select_vec)<2){
                  showNotification("At least two ions required for this calculation", type="error")
                  message("At least two ions required for this calculation")
                  return()
                }
                
                sm<-summarizePixels(overview_peaks_sel[select_vec,], stat=c(xic=input$display_mode), as="DataFrame")
                pData(overview_peaks_sel)$xic<-sm$xic
                
                label_txt=paste(input$display_mode, "mz(s)=", paste(ion, collapse=", "))
            
              }else if(input$display_mode=="ratio"){
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for ratio (mz1/mz2)", type="error")
                  message("Exactly two ions required for ratio (mz1/mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  overview_peaks_sel$xic <- (1 + mz1) / (1 + mz2)
                  
                  ion=round(ion, 4)
                  label_txt=paste("ratio of",ion[1],"/",ion[2])
                }
                
              } else if(input$display_mode=="subtract") {
                if(length(ion)!=2){
                  showNotification("Exactly two ions required for subtraction (mz1-mz2)", type="error")
                  message("Exactly two ions required for subtraction (mz1-mz2)")
                  return()
                } else {
                  mz1 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                  mz2 <- spectra(subsetFeatures(overview_peaks_sel, mz=ion[2]))[1,]
                  overview_peaks_sel$xic <- (mz1) - (mz2)
                  ion=round(ion, 4)
                  label_txt=paste("difference of",ion[1],"-",ion[2])
                  
                }
                
                
              }else if(input$display_mode=="multiply"){
                nelements=length(ion)
                xic <- 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[1]))[1,]
                
                for(i in 2:nelements){
                  mz2= 1+spectra(subsetFeatures(overview_peaks_sel, mz=ion[i]))[1,]
                  overview_peaks_sel$xic=(xic) * (mz2)
                  ion=round(ion, 4)
                  label_txt=paste(ion[1],"*",ion[2])
                  }
                }
            
  
              
              
              if(sum(is.na(pData(overview_peaks_sel)$xic))==length(overview_peaks_sel)) {
                showNotification("This calculation does not work!")
                return()
              }
              
              
              
              # smoothing_option <- if (input$smooth3 != "none") paste0(", smooth ='", input$smooth3,"'") else ""
              # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              # 
              plusminus=input$plusminus_viz3
              
            #   image_command <- paste("Cardinal::image(overview_peaks_sel, 'xic',
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                          enhance_option,
            #                          smoothing_option,",
            #                       scale=input$normalize3,
            #                       #superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            #   
              smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
              
              
              arg_list<-list(overview_peaks_sel,
                             'xic',
                             tolerance=round(plusminus,3), 
                             units='mz',
                             col=cpal(input$color3),
                             enhance=enhance_option,
                             smooth=smoothing_option,
                             scale=input$normalize3,
                             #superpose=input$superpose,
                             key=(input$colorkey3))
              
              if(input$dark_bg) {
                arg_list$style <- "dark"
              }
              
              if (requireNamespace("matter", quietly = TRUE)) {
                print(matter::as_facets(do.call(Cardinal::image, arg_list), labels=label_txt))
              } else {
                print(do.call(Cardinal::image, arg_list))
              }
              
              vizi_par(vp_orig)

            } else {
              
              
              
              mz_set=ion
              
              mz_range=range(mz(overview_peaks_sel))
              
              
              # smoothing_option <- if (input$smooth3 != "none") paste0(", smoothing ='", input$smooth3,"'") else ""
              # enhance_option <- if (input$contrast3 != "none") paste0(", enhance ='", input$contrast3,"'") else ""
              # 
              plusminus=input$plusminus_viz3
              
            #   image_command <- paste("Cardinal::image(overview_peaks_sel, 
            #                       mz=mz_set,
            #                       tolerance=round(plusminus,3), 
            #                       units='mz',
            #                       col=cpal(input$color3)",
            #                          enhance_option,
            #                          smoothing_option,",
            #                       scale=input$normalize3,
            #                       superpose=input$superpose,
            #                       key=(input$colorkey3),
            #                       #cex.axis=req(cex.axisp),
            #                       #cex.lab=cex.labp,
            #                       #cex.main=cex.mainp,
            #                       #cex.sub=cex.subp,
            #                       #mar=new_mar,
            #                       #mgp=new_mgp
            # )")
            #   
              
              
              smoothing_option <- if (input$smooth3 != "none")  input$smooth3 else NULL
              enhance_option <- if (input$contrast3 != "none")  input$contrast3 else NULL
              
              # Prepare arguments for Cardinal image
              arg_list <- list(overview_peaks_sel,
                             mz=mz_set,
                             tolerance=round(plusminus,3), 
                             units='mz',
                             col=cpal(input$color3),
                             enhance=enhance_option,
                             smooth=smoothing_option,
                             scale=input$normalize3,
                             superpose=input$superpose,
                             key=(input$colorkey3))
              
              if(input$dark_bg) {
                arg_list$style <- "dark"
              }
              
              # Check for clean rendering mode for custom ion
              if (!is.null(allInputs) && !is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
                # Use direct data extraction approach for better pixel control
                tryCatch({
                  # Apply median smoothing if requested
                  msset_clean <- overview_peaks_sel
                  if (isTRUE(input$smooth)) {
                    msset_clean <- tryCatch({
                      Cardinal::smoothSignal(msset_clean, method = "median")
                    }, error = function(e) {
                      cat("Smoothing failed, using original data:", e$message, "\n")
                      msset_clean
                    })
                  } else if (input$smooth3 != "none") {
                    msset_clean <- tryCatch({
                      Cardinal::smoothSignal(msset_clean, method = input$smooth3)
                    }, error = function(e) {
                      cat("Smoothing failed, using original data:", e$message, "\n")
                      msset_clean
                    })
                  }
                  
                  # Get intensity vector using our new function
                  iv <- tryCatch({
                    intensity_vector_for_mode(
                      msset_clean, 
                      "Custom Ion(s)", 
                      mz_values = mz_set, 
                      ppm = if(!is.null(input$ppm)) input$ppm else 10, 
                      norm = if(!is.null(input$norm)) input$norm else "None"
                    )
                  }, error = function(e) {
                    cat("Error in intensity_vector_for_mode:", e$message, "\n")
                    # Return default values to avoid errors
                    coords <- coord(msset_clean)
                    list(vals = rep(0, ncol(intensity(msset_clean))), coords = coords)
                  })
                  
                  # Create image matrix
                  img_mat <- vector_to_image(iv$vals, iv$coords)
                  
                  # Use the enhanced drawing function
                  q <- if(!is.null(input$clip)) {
                    pmin(pmax(input$clip, 0), 100) / 100
                  } else {
                    c(0.01, 0.99)  # Default if not set
                  }
                  
                  draw_clean_image_enhanced(
                    img_mat,
                    pal_name = if(!is.null(input$pal)) input$pal else "viridis",
                    clip_q = c(min(q), max(q)),
                    log10_scale = isTRUE(input$logScale)
                  )
                  
                  # Add title
                  if (length(mz_set) == 1) {
                    title(main = paste("m/z", round(mz_set, 4)))
                  } else if (length(mz_set) > 1) {
                    title(main = paste("Multiple m/z values"))
                  }
                  
                  if (input$colorkey3) {
                    # Add color key
                    pal_colors <- get_palette(if(!is.null(input$pal)) input$pal else "viridis")
                    n_colors <- min(10, length(pal_colors))
                    legend_colors <- pal_colors[seq(1, length(pal_colors), length.out = n_colors)]
                    legend("right", 
                           legend = rep("", length(legend_colors)), 
                           fill = legend_colors, 
                           border = NA,
                           bty = "n")
                  }
                  
                  cat("Clean MSI rendering completed for custom ion\n")
                }, error = function(e) {
                  cat("Error in clean MSI rendering for custom ion:", e$message, "\n")
                  # Fall back to standard rendering
                  print(do.call(Cardinal::image, arg_list))
                })
              } else {
                # Standard Cardinal image rendering
                print(do.call(Cardinal::image, arg_list))
              }
              
              vizi_par(vp_orig)

              
            }
          } else if (input$ion_viz3=="viz_first") {
            
            tol=0.05
            
              # Check for clean rendering mode
              # Re-enable clean rendering mode
              if (!is.null(allInputs) && !is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
                # Create clean MSI image with direct data extraction
                tryCatch({
                  # Get the first m/z value
                  selected_mz <- mz(overview_peaks_sel)[1]
                  cat("Using first m/z value:", selected_mz, "\n")
                  
                  # Apply median smoothing if requested
                  msset_clean <- overview_peaks_sel
                  if (isTRUE(input$smooth)) {
                    msset_clean <- tryCatch({
                      Cardinal::smoothSignal(msset_clean, method = "median")
                    }, error = function(e) {
                      cat("Smoothing failed, using original data:", e$message, "\n")
                      msset_clean
                    })
                  } else if (input$smooth3 != "none") {
                    msset_clean <- tryCatch({
                      Cardinal::smoothSignal(msset_clean, method = input$smooth3)
                    }, error = function(e) {
                      cat("Smoothing failed, using original data:", e$message, "\n")
                      msset_clean
                    })
                  }
                  
                  # Get intensity vector using our new function
                  iv <- tryCatch({
                    intensity_vector_for_mode(
                      msset_clean, 
                      "First Ion", 
                      mz_values = NULL, 
                      ppm = if(!is.null(input$ppm)) input$ppm else 10, 
                      norm = if(!is.null(input$norm)) input$norm else "None"
                    )
                  }, error = function(e) {
                    cat("Error in intensity_vector_for_mode:", e$message, "\n")
                    # Return default values to avoid errors
                    coords <- coord(msset_clean)
                    list(vals = rep(0, ncol(intensity(msset_clean))), coords = coords)
                  })
                  
                  # Create image matrix
                  img_mat <- vector_to_image(iv$vals, iv$coords)
                  
                  # Use the enhanced drawing function
                  q <- if(!is.null(input$clip)) {
                    pmin(pmax(input$clip, 0), 100) / 100
                  } else {
                    c(0.01, 0.99)  # Default if not set
                  }
                  
                  draw_clean_image_enhanced(
                    img_mat,
                    pal_name = if(!is.null(input$pal)) input$pal else "viridis",
                    clip_q = c(min(q), max(q)),
                    log10_scale = isTRUE(input$logScale)
                  )
                  
                  # Add title
                  title(main = paste("m/z", round(selected_mz, 4)))
                  
                  if (input$colorkey3) {
                    # Add color key
                    pal_colors <- get_palette(if(!is.null(input$pal)) input$pal else "viridis")
                    n_colors <- min(10, length(pal_colors))
                    legend_colors <- pal_colors[seq(1, length(pal_colors), length.out = n_colors)]
                    legend("right", 
                           legend = rep("", length(legend_colors)), 
                           fill = legend_colors, 
                           border = NA,
                           bty = "n")
                  }
                  
                  cat("Clean MSI rendering completed\n")
                }, error = function(e) {
                  cat("Error in clean MSI rendering:", e$message, "\n")
                  # Fall back to standard rendering
                  image_command <- Cardinal::image(overview_peaks_sel, 
                                      col=cpal(input$color3),
                                      scale=input$normalize3,
                                      key=(input$colorkey3))
                  print(image_command)
                })
              } else {
                # Standard Cardinal image rendering
                image_command <- Cardinal::image(overview_peaks_sel, 
                                    col=cpal(input$color3),
                                    scale=input$normalize3,
                                    key=(input$colorkey3))
                print(image_command)
              }            
              vizi_par(vp_orig)
            
          }
          
          # Apply histology overlay if images are uploaded
          # This must be done BEFORE dev.off()
          if (!is.null(allInputs$histology_upload) && 
              !is.null(allInputs$msi_upload)) {
            
            tryCatch({
              ensure_3d <- function(image) {
                if (length(dim(image)) == 2) {
                  image <- abind::abind(image, image, image, along = 3)
                }
                else if (length(dim(image)) > 3) {
                  image <- image[,,,1]
                }
                return(image)
              }
              
              alpha_display <- if (!is.null(allInputs) && !is.null(allInputs$alpha)) allInputs$alpha else 0.5
              cat("Applying histology overlay with alpha =", alpha_display, "\n")
              
              # Check if using sample data (helpful for testing)
              is_sample_data <- FALSE
              if (is.character(allInputs$histology_upload$name) && 
                  grepl("sample", tolower(allInputs$histology_upload$name))) {
                cat("Detected sample histology data in filename\n")
                is_sample_data <- TRUE
                
                # For sample data, try loading from current directory as fallback
                if (file.exists("sample_histology.png")) {
                  cat("Found sample_histology.png in current directory, will try as fallback if needed\n")
                }
              }
              
              # Load histology image with improved file type detection
              histology_path <- allInputs$histology_upload$datapath
              
              # More robust file type detection
              file_type <- tolower(tools::file_ext(histology_path))
              cat("Histology file type detected:", file_type, "\n")
              
              # Check if file exists and is readable
              if (!file.exists(histology_path)) {
                if (!is.null(allInputs$is_sample_data) && allInputs$is_sample_data && file.exists("sample_histology.png")) {
                  cat("Using sample_histology.png as fallback\n")
                  histology_path <- "sample_histology.png"
                  file_type <- "png"
                } else {
                  stop(paste("Histology file not found:", histology_path))
                }
              }
              
              # Better error handling for image loading
              tryCatch({
                if (file_type %in% c("png")) {
                  histology_image <- png::readPNG(histology_path)
                  cat("Loaded PNG image, dimensions:", paste(dim(histology_image), collapse="×"), "\n")
                } else if (file_type %in% c("jpg", "jpeg")) {
                  histology_image <- jpeg::readJPEG(histology_path)
                  cat("Loaded JPEG image, dimensions:", paste(dim(histology_image), collapse="×"), "\n")
                } else {
                  # Try to guess based on file content
                  cat("Attempting to determine file type from content...\n")
                  # Try PNG first
                  histology_image <- tryCatch(
                    png::readPNG(histology_path),
                    error = function(e) {
                      # If PNG fails, try JPEG
                      tryCatch(
                        jpeg::readJPEG(histology_path),
                        error = function(e2) {
                          stop(paste("Unsupported file type or corrupted image file:", file_type))
                        }
                      )
                    }
                  )
                  cat("Successfully loaded image by content inspection\n")
                }
              }, error = function(e) {
                stop(paste("Error loading histology image:", e$message))
              })
              
              # Fix dimensions properly depending on array structure
              if (length(dim(histology_image)) == 2) {
                # Grayscale image - convert to RGB
                histology_image <- array(rep(histology_image, 3), 
                                        dim = c(dim(histology_image), 3))
              } else if (length(dim(histology_image)) == 3) {
                if (dim(histology_image)[3] == 4) {
                  # Already has alpha channel - strip it
                  histology_image <- histology_image[,,1:3]
                } else if (dim(histology_image)[3] != 3) {
                  # Weird number of channels - force to RGB
                  histology_image <- histology_image[,,1:min(3, dim(histology_image)[3])]
                  # If less than 3 channels, duplicate the last one
                  if (dim(histology_image)[3] < 3) {
                    last_channel <- histology_image[,,dim(histology_image)[3]]
                    for (i in (dim(histology_image)[3]+1):3) {
                      histology_image <- abind::abind(histology_image, last_channel, along=3)
                    }
                  }
                }
              } else if (length(dim(histology_image)) > 3) {
                # Too many dimensions - flatten to 3D
                histology_image <- histology_image[,,,1]
                if (length(dim(histology_image)) == 2) {
                  histology_image <- array(rep(histology_image, 3), 
                                         dim = c(dim(histology_image), 3))
                }
              }
              
              # Add alpha channel with user-defined transparency
              alpha_value <- if (!is.null(allInputs) && !is.null(allInputs$alpha) && is.numeric(allInputs$alpha)) {
                allInputs$alpha
              } else {
                0.7  # Use a different default than 0.5 to show it's not hardcoded
              }
              cat("Using alpha transparency value:", alpha_value, "\n")
              
              alpha_channel <- array(alpha_value, 
                                     dim = c(dim(histology_image)[1], 
                                             dim(histology_image)[2]))
              histology_image <- abind::abind(histology_image, alpha_channel, along = 3)
              
              # Get base scale factors from UI sliders
              scalex <- allInputs$scalex
              scaley <- allInputs$scaley
              
              # Log current resolution and dimension values for debugging
              cat("Current resolution and dimension values in plot function:\n")
              cat(" - Histology μm/pixel:", 
                  if(!is.null(allInputs$histology_microns_per_pixel)) allInputs$histology_microns_per_pixel else "NULL", "\n")
              cat(" - MSI μm/pixel:", 
                  if(!is.null(allInputs$msi_microns_per_pixel)) allInputs$msi_microns_per_pixel else "NULL", "\n")
              cat(" - Histology dimensions (px):", 
                  if(!is.null(allInputs$histology_pixel_width)) allInputs$histology_pixel_width else "NULL", "×", 
                  if(!is.null(allInputs$histology_pixel_height)) allInputs$histology_pixel_height else "NULL", "\n")
              cat(" - MSI dimensions (px):", 
                  if(!is.null(allInputs$msi_pixel_width)) allInputs$msi_pixel_width else "NULL", "×", 
                  if(!is.null(allInputs$msi_pixel_height)) allInputs$msi_pixel_height else "NULL", "\n")
              cat(" - Auto-scale enabled:", 
                  if(!is.null(allInputs$auto_scale_resolution)) allInputs$auto_scale_resolution else "NULL", "\n")
              cat(" - Scale target:", 
                  if(!is.null(allInputs$scale_target)) allInputs$scale_target else "NULL", "\n")
              
              # Calculate combined scaling based on both pixel dimensions and resolution
              pixel_dimension_scale_x <- 1
              pixel_dimension_scale_y <- 1
              
              # Only apply pixel dimension scaling if all values are available
              if (!is.null(allInputs$histology_pixel_width) && !is.null(allInputs$msi_pixel_width) &&
                  !is.null(allInputs$histology_pixel_height) && !is.null(allInputs$msi_pixel_height) &&
                  !is.na(allInputs$histology_pixel_width) && !is.na(allInputs$msi_pixel_width) && 
                  !is.na(allInputs$histology_pixel_height) && !is.na(allInputs$msi_pixel_height) && 
                  allInputs$histology_pixel_width > 0 && allInputs$msi_pixel_width > 0 &&
                  allInputs$histology_pixel_height > 0 && allInputs$msi_pixel_height > 0) {
                  
                  # Calculate pixel dimension ratios
                  pixel_dimension_scale_x <- allInputs$histology_pixel_width / allInputs$msi_pixel_width
                  pixel_dimension_scale_y <- allInputs$histology_pixel_height / allInputs$msi_pixel_height
                  
                  cat("Pixel dimension scale factors - X:", pixel_dimension_scale_x, "Y:", pixel_dimension_scale_y, "\n")
              } else {
                  cat("Pixel dimension scaling not applied (missing or invalid dimensions)\n")
              }
              
              # Apply resolution-based scaling if enabled
              if (!is.null(allInputs$auto_scale_resolution) && allInputs$auto_scale_resolution && 
                  !is.null(allInputs$histology_microns_per_pixel) && 
                  !is.null(allInputs$msi_microns_per_pixel) && 
                  !is.na(allInputs$histology_microns_per_pixel) && 
                  !is.na(allInputs$msi_microns_per_pixel) && 
                  allInputs$histology_microns_per_pixel > 0 && 
                  allInputs$msi_microns_per_pixel > 0) {
                
                # Pick the target grid in µm/px based on selection
                if (allInputs$scale_target == "msi") {
                  target_um_per_px <- allInputs$msi_microns_per_pixel
                } else if (allInputs$scale_target == "histology") {
                  target_um_per_px <- allInputs$histology_microns_per_pixel
                } else { 
                  # "best_fit": geometric-mean grid is a sensible compromise
                  target_um_per_px <- sqrt(allInputs$msi_microns_per_pixel * allInputs$histology_microns_per_pixel)
                }
                
                # Calculate per-image scale factors (unitless, applied to pixel width/height)
                scale_histology <- allInputs$histology_microns_per_pixel / target_um_per_px
                scale_msi <- allInputs$msi_microns_per_pixel / target_um_per_px
                
                # For reference, calculate the old way too to log difference
                resolution_ratio <- allInputs$msi_microns_per_pixel / allInputs$histology_microns_per_pixel
                
                cat("Target µm/pixel:", target_um_per_px, "\n")
                cat("Resolution-based scale factors - Histology:", scale_histology, "MSI:", scale_msi, "\n")
                cat("Resolution ratio (MSI:Histology):", resolution_ratio, "\n")
                
                # Combine resolution scaling with pixel dimension scaling
                final_scale_x <- scalex * scale_histology * pixel_dimension_scale_x
                final_scale_y <- scaley * scale_histology * pixel_dimension_scale_y
                
                cat("Final combined scale factors - X:", final_scale_x, "Y:", final_scale_y, "\n")
                
                # Apply appropriate scaling to histology image
                scalex <- final_scale_x
                scaley <- final_scale_y
              } else {
                # Even if auto-scale is disabled, still apply pixel dimension scaling
                scalex <- scalex * pixel_dimension_scale_x
                scaley <- scaley * pixel_dimension_scale_y
                
                cat("Applied pixel dimension scaling only (auto-scale disabled)\n")
                cat("Final scale factors - X:", scalex, "Y:", scaley, "\n")
              }
              
              # Safety check for very large images (prevent memory issues)
              img_size <- prod(dim(histology_image)[1:2])
              if (img_size > 5000000) { # 5 million pixels threshold
                # Simple downsample for very large images
                downsample_factor <- sqrt(5000000 / img_size)
                new_dims <- round(dim(histology_image)[1:2] * downsample_factor)
                
                # Perform simple downsampling
                rows <- round(seq(1, dim(histology_image)[1], length.out = new_dims[1]))
                cols <- round(seq(1, dim(histology_image)[2], length.out = new_dims[2]))
                
                downsampled <- array(0, c(new_dims, dim(histology_image)[3]))
                for (i in 1:dim(histology_image)[3]) {
                  downsampled[,,i] <- histology_image[rows, cols, i]
                }
                histology_image <- downsampled
                cat("Large image downsampled to:", paste(new_dims, collapse="x"), "pixels\n")
              }
              
              # Ensure scale values are valid
              scalex <- ifelse(is.finite(scalex) && scalex > 0, scalex, 1)
              scaley <- ifelse(is.finite(scaley) && scaley > 0, scaley, 1)
              
              # Create and draw histology overlay
              histology_grob <- rasterGrob(histology_image, interpolate = TRUE)
              
              # We'll use normalized coordinates (npc) with appropriate scaling
              # This approach works reliably across different device sizes
              width_value <- 1 * scalex
              height_value <- 1 * scaley
              
              # Get the ratio between histology and MSI dimensions to inform scaling
              # but we'll still use normalized coordinates for robustness
              if (!is.null(allInputs$histology_pixel_width) && !is.null(allInputs$msi_pixel_width) &&
                  !is.null(allInputs$histology_pixel_height) && !is.null(allInputs$msi_pixel_height) &&
                  !is.na(allInputs$histology_pixel_width) && !is.na(allInputs$msi_pixel_width) && 
                  !is.na(allInputs$histology_pixel_height) && !is.na(allInputs$msi_pixel_height) && 
                  allInputs$histology_pixel_width > 0 && allInputs$msi_pixel_width > 0 &&
                  allInputs$histology_pixel_height > 0 && allInputs$msi_pixel_height > 0) {
                  
                  # Log calculated dimensions
                  cat("Using dimension ratios to inform scaling:\n")
                  cat("Histology dimensions:", allInputs$histology_pixel_width, "×", allInputs$histology_pixel_height, "\n")
                  cat("MSI dimensions:", allInputs$msi_pixel_width, "×", allInputs$msi_pixel_height, "\n")
                  cat("Scaling applied - X:", scalex, "Y:", scaley, "\n")
              }
              
              # Get the current device size to properly center the image
              current_dev <- dev.cur()
              dev_size <- dev.size("in")  # Get size in inches which is more reliable
              dev_width <- dev_size[1]
              dev_height <- dev_size[2]
              
              cat("Current device size:", dev_width, "×", dev_height, "inches\n")
              
              # Apply transformations to the grob using normalized coordinates (npc)
              # This is the most reliable approach across different devices
              # Get translation values with safe defaults
              translate_x <- if (!is.null(allInputs) && !is.null(allInputs$translate_x)) allInputs$translate_x else 0
              translate_y <- if (!is.null(allInputs) && !is.null(allInputs$translate_y)) allInputs$translate_y else 0
              rotate_val <- if (!is.null(allInputs) && !is.null(allInputs$rotate)) allInputs$rotate else 0
              
              histology_grob <- editGrob(
                histology_grob,
                vp = viewport(
                  x = unit(0.5, "npc") + unit(translate_x, "mm"),
                  y = unit(0.5, "npc") + unit(translate_y, "mm"),
                  angle = rotate_val,
                  width = unit(width_value, "npc"),  # Normalized coordinates with scaling
                  height = unit(height_value, "npc"), # Normalized coordinates with scaling
                  just = c("center", "center")
                )
              )
              
              # Draw the overlay with more detailed logging
              cat("Drawing histology overlay with transform parameters:\n")
              cat(" - Translation: X=", translate_x, "mm, Y=", translate_y, "mm\n")
              cat(" - Rotation:", rotate_val, "degrees\n")
              cat(" - Scale: X=", scalex, ", Y=", scaley, "\n")
              
              grid.draw(histology_grob)
              
            }, error = function(e) {
              # Basic error handling
              cat("Error in histology overlay:", e$message, "\n")
            })
          }
          
          dev.off()
          
          
          # Return a list containing the filename
          list(src = outfile,
               contentType = 'image/png',
               width = input$width_im,
               height = input$height_im,
               alt = "This is alternate text")
        }, deleteFile = TRUE)
       }) 
       
       
       if(spatialOnly==FALSE) {
         observe({  
          output$plot_ranges<- renderUI( {
            
            req(overview_peaks_sel)
            
            a<-Cardinal::plot(overview_peaks_sel)
            
  
                  sliderInput(ns("mass_range_plot"), 
                              label = p("m/z range for MS plot (X)"), 
                              min = round(a$channels$x$limits[1]),
                              max = round(a$channels$x$limits[2]), 
                              value = round(a$channels$x$limits), 
                              step = NULL,
                              round = TRUE)
  
            
          })
         })
         
         observe({
          output$plot_ranges2<- renderUI( {
            req(overview_peaks_sel)
            #browser()
            #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
            a<-Cardinal::plot(overview_peaks_sel, "mean")
            
            
           
            
            sliderInput (ns("int_range_plot"), 
                        label = p("intensity range for MS plot (Y)"), 
                        min = round(a$channels$y$limits[1],0),
                        max = round(a$channels$y$limits[2],0)*1.05, 
                        value = req(input$param_numeric), #round(a$par$ylim,0), 
                        step = a$channels$y$limits[2]/20,
                        round = TRUE
            )
            
                        
                        
            
          })
         })
         
         
         #set center of observed spectrum if using custom ion visualization
         observe({
           output$x_target <- renderUI({
             
             if(input$ion_viz3!="custom") {
               numericInput(ns("x_target"), "Center m/z value", value = NULL)
             } else {
               numericInput(ns("x_target"), "Center m/z value", value = mz_viz3())
             }
             
           })
         })
         
        
         
         
         
         observe( {
           
           output$slider_y_max <- renderUI({
             
             req(overview_peaks_sel)
            
             
             a<-Cardinal::plot(overview_peaks_sel, "mean")
             
             
             numericInput(ns("param_numeric"),
                          "Manual y-axis intensity max value",
                          min = round(a$channels$y$limits[1],0),
                          max = round(a$channels$y$limits[2],0),
                          value = round(a$channels$y$limits[2],0)
             )
           })
         })
               
          
         if(!is.null(overview_peaks_sel)) {
           
           updateSliderInput(session, ns("int_range_plot"), value = c(round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits,0) 
                                                                      ))
           updateNumericInput(session, ns("param_numeric"), value = round(Cardinal::plot(req(overview_peaks_sel))$channels$y$limits[2],0))
          }
        
        
        
         
         observe({
           output$plot4 <- renderImage( {
             
             
             req(overview_peaks_sel)
             req(input$mass_range_plot)
             
             # A temp file to save the output.
             # This file will be removed later by renderImage
             outfile <- tempfile(fileext = '.png')
             
             png(outfile, width = input$width_ms, height = input$height_ms)
             
            
             
             if(length(input$int_range_plot)==1) {
               ylim=c(0, input$int_range_plot)
             } else {
               ylim = input$int_range_plot
             }
             
             #change xlimits based on custom ion or not
             if(is.null(input$x_target) || is.na(input$x_target)){
               xlim=input$mass_range_plot
               
               overview_peaks_sel_plot<-overview_peaks_sel
               
             } else {
               xlim=c(input$x_target-input$x_tol, input$x_target+input$x_tol)
               
               #subsetFeatures to only include mz values within range
               overview_peaks_sel_plot<-subsetFeatures(overview_peaks_sel, mz > xlim[1] & mz < xlim[2])
               
             }
             
             if(!is.finite(xlim[1])){
               xlim=input$mass_range_plot
             }
             
             vp_orig <- vizi_par()
             if(input$spectrum_expand_fonts) {
               req(input$spectrum_axis_size)
               
               vizi_par(
                 cex.axis = req(input$spectrum_axis_size)/100,
                 cex.lab = input$spectrum_label_size/100,
                 cex.main = input$spectrum_label_size/100,
                 lwd = input$linewidth/100,
                 mar = c(0, 0, 1, 1)
               )
               
               
               #get margins
               #cur_mar<-par()$mar
               
               #new_mar<-c(cur_mar[1]+cex.labp/8, cur_mar[2]+cex.labp/2, cur_mar[3], cur_mar[4])
               
               #cur_mgp<-par()$mgp
               
               #new_mgp<-c(cur_mgp[1]+cex.labp/8, cur_mgp[2], 0)
               #new_mgp<-c(3+max(new_mar[1:2])/20, 1,0)
               
               #lwd=lwdp
                          
               
               
               
               
             } else {
               vizi_par(
                 cex.axis=1,
                 cex.lab=1,
                 cex.main=1,
                 cex.sub=1
               )
             }
             
             #browser()
             #overview_peaks_sel<-Cardinal::summarizeFeatures(overview_peaks_sel)
             
             p1<-Cardinal::plot(overview_peaks_sel_plot,
                                xlim=xlim,
                                ylim =ylim,
                                #cex.axis=req(cex.axisp),
                                #cex.lab=cex.labp,
                                #cex.main=cex.mainp,
                                #cex.sub=cex.subp,
                                #lwd=lwdp,
                                # mar=new_mar, 
                                # mgp=new_mgp, 
                                "mean",
                                annPeaks=input$show_mz,
                                free="y")
             print(p1)
             vizi_par(vp_orig)
             
             #check for ppm calc
             if(input$calc_ppm) {
               
               dat=overview_peaks_sel_plot
               x=mz(dat)
               targ_mz<-req(input$x_target)
               x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
               
               ppm_error<- round(1e6*(x_sel-targ_mz)/targ_mz, 2)
               
               p1_coord<-p1[[1]][[1]]$marks$peaks$encoding
               
               y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
               
               if(length(ppm_error)==0){
                 showNotification("No ppm error calculated, are there any peaks?", type="warning")
                 return()
               } else {
                print(text(x=x_sel, y=y_labs+ylim[2]*.25, req(ppm_error)))
               }
               
             }
             
             
             
             if(input$show_int) {
              
               
               ###
               p1_coord<-p1[[1]][[1]]$marks$peaks$encoding
               
               
               
               dat=overview_peaks_sel_plot
               x=mz(dat)
               
               x_sel<-subset(x, x>=xlim[1] & x<= xlim[2])
               
               
               
               y_labs<-p1_coord$y[p1_coord$x %in% x_sel]
               
               if(length(y_labs)==0){
                 showNotification("No intensities found, are there any peaks?", type="warning")
                 return()
               } else {
                 print(text(x=x_sel, y=y_labs+ylim[2]*.15, req(round(y_labs, 0))))
               }
               
             }
             
             

             
             dev.off()
             
             # Return a list containing the filename
             list(src = outfile,
                  contentType = 'image/png',
                  width = input$width_ms,
                  height = input$height_ms,
                  alt = "This is alternate text")
           }, deleteFile = TRUE)
           
         })
       }
  })
        
  
}

