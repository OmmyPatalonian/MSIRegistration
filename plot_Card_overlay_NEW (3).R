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
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::viridis(256)
        } else {
          rainbow(256)
        }
      },
      "Plasma" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::plasma(256)
        } else {
          heat.colors(256)
        }
      },
      "Inferno" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::inferno(256)
        } else {
          heat.colors(256)
        }
      },
      "Cividis" = {
        if (requireNamespace("viridis", quietly = TRUE)) {
          viridis::cividis(256)
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
                   column(12, sliderInput(ns("intensity_threshold"), "Noise threshold percentile", 
                                       min = 0, max = 0.2, value = 0.05, step = 0.01))
                 ),
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
          alpha_val <- allInputs$alpha
          resolution_trigger <- allInputs$resolution_change_trigger
          render_trigger <- allInputs$render_mode_change
          
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          png(outfile, width = input$width_im, height = input$height_im)
          
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
                        col=pals::alphabet())
                        
            if(input$dark_bg) {
              arg_list$style <- "dark"
            }
            
            plt_tmp<-do.call(Cardinal::image, arg_list)
            
            
            # Cardinal::image(overview_peaks_sel,
            #                         input$pdata_var_plot,
            #                         key=(input$colorkey3),
            #                         #superpose=input$superpose,
            #                         col=pals::alphabet())
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
            if (!is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
              # Create clean MSI image with direct data extraction for TIC
              tryCatch({
                cat("Using clean rendering mode for TIC visualization\n")
                
                # Extract grid coordinates directly
                coords <- coord(overview_peaks_sel)
                
                # Apply preprocessing for cleaner visualization
                msset_clean <- overview_peaks_sel
                if (input$smooth3 != "none") {
                  msset_clean <- tryCatch({
                    Cardinal::smoothSignal(msset_clean, method = input$smooth3)
                  }, error = function(e) {
                    cat("Smoothing failed, using original data:", e$message, "\n")
                    msset_clean
                  })
                }
                
                # For TIC, we want all features within the m/z range
                tol_value <- round(plusminus, 3)
                
                # Subset features based on m/z range
                sub_data <- Cardinal::subsetFeatures(msset_clean, 
                                                    mz >= mz_set - tol_value & 
                                                    mz <= mz_set + tol_value)
                
                # Calculate TIC for each pixel (sum of all intensities)
                intensities <- Cardinal::intensity(sub_data)
                pixel_values <- apply(intensities, 2, sum, na.rm = TRUE)
                
                # Remove NA and Inf values
                pixel_values[is.na(pixel_values) | is.infinite(pixel_values)] <- 0
                
                # Apply contrast enhancement if selected
                if (input$contrast3 == "histogram") {
                  # Histogram equalization for better contrast
                  pixel_values <- rank(pixel_values) / length(pixel_values) * max(pixel_values, na.rm = TRUE)
                } else if (input$contrast3 == "adaptive") {
                  # Simple adaptive contrast enhancement
                  q_min <- quantile(pixel_values, 0.05, na.rm = TRUE)
                  q_max <- quantile(pixel_values, 0.95, na.rm = TRUE)
                  pixel_values <- (pixel_values - q_min) / (q_max - q_min) * max(pixel_values, na.rm = TRUE)
                  pixel_values[pixel_values < 0] <- 0
                }
                
                # Filter low-intensity values (likely noise)
                nonzero_vals <- pixel_values[pixel_values > 0]
                if (length(nonzero_vals) > 0) {
                  intensity_threshold <- quantile(nonzero_vals, input$intensity_threshold, na.rm = TRUE)
                  pixel_values[pixel_values < intensity_threshold] <- 0
                  
                  # Generate colormap with selected palette
                  colormap <- cpal(input$color3)
                  
                  # Scale values to colormap range
                  scaled_values <- pixel_values
                  if (max(scaled_values, na.rm = TRUE) > 0) {
                    scaled_values <- (scaled_values - min(nonzero_vals, na.rm = TRUE)) / 
                                     (max(scaled_values, na.rm = TRUE) - min(nonzero_vals, na.rm = TRUE))
                    scaled_values[scaled_values < 0] <- 0
                    scaled_values[scaled_values > 1] <- 1
                  }
                  
                  # Map to color indices
                  color_indices <- round(scaled_values * (length(colormap) - 1)) + 1
                  color_indices[color_indices < 1 | is.na(color_indices)] <- 1
                  image_colors <- colormap[color_indices]
                  
                  # Set transparent color for zero values
                  image_colors[pixel_values == 0] <- NA
                  
                  # Render clean MSI image
                  plot.new()
                  plot.window(xlim = range(coords$x, na.rm = TRUE),
                             ylim = range(coords$y, na.rm = TRUE))
                  
                  # Adjust background color if needed
                  if (input$dark_bg) {
                    par(bg = "black")
                  }
                  
                  # Plot each pixel with appropriate color and size
                  points(coords$x, coords$y, 
                        pch = 15, 
                        col = image_colors,
                        cex = 1.5)  # Adjust pixel size as needed
                  
                  # Add axes and title with appropriate style
                  box()
                  title_color <- if(input$dark_bg) "white" else "black"
                  title(main = "Total Ion Current (TIC)", col.main = title_color)
                  
                  if (input$colorkey3) {
                    # Add color key with better legend formatting
                    n_colors <- min(10, length(unique(na.omit(image_colors))))
                    legend_colors <- colormap[seq(1, length(colormap), length.out = n_colors)]
                    legend("right", 
                           legend = rep("", length(legend_colors)), 
                           fill = legend_colors, 
                           border = NA,
                           bty = "n",
                           text.col = title_color)
                  }
                  
                  cat("Clean MSI rendering completed for TIC\n")
                } else {
                  cat("No non-zero intensity values found for TIC\n")
                  # Fall back to standard rendering
                  print(do.call(Cardinal::image, arg_list))
                }
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
              
              print(matter::as_facets(do.call(Cardinal::image, arg_list), labels=label_txt))
              
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
              if (!is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
                # Use direct data extraction approach for better pixel control
                tryCatch({
                  # Extract grid coordinates directly from MSI data
                  coords <- coord(overview_peaks_sel)
                  
                  # Identify m/z indices for selected ions
                  mz_indices <- sapply(mz_set, function(mz_value) {
                    which.min(abs(mz(overview_peaks_sel) - mz_value))
                  })
                  
                  # Apply preprocessing for cleaner visualization
                  msset_clean <- overview_peaks_sel
                  if (input$smooth3 != "none") {
                    msset_clean <- tryCatch({
                      Cardinal::smoothSignal(msset_clean, method = input$smooth3)
                    }, error = function(e) {
                      cat("Smoothing failed, using original data:", e$message, "\n")
                      msset_clean
                    })
                  }
                  
                  # Extract intensity values for each selected m/z
                  all_intensities <- list()
                  for (i in seq_along(mz_indices)) {
                    idx <- mz_indices[i]
                    # Get intensities for this m/z value with tolerance
                    mz_value <- mz(overview_peaks_sel)[idx]
                    sub_data <- Cardinal::subsetFeatures(msset_clean, 
                                                        mz >= mz_value - round(plusminus, 3) & 
                                                        mz <= mz_value + round(plusminus, 3))
                    
                    if (input$superpose && i > 1) {
                      # For superposition, combine with previous intensities
                      intensities <- all_intensities[[1]]
                    } else {
                      # Extract raw intensities
                      intensities <- Cardinal::intensity(sub_data)
                      
                      # If multiple features, combine according to display mode
                      if (nrow(sub_data) > 1) {
                        intensities <- apply(intensities, 2, mean, na.rm = TRUE)
                      } else {
                        intensities <- as.vector(intensities)
                      }
                    }
                    all_intensities[[i]] <- intensities
                  }
                  
                  # Combine intensities based on display mode or superposition
                  if (input$superpose && length(all_intensities) > 1) {
                    # Average intensities for superposed ions
                    pixel_values <- Reduce("+", all_intensities) / length(all_intensities)
                  } else if (length(all_intensities) == 1) {
                    pixel_values <- all_intensities[[1]]
                  } else {
                    pixel_values <- all_intensities[[1]] # Default to first ion if issues
                  }
                  
                  # Remove NA and Inf values
                  pixel_values[is.na(pixel_values) | is.infinite(pixel_values)] <- 0
                  
                  # Apply contrast enhancement if selected
                  if (input$contrast3 == "histogram") {
                    # Histogram equalization for better contrast
                    n_breaks <- 100
                    pixel_values <- rank(pixel_values) / length(pixel_values) * max(pixel_values, na.rm = TRUE)
                  } else if (input$contrast3 == "adaptive") {
                    # Simple adaptive contrast enhancement
                    q_min <- quantile(pixel_values, 0.05, na.rm = TRUE)
                    q_max <- quantile(pixel_values, 0.95, na.rm = TRUE)
                    pixel_values <- (pixel_values - q_min) / (q_max - q_min) * max(pixel_values, na.rm = TRUE)
                    pixel_values[pixel_values < 0] <- 0
                  }
                  
                  # Normalize values to color range
                  if (any(pixel_values > 0)) {
                    # Generate colormap with appropriate palette
                    colormap <- cpal(input$color3)
                    
                    # Map intensities to colors, ensuring proper scaling
                    nonzero_vals <- pixel_values[pixel_values > 0]
                    if (length(nonzero_vals) > 0) {
                      # Filter very low values (likely noise)
                      intensity_threshold <- quantile(nonzero_vals, input$intensity_threshold, na.rm = TRUE)
                      pixel_values[pixel_values < intensity_threshold] <- 0
                      
                      # Scale values to colormap range
                      scaled_values <- pixel_values
                      if (max(scaled_values, na.rm = TRUE) > 0) {
                        scaled_values <- (scaled_values - min(nonzero_vals, na.rm = TRUE)) / 
                                        (max(scaled_values, na.rm = TRUE) - min(nonzero_vals, na.rm = TRUE))
                        scaled_values[scaled_values < 0] <- 0
                        scaled_values[scaled_values > 1] <- 1
                      }
                      
                      # Map to color indices
                      color_indices <- round(scaled_values * (length(colormap) - 1)) + 1
                      color_indices[color_indices < 1 | is.na(color_indices)] <- 1
                      image_colors <- colormap[color_indices]
                      
                      # Set transparent color for zero values
                      image_colors[pixel_values == 0] <- NA
                      
                      # Render clean MSI image
                      plot.new()
                      plot.window(xlim = range(coords$x, na.rm = TRUE),
                                ylim = range(coords$y, na.rm = TRUE))
                      
                      # Plot each pixel with appropriate color and size
                      points(coords$x, coords$y, 
                            pch = 15, 
                            col = image_colors,
                            cex = 1.5)  # Adjust pixel size as needed
                      
                      # Add axes and title
                      box()
                      if (length(mz_set) == 1) {
                        title(main = paste("m/z", round(mz_set, 4)))
                      } else if (length(mz_set) > 1) {
                        title(main = paste("Multiple m/z values"))
                      }
                      
                      if (input$colorkey3) {
                        # Add color key with better legend formatting
                        n_colors <- min(10, length(unique(image_colors)))
                        legend_colors <- colormap[seq(1, length(colormap), length.out = n_colors)]
                        legend("right", 
                              legend = rep("", length(legend_colors)), 
                              fill = legend_colors, 
                              border = NA,
                              bty = "n")
                      }
                      
                      cat("Clean MSI rendering completed for custom ion\n")
                    } else {
                      cat("No non-zero intensity values found for selected m/z\n")
                      # Fall back to standard rendering if no valid intensities
                      print(do.call(Cardinal::image, arg_list))
                    }
                  } else {
                    cat("No intensity values above zero for selected m/z\n")
                    # Fall back to standard rendering if no valid intensities
                    print(do.call(Cardinal::image, arg_list))
                  }
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
              if (!is.null(allInputs$render_mode) && allInputs$render_mode == "clean") {
                # Create clean MSI image with direct data extraction
                tryCatch({
                  # Get the first m/z value
                  selected_mz <- mz(overview_peaks_sel)[1]
                  cat("Using first m/z value:", selected_mz, "\n")
                  
                  # Extract grid coordinates directly
                  coords <- coord(overview_peaks_sel)
                  
                  # Apply preprocessing for cleaner visualization
                  msset_clean <- overview_peaks_sel
                  if (input$smooth3 != "none") {
                    msset_clean <- tryCatch({
                      Cardinal::smoothSignal(msset_clean, method = input$smooth3)
                    }, error = function(e) {
                      cat("Smoothing failed, using original data:", e$message, "\n")
                      msset_clean
                    })
                  }
                  
                  # Extract intensity values for the selected m/z
                  tol_value <- if(!is.null(input$plusminus_viz3)) input$plusminus_viz3 else 0.05
                  
                  # Subset features based on m/z tolerance
                  sub_data <- Cardinal::subsetFeatures(msset_clean, 
                                                      mz >= selected_mz - tol_value & 
                                                      mz <= selected_mz + tol_value)
                  
                  # Extract intensities
                  intensities <- Cardinal::intensity(sub_data)
                  
                  # If multiple features match, combine them
                  if (nrow(sub_data) > 1) {
                    pixel_values <- apply(intensities, 2, mean, na.rm = TRUE)
                  } else {
                    pixel_values <- as.vector(intensities)
                  }
                  
                  # Remove NA and Inf values
                  pixel_values[is.na(pixel_values) | is.infinite(pixel_values)] <- 0
                  
                  # Apply contrast enhancement if selected
                  if (input$contrast3 == "histogram") {
                    # Histogram equalization for better contrast
                    pixel_values <- rank(pixel_values) / length(pixel_values) * max(pixel_values, na.rm = TRUE)
                  } else if (input$contrast3 == "adaptive") {
                    # Simple adaptive contrast enhancement
                    q_min <- quantile(pixel_values, 0.05, na.rm = TRUE)
                    q_max <- quantile(pixel_values, 0.95, na.rm = TRUE)
                    pixel_values <- (pixel_values - q_min) / (q_max - q_min) * max(pixel_values, na.rm = TRUE)
                    pixel_values[pixel_values < 0] <- 0
                  }
                  
                  # Filter low-intensity values (likely noise)
                  nonzero_vals <- pixel_values[pixel_values > 0]
                  if (length(nonzero_vals) > 0) {
                    intensity_threshold <- quantile(nonzero_vals, input$intensity_threshold, na.rm = TRUE)
                    pixel_values[pixel_values < intensity_threshold] <- 0
                    
                    # Generate colormap with selected palette
                    colormap <- cpal(input$color3)
                    
                    # Scale values to colormap range
                    scaled_values <- pixel_values
                    if (max(scaled_values, na.rm = TRUE) > 0) {
                      scaled_values <- (scaled_values - min(nonzero_vals, na.rm = TRUE)) / 
                                       (max(scaled_values, na.rm = TRUE) - min(nonzero_vals, na.rm = TRUE))
                      scaled_values[scaled_values < 0] <- 0
                      scaled_values[scaled_values > 1] <- 1
                    }
                    
                    # Map to color indices
                    color_indices <- round(scaled_values * (length(colormap) - 1)) + 1
                    color_indices[color_indices < 1 | is.na(color_indices)] <- 1
                    image_colors <- colormap[color_indices]
                    
                    # Set transparent color for zero values
                    image_colors[pixel_values == 0] <- NA
                    
                    # Render clean MSI image
                    plot.new()
                    plot.window(xlim = range(coords$x, na.rm = TRUE),
                               ylim = range(coords$y, na.rm = TRUE))
                    
                    # Plot each pixel with appropriate color and size
                    points(coords$x, coords$y, 
                          pch = 15, 
                          col = image_colors,
                          cex = 1.5)  # Adjust pixel size as needed
                    
                    # Add axes and title
                    box()
                    title(main = paste("m/z", round(selected_mz, 4)))
                    
                    if (input$colorkey3) {
                      # Add color key with better legend formatting
                      n_colors <- min(10, length(unique(na.omit(image_colors))))
                      legend_colors <- colormap[seq(1, length(colormap), length.out = n_colors)]
                      legend("right", 
                             legend = rep("", length(legend_colors)), 
                             fill = legend_colors, 
                             border = NA,
                             bty = "n")
                    }
                    
                    cat("Clean MSI rendering completed\n")
                  } else {
                    cat("No non-zero intensity values found for selected m/z\n")
                    # Fall back to standard rendering
                    image_command <- Cardinal::image(overview_peaks_sel, 
                                      col=cpal(input$color3),
                                      scale=input$normalize3,
                                      key=(input$colorkey3))
                    print(image_command)
                  }
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
              
              cat("Applying histology overlay with alpha =", allInputs$alpha, "\n")
              
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
              alpha_value <- if (!is.null(allInputs$alpha) && is.numeric(allInputs$alpha)) {
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
              histology_grob <- editGrob(
                histology_grob,
                vp = viewport(
                  x = unit(0.5, "npc") + unit(allInputs$translate_x, "mm"),
                  y = unit(0.5, "npc") + unit(allInputs$translate_y, "mm"),
                  angle = allInputs$rotate,
                  width = unit(width_value, "npc"),  # Normalized coordinates with scaling
                  height = unit(height_value, "npc"), # Normalized coordinates with scaling
                  just = c("center", "center")
                )
              )
              
              # Draw the overlay with more detailed logging
              cat("Drawing histology overlay with transform parameters:\n")
              cat(" - Translation: X=", allInputs$translate_x, "mm, Y=", allInputs$translate_y, "mm\n")
              cat(" - Rotation:", allInputs$rotate, "degrees\n")
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
             
             vp_orig<-vizi_par()
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
                  width = input$width,
                  height = input$height,
                  alt = "This is alternate text")
           }, deleteFile = TRUE)
           
         })
       }
  })
        
  
}

