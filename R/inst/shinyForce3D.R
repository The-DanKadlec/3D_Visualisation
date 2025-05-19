# Load required libraries
library(animation)
library(tidyr)
library(ggplot2)
library(shinyjs)
library(dplyr)
library(rgl)
library(shinyRGL)   
library(readr)
library(grid)
library(magick)
library(webshot2)
library(plot3D)
library(shinyWidgets)
library(stringr)

# Define UI for the app
ui <- fluidPage(
  useShinyjs(),
  tags$div(
    style = "text-align: left;",  
    tags$h1("3D Visualisation of Force Vector"),  
    tags$h4("Kadlec & Vial et al., in preperation", style = "font-size: 14px; color: grey;")
  ),
  sidebarLayout(
    sidebarPanel(
      div(style = "display: flex; align-items: center; gap: 10px;",
          div(style = "flex-grow: 1;",
              fileInput(
                "file1",
                label = tagList(
                  tags$strong("Dataset A"), 
                  tags$span(
                    style  = "font-weight: normal; margin-left: 8px; color: #555;",
                    uiOutput("trialInfo1", inline = TRUE)
                  )
                ),
                accept = ".csv",
                width  = "100%"
              ))
      ),
      
      div(
        style = "display: flex; align-items: center; gap: 10px;",
        div(style = "flex-grow: 1;",
            fileInput(
              "file2",
              label = tagList(
                tags$strong("Dataset B"),
                tags$span(
                  style = "font-weight: normal; margin-left: 8px; color: #555;",
                  uiOutput("trialInfo2", inline = TRUE)
                )
              ),
              accept = ".csv",
              width = "100%"
            )
        ),
        uiOutput("removeFile2UI")
      ),
      
      tags$hr(),
      tags$h4(strong("3D Plot Options")),  # <-- Add this line
      checkboxInput("showLines1", "Show Butterfly for First Dataset", FALSE),
      checkboxInput("showLines2", "Show Butterfly for Second Dataset", FALSE),
      checkboxInput("showPeakResultantLine", "Show vertical line at peak Fxyz", value = FALSE),
      sliderInput(
        "radiusMult",
        "Sphere‑size multiplier:",
        min   = 0.005,
        max   = 0.05,
        value = 0.015,
        step  = 0.005
      ),
      checkboxInput("showStaticShadows", "Show static‐3D shadows", TRUE),
      tags$h4(strong("2D Plot Options")),
      radioButtons(
        inputId = "plotType2D",
        label   = "2D Plot Mode:",
        choices = c("Mean ± SD", "All Trials"),
        selected = "Mean ± SD",
        inline  = TRUE
      ),
      fluidRow(
        column(
          6,
          checkboxGroupInput(
            "displayComponents",
            "Display curves:",
            choices  = c("Fx", "Fy", "Fz", "Resultant"),
            selected = c("Fx", "Fy", "Fz", "Resultant")
          )
        ),
        column(
          6,
          conditionalPanel(
            condition = "input.plotType2D == 'Mean ± SD'",
            checkboxGroupInput(
              "showPeaks",
              "Show peak of:",
              choices  = c("Fx", "Fy", "Fz", "Resultant"),
              selected = NULL
            )
          )
        )  
      ),
      tags$h4(strong("Pedotti Plot Options")),
      selectInput("pedottiComp1", "X-Axis Component:", choices = c("Fx", "Fy", "Fz"), selected = "Fy"),
      selectInput("pedottiComp2", "Y-Axis Component:", choices = c("Fx", "Fy", "Fz"), selected = "Fz"),
    ),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Dynamic 3D Plot",  
                 tagList(
                   tags$h4("3D Force Vector", style = "margin-top: 10px; font-weight: bold;"),
                   rglwidgetOutput("plot3D", width = "800px", height = "800px"),
                   fluidRow(
                     column(6, plotOutput("colorLegendDynamic", height = "65px", width = "400px")),
                     column(6, plotOutput("colorLegendDynamic2", height = "65px", width = "400px"))
                   ),
                   tags$hr(style = "border-top: 1px solid #999; margin: 20px 0;"), 
                   fluidRow(
                     column(6, tags$h4("3D Force Vector GIF")),
                     column(3, actionButton("createGIF", "Create GIF")),
                     column(3, downloadButton("download3DGIF", "Download GIF")),
                   ),
                   uiOutput("gifDisplay"))
        ),
        tabPanel("Static 3D Plot", 
                 tagList(
                   tags$h4("3D Force Vector", style = "margin-top: 10px; font-weight: bold;"),
                   plotOutput("staticRow", height = "600px", width = "800px"),
                   fluidRow(
                     column(6, plotOutput("colorLegend", height = "65px", width = "400px")),
                     column(6, plotOutput("colorLegend2", height = "65px", width = "400px"))
                   )
                 )
        ),
        tabPanel("2D Plot",
                 tagList(
                   tags$h4("Force Vector (Fx, Fy, Fz, Resultant)", style = "margin-top: 10px; font-weight: bold;"),
                   plotOutput("combinedPlot", height = "550px", width = "800px"),
                   tags$hr(style = "border-top: 1px solid #999; margin: 20px 0;"),
                   # Fluid row for the polar header, toggle, and download buttons
                   fluidRow(
                     column(
                       width = 6,
                       # Wrap header + toggle together
                       tagList(
                         tags$h4("Force Vector (r, θ, φ)", style = "margin-top: 10px; font-weight: bold;"),
                         checkboxInput("polarInDegrees", "Display angles in degrees", FALSE)
                       )
                     ),
                     column(
                       width = 6,
                       # right‐align buttons if you like
                       div(style = "display: flex; justify-content: flex-end; gap: 10px;",
                           downloadButton("downloadPolarA", "Download Dataset A"),
                           downloadButton("downloadPolarB", "Download Dataset B")
                       )
                     )
                   ),
                   plotOutput("polarPlot", height = "750px", width = "800px")
                 )
        ),
        tabPanel("Pedotti Plot",
                 fluidRow(
                   column(12, tags$h4("Pedotti Plot", style = "font-weight: bold;"),
                          plotOutput("pedottiPlot", height = "800px", width = "800px"))
                 ),
                 fluidRow(
                   column(6, plotOutput("colorLegendPedottiA", height = "65px", width = "400px")),
                   column(6, plotOutput("colorLegendPedottiB", height = "65px", width = "400px"))
                 ),
                 tags$hr(style = "border-top: 1px solid #999; margin: 20px 0;"),
                 fluidRow(
                   column(6, tags$h4("Animated Pedotti Plot", style = "font-weight: bold;")),
                   column(3, actionButton("animatePedotti", "Animate Pedotti Plot")),
                   column(3, downloadButton("downloadPedottiGIF", "Download GIF"))
                 ),
                 fluidRow(
                   column(12, plotOutput("pedottiAnimationPlot", height = "800px", width = "800px"))
                 )
        ),
      ) # closes tabsetPanel
    )   # closes mainPanel
  ))

# Define server logic
server <- function(input, output, session) {
  
  pedottiFrame <- reactiveVal(1)
  pedottiAnimating <- reactiveVal(FALSE)
  # Reactive file uploads
  file1_loaded <- reactiveVal(NULL)
  file2_loaded <- reactiveVal(NULL)
  
  haveData <- reactive({ !is.null(file1_loaded()) || !is.null(file2_loaded()) })
  
  observeEvent(input$file1, { file1_loaded(input$file1) })
  observeEvent(input$file2, { file2_loaded(input$file2) })
  observeEvent(input$removeFile2, { file2_loaded(NULL); reset("file2") })
  
  dataInput1 <- reactive({
    req(file1_loaded())
    df <- read_csv(file1_loaded()$datapath)
    # assume columns in exact order Fx Fy Fz | Fx Fy Fz | …
    ntr <- ncol(df) / 3
    if (ntr %% 1 != 0) {
      stop("File1 must have a multiple of 3 columns (Fx, Fy, Fz each trial).")
    }
    ntr <- as.integer(ntr)
    # compute trial list (if you ever want them individually):
    trials <- lapply(seq_len(ntr), function(i) {
      slice <- df[, ((i-1)*3 + 1):(i*3)]
      colnames(slice) <- c("Fx","Fy","Fz")
      slice
    })
    # compute mean across trials
    mean_df <- data.frame(
      Fx = rowMeans(sapply(trials, `[[`, "Fx")),
      Fy = rowMeans(sapply(trials, `[[`, "Fy")),
      Fz = rowMeans(sapply(trials, `[[`, "Fz"))
    )
    list(ntrials = ntr, mean = mean_df, trials = trials)
  })
  
  output$trialInfo1 <- renderUI({
    inp <- dataInput1()
    txt <- if (inp$ntrials == 1) {
      "1 trial"
    } else {
      paste0(inp$ntrials, " trials → plotting mean")
    }
    tags$div(
      strong("Trial Count:"), 
      txt
    )
  })
  
  dataInput2 <- reactive({
    req(file2_loaded())
    df <- read_csv(file2_loaded()$datapath)
    # assume columns in exact order Fx Fy Fz | Fx Fy Fz | …
    ntr <- ncol(df) / 3
    if (ntr %% 1 != 0) {
      stop("File1 must have a multiple of 3 columns (Fx, Fy, Fz each trial).")
    }
    ntr <- as.integer(ntr)
    # compute trial list (if you ever want them individually):
    trials <- lapply(seq_len(ntr), function(i) {
      slice <- df[, ((i-1)*3 + 1):(i*3)]
      colnames(slice) <- c("Fx","Fy","Fz")
      slice
    })
    # compute mean across trials
    mean_df <- data.frame(
      Fx = rowMeans(sapply(trials, `[[`, "Fx")),
      Fy = rowMeans(sapply(trials, `[[`, "Fy")),
      Fz = rowMeans(sapply(trials, `[[`, "Fz"))
    )
    list(ntrials = ntr, mean = mean_df, trials = trials)
  })
  
  output$trialInfo2 <- renderUI({
    req(file2_loaded())
    inp <- dataInput2()
    txt <- if (inp$ntrials == 1) {
      "1 trial"
    } else {
      paste0(inp$ntrials, " trials → plotting mean")
    }
    tags$div(
      strong("Trial Count:"), 
      txt
    )
  })
  
  output$removeFile2UI <- renderUI({
    if (is.null(input$file2)) return(NULL)
    actionButton(
      "removeFile2",
      label = "❌",
      style = "color: red; margin-top: 5px;",
      title = "Remove Dataset B"
    )
  })
  
  # Helper: compute resultants, axis limits, etc.
  params <- reactive({
    inp1 <- if (!is.null(file1_loaded())) dataInput1() else NULL
    inp2 <- if (!is.null(file2_loaded())) dataInput2() else NULL
    
    df1 <- if (!is.null(inp1)) inp1$mean else NULL
    df2 <- if (!is.null(inp2)) inp2$mean else NULL
    
    res1 <- if (!is.null(df1)) sqrt(df1$Fx^2 + df1$Fy^2 + df1$Fz^2) else numeric(0)
    res2 <- if (!is.null(df2)) sqrt(df2$Fx^2 + df2$Fy^2 + df2$Fz^2) else numeric(0)
    
    global_max <- max(c(res1, res2), na.rm=TRUE)
    if(global_max == 0) global_max <- 1
    
    all_z  <- c(if(!is.null(df1)) df1$Fz, if(!is.null(df2)) df2$Fz)
    peak_z <- if(length(all_z)) max(abs(all_z), na.rm=TRUE) else 1
    
    list(
      df1        = df1,
      df2        = df2,
      res1       = res1,
      res2       = res2,
      global_max = global_max,
      axis_len   = peak_z * 1.1,
      half_ax    = (peak_z * 1.1) / 2,
      stance_n   = if (!is.null(df1)) nrow(df1) else if (!is.null(df2)) nrow(df2) else 101
    )
  })
  
  output$download3DGIF <- downloadHandler(
    filename = function() {
      paste0("3D_animation_", Sys.Date(), ".gif")
    },
    content = function(file) {
      gif_file <- file.path("www", "3D_animation.gif")
      req(file.exists(gif_file))
      file.copy(gif_file, file)
    }
  )
  
  # Helper function to convert Cartesian to Polar for a list of trials
  convert_to_polar <- function(trials) {
    trial_names <- paste0("Trial_", seq_along(trials))
    results <- lapply(seq_along(trials), function(i) {
      trial <- trials[[i]]
      with(trial, {
        r <- sqrt(Fx^2 + Fy^2 + Fz^2)
        θ <- acos(pmin(pmax(Fz / r, -1), 1))
        φ <- atan2(Fy, Fx)
        data.frame(r, θ, φ)
      })
    })
    
    # Combine and transpose: each trial as columns, rows = time points
    r_mat <- do.call(cbind, lapply(results, function(x) x$r))
    θ_mat <- do.call(cbind, lapply(results, function(x) x$θ))
    φ_mat <- do.call(cbind, lapply(results, function(x) x$φ))
    
    colnames(r_mat) <- paste0("Trial", seq_along(results), "_r")
    colnames(θ_mat) <- paste0("Trial", seq_along(results), "_theta")
    colnames(φ_mat) <- paste0("Trial", seq_along(results), "_phi")
    
    out <- cbind(r_mat, θ_mat, φ_mat)
    out <- as.data.frame(out)
    out <- cbind(Stance = seq_len(nrow(out)), out)
    return(out)
  }
  
  make_polar_download <- function(trials, file, in_degrees) {
    # 1) compute polar
    pd <- convert_to_polar(trials)
    
    # 2) convert angles if needed
    unit_suffix <- if (in_degrees) "deg" else "rad"
    if (in_degrees) {
      pd <- pd %>%
        mutate(across(matches("_(theta|phi)$"), ~ . * 180 / pi))
    }
    
    # 3) rename columns: append suffix to r, theta, phi
    names(pd) <- names(pd) %>%
      str_replace("_r$",     paste0("_r_",     unit_suffix)) %>%
      str_replace("_theta$", paste0("_theta_", unit_suffix)) %>%
      str_replace("_phi$",   paste0("_phi_",   unit_suffix))
    
    # 4) reorder: Stance first, then for each trial its r/theta/phi trio
    #    find all trial prefixes
    prefixes <- names(pd) %>%
      grep("^Trial\\d+_r_", ., value = TRUE) %>%
      str_remove("_r_.*$")
    # build the desired order
    ordered_cols <- c(
      "Stance",
      unlist(lapply(prefixes, function(tr) {
        paste0(tr, c("_r_", "_theta_", "_phi_"), unit_suffix)
      }))
    )
    pd <- pd[, ordered_cols, drop = FALSE]
    
    # 5) write
    write.csv(pd, file, row.names = FALSE)
  }
  
  output$downloadPolarA <- downloadHandler(
    filename = function() {
      suf <- if (input$polarInDegrees) "deg" else "rad"
      paste0("polar_coords_dataset_A_", suf, "_", Sys.Date(), ".csv")
    },
    contentType = "text/csv",
    content = function(file) {
      req(file1_loaded())
      make_polar_download(dataInput1()$trials, file, input$polarInDegrees)
    }
  )
  
  output$downloadPolarB <- downloadHandler(
    filename = function() {
      suf <- if (input$polarInDegrees) "deg" else "rad"
      paste0("polar_coords_dataset_B_", suf, "_", Sys.Date(), ".csv")
    },
    contentType = "text/csv",
    content = function(file) {
      req(file2_loaded())
      make_polar_download(dataInput2()$trials, file, input$polarInDegrees)
    }
  )
  
  # Start animation loop
  observeEvent(input$animatePedotti, {
    pedottiFrame(0)
    pedottiAnimating(TRUE)
  })
  
  observe({
    invalidateLater(50, session)
    isolate({
      if (pedottiAnimating()) {
        i <- pedottiFrame()
        if (i < params()$stance_n) {
          pedottiFrame(i + 1)
        } else {
          pedottiFrame(0)  # Loop forever
        }
      }
    })
  })
  
  pedottiFrame <- reactiveVal(0)
  output$pedottiAnimationPlot <- renderPlot({
    req(haveData())
    comp1 <- input$pedottiComp1
    comp2 <- input$pedottiComp2
    validate(need(comp1 != comp2, "Please select two different components."))
    
    p <- params()
    step <- pedottiFrame()
    df_list <- list()
    
    if (!is.null(p$df1)) {
      df1 <- p$df1 %>%
        mutate(Resultant = sqrt(Fx^2 + Fy^2 + Fz^2),
               Dataset = "A",
               Time = row_number(),
               Color = colorRampPalette(c("orange", "red", "purple"))(n())) %>%
        select(Fx, Fy, Fz, Resultant, Dataset, Time, Color) %>%
        filter(Time > 0, Time <= step)
      df_list[[1]] <- df1
    }
    
    if (!is.null(p$df2)) {
      df2 <- p$df2 %>%
        mutate(Resultant = sqrt(Fx^2 + Fy^2 + Fz^2),
               Dataset = "B",
               Time = row_number(),
               Color = colorRampPalette(c("cadetblue1", "cyan3", "blue"))(n())) %>%
        select(Fx, Fy, Fz, Resultant, Dataset, Time, Color) %>%
        filter(Time > 0, Time <= step)
      df_list[[2]] <- df2
    }
    
    full_data <- bind_rows(df_list)
    max_res <- params()$global_max
    arrow_data <- full_data %>%
      select(X = all_of(comp1), Y = all_of(comp2), Color)
    
    ggplot() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
      geom_segment(
        data = arrow_data,
        aes(x = 0, y = 0, xend = X, yend = Y, color = Color),
        arrow = arrow(length = unit(0.08, "inches")),
        alpha = 0.85
      ) +
      coord_fixed(xlim = c(-(max_res / 6), max_res), ylim = c(-(max_res / 6), max_res)) +
      scale_color_identity(guide = "none") +
      labs(x = comp1, y = comp2) +
      theme_minimal(base_size = 16)
  })
  
  ################################################################################
  
  output$pedottiPlot <- renderPlot({
    req(haveData())
    comp1 <- input$pedottiComp1
    comp2 <- input$pedottiComp2
    validate(need(comp1 != comp2, "Please select two different components."))
    
    output$colorLegendPedottiA <- renderPlot({
      req(file1_loaded())
      df1 <- params()$df1
      n <- nrow(df1)
      cols <- colorRampPalette(c("orange", "red", "purple"))(n)
      par(mar = c(2, 1, 2, 1))  # small margins
      image(1:n, 1, as.matrix(1:n), col = cols, axes = FALSE, useRaster = TRUE)
      title("0–100% Stance (A)", line = 1, cex.main = 1.1)
    })
    
    output$colorLegendPedottiB <- renderPlot({
      req(file2_loaded())
      df2 <- params()$df2
      n2 <- nrow(df2)
      cols2 <- colorRampPalette(c("cadetblue1", "cyan3", "blue"))(n2)
      par(mar = c(2, 1, 2, 1))  # small margins
      image(1:n2, 1, as.matrix(1:n2), col = cols2, axes = FALSE, useRaster = TRUE)
      title("0–100% Stance (B)", line = 1, cex.main = 1.1)
    })
    
    p <- params()
    df_list <- list()
    
    if (!is.null(p$df1)) {
      df1 <- p$df1 %>%
        mutate(
          Resultant = sqrt(Fx^2 + Fy^2 + Fz^2),
          Dataset = "A",
          Time = row_number(),
          Color = colorRampPalette(c("orange", "red", "purple"))(n())
        ) %>%
        select(Fx, Fy, Fz, Resultant, Dataset, Time, Color)
      df_list[[1]] <- df1
    }
    
    if (!is.null(p$df2)) {
      df2 <- p$df2 %>%
        mutate(
          Resultant = sqrt(Fx^2 + Fy^2 + Fz^2),
          Dataset = "B",
          Time = row_number(),
          Color = colorRampPalette(c("cadetblue1", "cyan3", "blue"))(n())
        ) %>%
        select(Fx, Fy, Fz, Resultant, Dataset, Time, Color)
      df_list[[2]] <- df2
    }
    
    full_data <- bind_rows(df_list)
    max_res <- max(full_data$Resultant, na.rm = TRUE)
    
    arrow_data <- full_data %>%
      select(X = all_of(comp1), Y = all_of(comp2), Color, Time, Dataset)
    
    ggplot() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
      geom_segment(
        data = arrow_data,
        aes(x = 0, y = 0, xend = X, yend = Y, color = Color),
        arrow = arrow(length = unit(0.08, "inches")),
        alpha = 0.85
      ) +
      coord_fixed(xlim = c(-(max_res / 6), max_res), ylim = c(-(max_res / 6), max_res)) +
      scale_color_identity(guide = "none") +
      labs(x = comp1, y = comp2) +
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 16)
      )
  })
  
  output$downloadPedottiGIF <- downloadHandler(
    filename = function() {
      paste0("pedotti_diagram_", Sys.Date(), ".gif")
    },
    content = function(file) {
      req(haveData())
      comp1 <- isolate(input$pedottiComp1)
      comp2 <- isolate(input$pedottiComp2)
      validate(need(comp1 != comp2, "Cannot generate GIF with identical axes."))
      
      p <- params()
      n_frames <- p$stance_n
      df1 <- p$df1
      df2 <- p$df2
      
      # Create temp directory for frames
      tmpdir <- tempdir()
      frame_paths <- file.path(tmpdir, sprintf("frame%03d.png", 1:n_frames))
      
      for (i in seq_len(n_frames)) {
        step <- i
        
        df_list <- list()
        if (!is.null(df1)) {
          df1_step <- df1 %>%
            mutate(Resultant = sqrt(Fx^2 + Fy^2 + Fz^2),
                   Dataset = "A",
                   Time = row_number(),
                   Color = colorRampPalette(c("orange", "red", "purple"))(n())) %>%
            select(Fx, Fy, Fz, Resultant, Dataset, Time, Color) %>%
            filter(Time <= step)
          df_list[[1]] <- df1_step
        }
        if (!is.null(df2)) {
          df2_step <- df2 %>%
            mutate(Resultant = sqrt(Fx^2 + Fy^2 + Fz^2),
                   Dataset = "B",
                   Time = row_number(),
                   Color = colorRampPalette(c("cadetblue1", "cyan3", "blue"))(n())) %>%
            select(Fx, Fy, Fz, Resultant, Dataset, Time, Color) %>%
            filter(Time <= step)
          df_list[[2]] <- df2_step
        }
        
        full_data <- bind_rows(df_list)
        max_res <- p$global_max
        
        arrow_data <- full_data %>%
          select(X = all_of(comp1), Y = all_of(comp2), Color)
        
        p_frame <- ggplot() +
          geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
          geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
          geom_segment(
            data = arrow_data,
            aes(x = 0, y = 0, xend = X, yend = Y, color = Color),
            arrow = arrow(length = unit(0.08, "inches")),
            alpha = 0.85
          ) +
          coord_fixed(xlim = c(-(max_res / 6), max_res), ylim = c(-(max_res / 6), max_res)) +
          scale_color_identity(guide = "none") +
          labs(x = comp1, y = comp2) +
          theme_minimal(base_size = 16) +                # <-- minimal background
          theme(
            panel.background = element_rect(fill = "white", color = NA),  # <-- white background
            plot.background = element_rect(fill = "white", color = NA),
            #axis.line = element_line(color = "black"),
            axis.title = element_text(color = "black")
          )
        
        ggsave(filename = frame_paths[i], plot = p_frame, width = 8, height = 8, units = "in", dpi = 100)
      }
      
      # Combine frames into a GIF
      imgs <- magick::image_read(frame_paths)
      gif <- magick::image_animate(imgs, fps = 10)
      magick::image_write(gif, path = file)
      
      # Clean up
      file.remove(frame_paths)
    }
  )
  ################################################################################
  # Dynamic 3D Plot
  output$plot3D <- renderRglwidget({
    req(haveData())
    p   <- params()
    df1 <- p$df1; df2 <- p$df2
    
    #max_radius <- p$axis_len * 0.015
    max_radius <- p$axis_len * input$radiusMult
    max_res1   <- if(length(p$res1)) max(p$res1, na.rm=TRUE) else 0
    max_res2   <- if(length(p$res2)) max(p$res2, na.rm=TRUE) else 0
    radius_A   <- (max_res1 / p$global_max) * max_radius
    radius_B   <- (max_res2 / p$global_max) * max_radius
    
    open3d(useNULL=TRUE)
    aspect3d(1, 1, 1)
    plot3d(NA,
           xlim = c(-p$half_ax, p$half_ax),
           ylim = c(-p$half_ax, p$half_ax),
           zlim = c(0, p$axis_len),
           xlab = "Fx", ylab = "Fy", zlab = "Fz",
           type = "n")
    
    if(!is.null(df1)) {
      cols1 <- colorRampPalette(c("orange","red","purple"))(nrow(df1))
      spheres3d(df1$Fx, df1$Fy, df1$Fz, col=cols1, radius=radius_A)
      if(input$showLines1) {
        for(i in seq_len(nrow(df1))) {
          lines3d(c(0, df1$Fx[i]), c(0, df1$Fy[i]), c(0, df1$Fz[i]), col=cols1[i], lwd=1)
        }
      }
      if(input$showPeakResultantLine && max_res1>0) {
        i1 <- which.max(p$res1)
        lines3d(c(0, df1$Fx[i1]), c(0, df1$Fy[i1]), c(0, df1$Fz[i1]), col="black", lwd=4)
      }
    }
    
    if(!is.null(df2)) {
      cols2 <- colorRampPalette(c("cadetblue1","cyan3","blue"))(nrow(df2))
      spheres3d(df2$Fx, df2$Fy, df2$Fz, col=cols2, radius=radius_B)
      if(input$showLines2) {
        for(i in seq_len(nrow(df2))) {
          lines3d(c(0, df2$Fx[i]), c(0, df2$Fy[i]), c(0, df2$Fz[i]), col=cols2[i], lwd=1)
        }
      }
      if(input$showPeakResultantLine && max_res2>0) {
        i2 <- which.max(p$res2)
        lines3d(c(0, df2$Fx[i2]), c(0, df2$Fy[i2]), c(0, df2$Fz[i2]), col="black", lwd=4)
      }
    }
    
    rglwidget()
  })
  
  ################################################################################
  # Static 3D Plot with true 3-plane shadows
  output$staticRow <- renderPlot({
    req(haveData())
    p   <- params()
    df1 <- p$df1; df2 <- p$df2
    
    # cube limits
    xlim   <- c(-p$half_ax, p$half_ax)
    ylim   <- c(-p$half_ax, p$half_ax)
    zlim   <- c(0,         p$axis_len)
    floor_z <- zlim[1]   # = 0
    
    # precompute the “floor”, “back” and “side” constants for each view:
    back_y_0  <- -p$half_ax   # for θ=0, y = -half_ax
    side_x_0  <- -p$half_ax   # for θ=0, x = -half_ax
    
    back_y_90 <- -p$half_ax   # for θ=90, back‐wall is still the “far” edge in y?
    side_x_90 <-  p$half_ax   # for θ=90, side‐wall is the +x face
    
    # size scaling (as before)
    base_cex <- c(0.5,2.0)
    normA    <- if(length(p$res1)) max(p$res1)/p$global_max else 1
    normB    <- if(length(p$res2)) max(p$res2)/p$global_max else 1
    sr       <- input$radiusMult / 0.015
    cexA     <- base_cex[1] + diff(base_cex)*normA*sr
    cexB     <- base_cex[1] + diff(base_cex)*normB*sr
    
    par(mfrow=c(1,2), mar=c(2,2,1,1), mgp=c(1.5,0.5,0), pty="s")
    
    for(theta in c(0, 90)) {
      par(cex.axis=1.5, cex.lab=1.5)
      shadow_col <- "#99999933"
      
      # draw dataset A
      if (!is.null(df1)) {
        cols1 <- colorRampPalette(c("orange","red","purple"))(nrow(df1))
        scatter3D(df1$Fx,df1$Fy,df1$Fz,
                  theta=theta, phi=15,
                  xlim=xlim, ylim=ylim, zlim=zlim,
                  pch=16, cex=cexA, col=cols1, colvar=NULL, colkey=FALSE,
                  xlab="Fx", ylab="Fy", zlab="Fz", bty="b", ticktype="detailed")
        if(input$showLines1){
          segments3D(
            x0 = rep(0, nrow(df1)), y0 = rep(0, nrow(df1)), z0 = rep(0, nrow(df1)),
            x1 = df1$Fx,           y1 = df1$Fy,           z1 = df1$Fz,
            col = cols1, add = TRUE
          )
        }
        if(input$showPeakResultantLine){
          i1 <- which.max(p$res1)
          segments3D(0,0,0,
                     df1$Fx[i1], df1$Fy[i1], df1$Fz[i1],
                     col="black", lwd=3, add=TRUE)
        }
        
        # **now the shadows**, split by θ
        if (theta == 0) {
          if (input$showStaticShadows) {
            # LEFT panel shadows
            points3D(df1$Fx, df1$Fy, rep(floor_z, nrow(df1)),
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(df1$Fx, rep(-back_y_0,  nrow(df1)), df1$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(rep(-side_x_0, nrow(df1)), df1$Fy, df1$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
          }
        } else {
          if (input$showStaticShadows) {
            # RIGHT panel shadows
            points3D(df1$Fx, df1$Fy, rep(floor_z, nrow(df1)),
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(df1$Fx, rep(-back_y_90,  nrow(df1)), df1$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(rep(-side_x_90, nrow(df1)), df1$Fy, df1$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
          }}
      }
      
      # draw dataset B (same split)
      if (!is.null(df2)) {
        cols2 <- colorRampPalette(c("cadetblue1","cyan3","blue"))(nrow(df2))
        points3D(df2$Fx, df2$Fy, df2$Fz,
                 pch=16, cex=cexB, col=cols2, colvar=NULL, colkey=FALSE, add=TRUE)
        if(input$showLines2){
          segments3D(
            x0 = rep(0, nrow(df2)), y0 = rep(0, nrow(df2)), z0 = rep(0, nrow(df2)),
            x1 = df2$Fx,           y1 = df2$Fy,           z1 = df2$Fz,
            col = cols2, add = TRUE
          )
        }
        if(input$showPeakResultantLine){
          i2 <- which.max(p$res2)
          segments3D(0,0,0,
                     df2$Fx[i2], df2$Fy[i2], df2$Fz[i2],
                     col="black", lwd=3, add=TRUE)
        }
        
        if (theta == 0) {
          if (input$showStaticShadows) {
            points3D(df2$Fx, df2$Fy, rep(floor_z, nrow(df2)),
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(df2$Fx, rep(-back_y_0,  nrow(df2)), df2$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(rep(-side_x_0, nrow(df2)), df2$Fy, df2$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
          }
        } else {
          if (input$showStaticShadows) {
            # RIGHT
            points3D(df2$Fx, df2$Fy, rep(floor_z, nrow(df2)),
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(df2$Fx, rep(-back_y_90,  nrow(df2)), df2$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
            points3D(rep(-side_x_90, nrow(df2)), df2$Fy, df2$Fz,
                     col=shadow_col, pch=16, cex=0.6, add=TRUE)
          }}
      }
    }
  })
  
  ################################################################################
  # Combined 2D Plot
  output$combinedPlot <- renderPlot({
    req(haveData())
    
    trial_list_to_df <- function(trials, dataset_label) {
      bind_rows(lapply(seq_along(trials), function(i) {
        trials[[i]] %>%
          mutate(Index = row_number(),
                 Trial = paste0("Trial_", i),
                 Dataset = dataset_label)
      }))
    }
    
    data_list <- list()
    if (!is.null(file1_loaded())) {
      d1 <- dataInput1()
      data_list[[1]] <- trial_list_to_df(d1$trials, "A")
    }
    if (!is.null(file2_loaded())) {
      d2 <- dataInput2()
      data_list[[2]] <- trial_list_to_df(d2$trials, "B")
    }
    
    full_data <- bind_rows(data_list) %>%
      mutate(Resultant = sqrt(Fx^2 + Fy^2 + Fz^2))
    
    max_y <- max(full_data$Resultant, na.rm = TRUE) * 1.10
    
    long_data <- full_data %>%
      pivot_longer(cols = c(Fx, Fy, Fz, Resultant),
                   names_to = "Component",
                   values_to = "Value") %>%
      filter(Component %in% input$displayComponents) %>%
      mutate(
        Component = factor(Component, levels = c("Fx", "Fy", "Fz", "Resultant")),
        Trial = factor(Trial)
      )
    
    if (input$plotType2D == "All Trials") {
      p <- ggplot(long_data, aes(Index, Value,
                                 color = Component,
                                 linetype = Dataset,
                                 group = interaction(Dataset, Trial, Component))) +
        geom_line(alpha = 0.7)
    } else {
      summary_data <- long_data %>%
        group_by(Dataset, Component, Index) %>%
        summarise(
          mean = mean(Value, na.rm = TRUE),
          sd   = sd(Value, na.rm = TRUE),
          .groups = "drop"
        )
      
      p <- ggplot(summary_data, aes(x = Index, y = mean, color = Component, linetype = Dataset)) +
        geom_line(size = 1.2) +
        geom_ribbon(
          data = summary_data,
          aes(x = Index, ymin = mean - sd, ymax = mean + sd,
              group = interaction(Component, Dataset), fill = Component),
          alpha = 0.2,
          inherit.aes = FALSE
        ) +
        scale_color_manual(values = c(Fx = "blue", Fy = "purple", Fz = "red", Resultant = "black"), drop = FALSE) +
        scale_fill_manual(values = c(Fx = "blue", Fy = "purple", Fz = "red", Resultant = "black"), drop = FALSE)
    }
    
    p <- p +
      scale_color_manual(values = c(Fx = "blue", Fy = "purple", Fz = "red", Resultant = "black"), drop = FALSE) +
      scale_linetype_manual(values = c("A" = "solid", "B" = "dashed"), drop = FALSE) +
      geom_hline(yintercept = 0, size = 0.5) +
      labs(title = NULL, x = "Stance (%)", y = "GRF (N)") +
      theme_classic(base_size = 14) +
      theme(
        legend.position     = "bottom",
        legend.box          = "vertical",
        legend.title        = element_blank(),
        legend.text         = element_text(size = 16),
        legend.key.width    = unit(1.2, "cm"),
        legend.key.height   = unit(1, "lines"),
        strip.placement     = "outside",
        strip.text.y.left   = element_text(angle = 0, size = 14, face = "bold"),
        panel.spacing.y     = unit(0.5, "cm")
      ) +
      ylim(-1.10 * max(abs(full_data$Fy), na.rm = TRUE),
           1.10 * max(full_data$Resultant, na.rm = TRUE))
    
    if (input$plotType2D == "Mean ± SD" && length(input$showPeaks)) {
      peak_data <- long_data %>%
        filter(Component %in% input$showPeaks) %>%
        group_by(Dataset, Component, Index) %>%
        summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop") %>%
        group_by(Dataset, Component) %>%
        slice_max(order_by = abs(mean_val), n = 1) %>%
        rename(Value = mean_val)
      
      p <- p + geom_segment(data = peak_data,
                            aes(x = Index, xend = Index, y = 0, yend = Value,
                                color = Component, linetype = Dataset),
                            linetype = "dotted", size = 0.8, show.legend = FALSE)
    }
    
    p
  })
  
  ######################POLAR
  output$polarPlot <- renderPlot({
    req(haveData())
    deg <- input$polarInDegrees
    
    # Combine all trials from both datasets
    data_list <- list()
    if (!is.null(file1_loaded())) {
      d1 <- dataInput1()
      data_list[[1]] <- bind_rows(lapply(seq_along(d1$trials), function(i) {
        d1$trials[[i]] %>%
          mutate(Index = row_number(), Trial = paste0("Trial_", i), Dataset = "A")
      }))
    }
    if (!is.null(file2_loaded())) {
      d2 <- dataInput2()
      data_list[[2]] <- bind_rows(lapply(seq_along(d2$trials), function(i) {
        d2$trials[[i]] %>%
          mutate(Index = row_number(), Trial = paste0("Trial_", i), Dataset = "B")
      }))
    }
    full_data <- bind_rows(data_list)
    
    # Convert to spherical coordinates
    polar_data <- full_data %>%
      mutate(
        r   = sqrt(Fx^2 + Fy^2 + Fz^2),
        θ   = acos(pmin(pmax(Fz / r, -1), 1)),  # Clamp for stability
        φ   = atan2(Fy, Fx),
        Dataset = factor(Dataset)
      ) %>%
      group_by(Dataset, Index) %>%
      summarise(
        r = mean(r, na.rm = TRUE),
        θ = mean(θ, na.rm = TRUE),
        φ = mean(φ, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Define the labeller first
    facet_labeller <- as_labeller(c(
      r  = "Magnitude (N)",
      θ  = "Inclination (rad)",
      φ  = "Azimuth (rad)"
    ))
    
    # Convert to degrees if requested
    if (deg) {
      polar_data <- polar_data %>%
        mutate(
          θ = θ * 180 / pi,
          φ = φ * 180 / pi
        )
    }
    
    # 3) now pivot once
    polar_long <- polar_data %>%
      pivot_longer(
        cols = c(r, θ, φ),
        names_to  = "Component",
        values_to = "Value"
      ) %>%
      mutate(Component = factor(Component, levels = c("r", "θ", "φ")))
    
    # Axis labels depending on units
    deg_label <- if (deg) " (°)" else " (rad)"
    
    # Create individual plots
    plot_r <- ggplot(filter(polar_long, Component == "r"), aes(Index, Value, color = Dataset, linetype = Dataset)) +
      geom_line(size = 1.2) +
      labs(y = "Magnitude (N)", x = NULL) +
      scale_color_manual(values = c("A" = "black", "B" = "green")) +
      scale_linetype_manual(values = c("A" = "solid", "B" = "dashed")) +
      theme_classic(base_size = 14) +
      theme(
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    plot_theta <- ggplot(filter(polar_long, Component == "θ"), aes(Index, Value, color = Dataset, linetype = Dataset)) +
      geom_line(size = 1.2) +
      labs(y = paste0("Inclination", deg_label), x = NULL) +
      scale_color_manual(values = c("A" = "black", "B" = "green")) +
      scale_linetype_manual(values = c("A" = "solid", "B" = "dashed")) +
      theme_classic(base_size = 14) +
      theme(
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none"
      )
    
    plot_phi <- ggplot(filter(polar_long, Component == "φ"), aes(Index, Value, color = Dataset, linetype = Dataset)) +
      geom_line(size = 1.2) +
      labs(y = paste0("Azimuth", deg_label), x = "Stance (%)") +
      scale_color_manual(values = c("A" = "black", "B" = "green")) +
      scale_linetype_manual(values = c("A" = "solid", "B" = "dashed")) +
      theme_classic(base_size = 14) +
      theme(
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 16), 
        legend.key.width    = unit(1.2, "cm")
      )
    
    # Combine using patchwork
    library(patchwork)
    combined_plot <- plot_r / plot_theta / plot_phi + plot_layout(heights = c(1, 1, 1))
    
    # Render in output$polarPlot or wherever you want
    print(combined_plot)
  })
  
  ################################################################################
  # Color legends
  output$colorLegendDynamic <- renderPlot({
    req(file1_loaded())
    df1 <- params()$df1
    n <- nrow(df1)
    cols <- colorRampPalette(c("orange", "red", "purple"))(n)
    
    # Use layout to control spacing instead of manual par reset
    par(mar = c(2, 1, 2, 1))  # bottom, left, top, right
    image(1:n, 1, as.matrix(1:n), col = cols, axes = FALSE, useRaster = TRUE)
    mtext("0–100% Stance (A)", side = 3, line = 0.5, cex = 1.2)
  })
  
  output$colorLegendDynamic2 <- renderPlot({
    req(file2_loaded())
    df2 <- params()$df2
    n2 <- nrow(df2)
    cols2 <- colorRampPalette(c("cadetblue1", "cyan3", "blue"))(n2)
    
    par(mar = c(2, 1, 2, 1))
    image(1:n2, 1, as.matrix(1:n2), col = cols2, axes = FALSE, useRaster = TRUE)
    mtext("0–100% Stance (B)", side = 3, line = 0.5, cex = 1.2)
  })
  
  output$colorLegend <- renderPlot({
    req(file1_loaded())
    df1 <- params()$df1
    n <- nrow(df1)
    cols <- colorRampPalette(c("orange", "red", "purple"))(n)
    
    # Use layout to control spacing instead of manual par reset
    par(mar = c(2, 1, 2, 1))  # bottom, left, top, right
    image(1:n, 1, as.matrix(1:n), col = cols, axes = FALSE, useRaster = TRUE)
    mtext("0–100% Stance (A)", side = 3, line = 0.5, cex = 1.2)
  })
  
  output$colorLegend2 <- renderPlot({
    req(file2_loaded())
    df2 <- params()$df2
    n2 <- nrow(df2)
    cols2 <- colorRampPalette(c("cadetblue1", "cyan3", "blue"))(n2)
    
    par(mar = c(2, 1, 2, 1))
    image(1:n2, 1, as.matrix(1:n2), col = cols2, axes = FALSE, useRaster = TRUE)
    mtext("0–100% Stance (B)", side = 3, line = 0.5, cex = 1.2)
  })
  
  ################################################################################
  # GIF creator (inside server, only once!)
  observeEvent(input$createGIF, {
    req(haveData())
    p   <- params()
    df1 <- p$df1; df2 <- p$df2
    
    # a) bubble radii
    #max_radius <- p$axis_len * 0.015
    max_radius <- p$axis_len * input$radiusMult
    max_res1   <- if(length(p$res1)) max(p$res1, na.rm=TRUE) else 0
    max_res2   <- if(length(p$res2)) max(p$res2, na.rm=TRUE) else 0
    radius_A   <- max_radius
    radius_B   <- max_radius
    
    # b) folders
    if(!dir.exists("www")) dir.create("www")
    if(!dir.exists("frames")) dir.create("frames")
    
    # c) headless scene
    # 2. headless RGL + **plot your data directly** (just like your older code)
    open3d(windowRect=c(0,0,1000,1000))
    bg3d("white")   
    aspect3d(1, 1, 1)  # equal scaling
    par3d(userMatrix = diag(4))  # reset any previous rotation
    
    # plot Dataset A as coloured spheres
    cols1 <- colorRampPalette(c("orange","red","purple"))(nrow(df1))
    plot3d(df1$Fx, df1$Fz, df1$Fy,
           col   = cols1,
           radius= radius_A,
           type  = "s",
           xlim  = c(-p$half_ax,p$half_ax),
           ylim  = c(0,p$axis_len),
           zlim  = c(-p$half_ax,p$half_ax),
           xlab  = "Fx", ylab = "Fy", zlab = "Fz")
    
    
    # 5) add butterfly & peak for A
    if (input$showLines1) {
      for (i in seq_len(nrow(df1))) {
        lines3d(c(0,df1$Fx[i]), c(0,df1$Fz[i]), c(0,df1$Fy[i]), col=cols1[i], lwd=1)
      }
    }
    if (input$showPeakResultantLine && max_res1>0) {
      i1 <- which.max(p$res1)
      lines3d(c(0,df1$Fx[i1]), c(0,df1$Fz[i1]), c(0,df1$Fy[i1]), col="black", lwd=4)
    }
    
    # now add Dataset B as spheres, same logic
    if (!is.null(df2)) {
      cols2 <- colorRampPalette(c("cadetblue1","cyan3","blue"))(nrow(df2))
      spheres3d(df2$Fx, df2$Fz, df2$Fy,
                col    = cols2,
                radius = radius_B)
      
      # 7) add butterfly & peak for B
      if (input$showLines2) {
        for (i in seq_len(nrow(df2))) {
          lines3d(c(0,df2$Fx[i]), c(0,df2$Fz[i]), c(0,df2$Fy[i]), col=cols2[i], lwd=1)
        }
      }
      if (input$showPeakResultantLine && max_res2>0) {
        i2 <- which.max(p$res2)
        lines3d(c(0,df2$Fx[i2]), c(0,df2$Fz[i2]), c(0,df2$Fy[i2]), col="black", lwd=4)
      }
    }
    
    # e) spin & capture
    view3d(theta =   0, phi = 15, zoom = 0.75)  
    n_frames <- 100
    # first frame at theta = 0
    rgl.snapshot(sprintf("frames/frame%03d.png", 0))
    for(i in 1:(n_frames-1)) {
      view3d(theta = i*(360/n_frames), phi = 15, zoom = 0.75)
      rgl.snapshot(sprintf("frames/frame%03d.png", i))
    }
    
    
    close3d()
    
    # f) stitch
    imgs <- image_read(list.files("frames",pattern="png$",full.names=TRUE))
    gif  <- image_animate(imgs, fps=10)
    image_write(gif, path=file.path("www","3D_animation.gif"))
    unlink("frames", recursive=TRUE)
    
    showNotification("GIF created!", type="message")
    output$gifDisplay <- renderUI({
      gif_path <- "3D_animation.gif"
      if (file.exists(file.path("www", gif_path))) {
        tags$img(src = paste0(gif_path, "?", as.numeric(Sys.time())), width = "800px")
      } else {
        tags$p("No GIF created yet.")
      }
    })
  })
  
}

# Ensure 'www' directory exists
if (!dir.exists("www")) { dir.create("www") }

shinyApp(ui, server)


#Packages used:
#Jackson, J. (2017). animation: A Package for Animations in R. R package version 2.10.0
#Wickham, H. (2020). tidyr: Tidy Messy Data. R package version 1.3.5
#Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer. R package version 3.4.4
#Attali, D. (2016). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 2.1.6
#Wickham, H. et al. (2023). dplyr: A Grammar of Data Manipulation. R package version 1.1.4
#Murrell, P. (2024). rgl: 3D Visualization Using OpenGL. R package version 0.111.2
#Peterson, B. (2015). shinyRGL: Interactive 3D Scatterplots in Shiny. R package version 0.94
#Wickham, H. (2023). readr: Read Rectangular Text Data. R package version 2.2.1
#R Core Team (2024). grid: The Grid Graphics Package. 
#Ooms, J. & the magick team (2024). magick: Advanced Graphics and Image-processing in R. R package version 2.9.3
#Ooms, J. (2024). webshot2: Take Screenshots of Web Pages. R package version 1.0.0
#Soetaert, K. (2013). plot3D: Plotting Multi-Dimensional Data. R package version 1.5.2
#Chang, W. (2023). shinyWidgets: Custom Inputs Widgets for Shiny. R package version 0.8.0
#Pedersen, T. L. (2023). patchwork: The Composer of Plots. R package version 1.2.1