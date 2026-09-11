
#.libPaths(c("/datastore/lbcfs/labs/proteomics/Rserver/QC_Fig_Generator/lib", .libPaths()))


library(tidyverse)
library(shiny)
library(DT)
library(markdown)


# Functions

xcalibur = function(ID, method, numSamples, techReps, numPools, StartingVialPosition, Inj_Vol,
                    conditions = NULL,   # NEW ARG
                    randomize = c("None", "All", "Block"), 
                    comment = NULL){
  
  randomize <- match.arg(randomize)
  
  # split method
  
  sep = strsplit(method, "-")[[1]]
  
  instrument = sep[1]
  
  method2 = sep[2]
  
  # -----Select sample method -----
  if(method == "Astral-DIA"){
    sample_method = "C:\\Xcalibur\\methods\\Neo\\DIA_3mz_IO25cm_30min"
  }
  if(method == "Astral-DDA"){
    sample_method = "C:\\Xcalibur\\methods\\Neo\\DDA\\DDA_0,6s_MS1-300AGC-5ms_IO25cm_30min"
  }
  if(method == "OE480-DDA"){
    sample_method = "C:\\Xcalibur\\methods\\DDA Global\\PC-DDA_2s_90min"
  }
  
  # ----- Select blank method -----
  
  if(instrument == "Astral"){
    blank_method = "C:\\Xcalibur\\methods\\Neo\\FullScan_IO25cm_15min"
  }
  if(instrument == "OE480"){
    blank_method = "C:\\Xcalibur\\methods\\Blanks&Standards\\Blank-DDA_UDP_flush_NB"
  }
  
  # ----- Column selection -----
  if(instrument == "Astral"){
    column = "IO-25cm"
  }
  if(instrument == "OE480"){
    column = "IO-15cm"
  }
  
  date = format(Sys.Date(), "%Y%m%d")
  
  # ----- Normalize conditions (treat empty as NULL) -----
  if (is.null(conditions) || length(conditions) == 0 || all(is.na(conditions))) {
    conditions <- NULL
  }
  
  # ----- Validate conditions ONLY if Block -----
  if (randomize == "Block") {
    if (is.null(conditions)) {
      stop("Block randomization requires conditions.")
    }
    if (length(conditions) != numSamples) {
      stop("Length of conditions must equal numSamples")
    }
  } else {
    # Ignore conditions for other modes
    conditions <- NULL
  }  
  
  # ----- Base names -----
  base_name = paste(date, instrument, column, ID, method2, sep = "_")
  
  samp_ids = seq_len(numSamples)
  
  samp_base_name = paste(base_name, samp_ids, sep = "_")
  
  # ----- Technical replicates -----
  
  if(techReps > 0){
    samp_base_name = paste0(
      rep(samp_base_name, each = techReps),
      "_R", rep(1:techReps, times = length(samp_base_name))
    )
    
    if(!is.null(conditions)){
      conditions = rep(conditions, each = techReps)
    }
  }

  
  # ----- Create working dataframe -----
  df_samples = data.frame(
    name = samp_base_name,
    condition = if(is.null(conditions)) NA else conditions,
    stringsAsFactors = FALSE
  )
  
  # ----- Randomization -----
  
  # Case 1: All
  if(randomize == "All"){
    df_samples = df_samples[sample(nrow(df_samples)), ]
  }
  
  # Case 2: Block
  if(randomize == "Block"){
    
    if(is.null(conditions)){
      stop("Block randomization requires conditions.")
    }
    
    # Split by condition
    split_cond <- split(df_samples, df_samples$condition)
    
    # Shuffle each condition independently
    split_cond <- lapply(split_cond, function(df){
      df[sample(nrow(df)), , drop = FALSE]
    })
    
    blocks <- list()
    
    repeat {
      
      # Check if any condition still has samples left
      remaining <- sapply(split_cond, nrow)
      if(all(remaining == 0)) break
      
      block_rows <- list()
      
      # Take ONE sample from EACH condition (if available)
      for(cond in names(split_cond)){
        if(nrow(split_cond[[cond]]) > 0){
          
          # Randomly pick ONE row from remaining pool
          idx <- sample(seq_len(nrow(split_cond[[cond]])), 1)
          
          block_rows[[cond]] <- split_cond[[cond]][idx, , drop = FALSE]
          
          # Remove that row (no replacement)
          split_cond[[cond]] <- split_cond[[cond]][-idx, , drop = FALSE]
        }
      }
      
      block <- do.call(rbind, block_rows)
      
      # Randomize order within block
      block <- block[sample(nrow(block)), , drop = FALSE]
      
      blocks[[length(blocks) + 1]] <- block
    }
    
    df_samples <- do.call(rbind, blocks)
  }
  
  samp_base_name = df_samples$name
  
  # ----- Split blocks (post-randomization) -----
  
  if(randomize == "Block"){
    splt_base_names = lapply(blocks, function(b) b$name)
  } else {
    splt_base_names = list(samp_base_name)  # treat as one big block
  }
  
  
  out = c()
  
  # ----- Pools -----
  
  pool_base_names <- character(0)
  
  if (numPools >= 1) {
    pool_name <- paste(date, instrument, column, ID, method2, "Pool", sep = "_")
    pool_base_names <- paste0(pool_name, "_R", seq_len(numPools))
  }
  
  # ----- Blank generator -----
  
  blank_counter <- 1
  
  make_blank <- function() {
    nm <- paste(date, instrument, column, "blank", blank_counter, sep = "_")
    blank_counter <<- blank_counter + 1
    nm
  }
  
  # ----- Track block membership for samples -----
  
  sample_block_id <- rep(seq_along(splt_base_names), lengths(splt_base_names))
  
  append_samples_with_block_blanks <- function(out, start_idx, end_idx) {
    
    if (start_idx > end_idx) {
      return(out)
    }
    
    idxs <- start_idx:end_idx
    
    for (j in seq_along(idxs)) {
      
      idx <- idxs[j]
      
      # Add sample
      out <- c(out, all_samples[idx])
      
      # If this sample and the next sample are in different blocks,
      # insert a blank between blocks
      if (randomize == "Block" && idx < length(all_samples)) {
        if (sample_block_id[idx] != sample_block_id[idx + 1]) {
          out <- c(out, make_blank())
        }
      }
    }
    
    return(out)
  }
  
  
  # ----- Build output sequence with correctly spaced pools -----
  
  all_samples <- unlist(splt_base_names, use.names = FALSE)
  n_samples <- length(all_samples)
  
  
  # Always start with two blanks
  out <- c(out, make_blank(), make_blank())
  
  # --------- Case 1: No pools ---------
  if (numPools == 0) {
    
    out <- append_samples_with_block_blanks(out, 1, n_samples)
    out <- c(out, make_blank(), make_blank())
  }
  
  # --------- Case 2: One pool ---------
  if (numPools == 1) {
    
    out <- c(out,
             pool_base_names[1],
             make_blank(),
             all_samples,
             make_blank(),
             make_blank())
  }
  
  # --------- Case 3: Two or more pools ---------
  if (numPools >= 2) {
    
    # First pool always goes at the beginning
    out <- c(out,
             pool_base_names[1],
             make_blank())
    
    # If exactly 2 pools:
    # Pool 1 at beginning, Pool 2 at end
    if (numPools == 2) {
      
      out <- append_samples_with_block_blanks(out, 1, n_samples)
      
      out <- c(out,
               make_blank(),
               pool_base_names[2],
               make_blank(),
               make_blank())
      
    }
    # If 3 or more pools:
    # Pool 1 at beginning, last pool at end,
    # middle pools evenly spaced within samples
    if (numPools >= 3) {
      
      middle_pools <- pool_base_names[2:(numPools - 1)]
      last_pool <- pool_base_names[numPools]
      
      n_middle <- length(middle_pools)
      
      # Safety check:
      # Need at least enough sample positions to place middle pools between sample chunks
      if (n_middle > n_samples - 1) {
        stop("Too many pools requested for the number of samples. Reduce numPools or increase numSamples.")
      }
      
      # Evenly spaced cut points.
      # Example: 20 samples, 4 pools -> middle pools after sample 7 and sample 14
      cut_points <- ceiling(seq_along(middle_pools) * n_samples / (numPools - 1))
      
      current_start <- 1
      
      for (i in seq_along(middle_pools)) {
        
        current_end <- cut_points[i]
        
        # Add samples up to the cut point
        if (current_start <= current_end) {
          out <- append_samples_with_block_blanks(out, current_start, current_end)
        }
        
        # Add blank, pool, blank
        out <- c(out,
                 make_blank(),
                 middle_pools[i],
                 make_blank())
        
        current_start <- current_end + 1
      }
      
      # Add remaining samples after the last middle pool
      if (current_start <= n_samples) {
        out <- append_samples_with_block_blanks(out, current_start, n_samples)
      }
      
      # Final pool at the end
      out <- c(out,
               make_blank(),
               last_pool,
               make_blank(),
               make_blank())
    }
  }

  
    # ----- Vial position assignment (NEW LOGIC) -----
  
  assign_positions <- function(instrument, out, sample_names, pool_names, StartingVialPosition) {
    
    # Format helper
    format_pos <- function(plate, row, col) {
      if (instrument == "Astral") {
        return(paste0(plate, ":", row, col))
      } else {
        return(paste0(plate, row, col))
      }
    }
    
    # Fixed blank position
    blank_pos <- if (instrument == "Astral") "R:A1" else "RA1"
    
    # Plate structure
    if (instrument == "Astral") {
      plates <- c("R", "Y", "B", "G")
      rows   <- LETTERS[1:6]  # A-F
      cols   <- 1:9
    } else {
      plates <- c("R", "G", "B")
      rows   <- LETTERS[1:5]  # A-E
      cols   <- 1:8
    }
    
    # Generate ordered sample positions
    sample_positions <- c()
    
    for (p in plates) {
      
      plate_rows <- if (p == "R") {
        rows[-1]     # remove A from red only
      } else {
        rows
      }
      
      for (r in plate_rows) {
        for (c in cols) {
          sample_positions <- c(sample_positions,
                                format_pos(p, r, c))
        }
      }
    }
  
    start_idx <- match(StartingVialPosition, sample_positions)
    
    if (is.na(start_idx)) {
      stop(
        paste(
          "StartingVialPosition",
          StartingVialPosition,
          "is not a valid position"
        )
      )
    }
    
    # Reorder positions so assignment starts there
    sample_positions <- c(
      sample_positions[start_idx:length(sample_positions)],
      sample_positions[1:(start_idx - 1)]
    )
    
    
    
    # Assign positions by category
    positions <- character(length(out))
    
    sample_idx <- 1
    pool_idx   <- 1
    
    # Total counts
    n_samples <- length(sample_names)
    
    for (i in seq_along(out)) {
      
      current <- out[i]
      
      if (grepl("blank", current)) {
        positions[i] <- blank_pos
        
      } else if (current %in% sample_names) {
        
        # Extract original sample number (before _R replicate if present)
        sample_num <- as.numeric(sub(".*_([0-9]+)(_R[0-9]+)?$", "\\1", current))
        
        positions[i] <- sample_positions[sample_num]
        
        
      } else if (current %in% pool_names) {
        # Pools come AFTER all samples
        positions[i] <- sample_positions[n_samples + 1]
        
      } else {
        positions[i] <- NA
      }
    }
    
    return(positions)
  }
    
  
  # Identify sample + pool names
  sample_names <- samp_base_name
  
  pool_names <- if (numPools >= 2) pool_base_names else character(0)
  
  positions <- assign_positions(instrument, out, sample_names, pool_names, StartingVialPosition)
  
  # ----- Output -----
  path = paste0("D:\\data", format(Sys.Date(), "%Y_%m_%B"))
  
  # Create method assignment vector
  inst_methods <- ifelse(
    grepl("blank", out, ignore.case = TRUE),
    blank_method,
    sample_method
  )
  
  df = tibble(
    `File Name` = out,
    Path = path,
    `Inst Meth` = inst_methods,
    Position = positions,
    `Inj Vol` = rep(Inj_Vol, length(out)),
    Comment = comment
  )
  
  return(df)
}


ui <- navbarPage("Xcalibur Sequence Generator",
                 tabPanel("Sample Information",
                          sidebarPanel(
                            h3("Input Sample Information"),
                            textInput("uniqueID", label = "Input Unique Identifier:", value = ""),
                            radioButtons("method", label = "Select Method:", 
                                         choices = c("Astral-DIA", "Astral-DDA", "OE480-DDA"), selected = "Astral-DIA"),
                            numericInput("nSamples", label = "Number of Samples (not including pools)", value = 2, min = 2),
                            numericInput("nTechReps", label = "Number of technical injections (0 = none):", value = 0, min = 2),
                            numericInput("nPools", label = "Number of Pools:", value = 0, min = 2),
                            textInput("StartingVialPosition", label = "Starting Vial Position: (Color-position)", value = "R:B1"),
                            numericInput("Inj_Volume", label = "Injection volume for all samples:", value = 1),
                            textInput("comment", label = "Add comment for all samples:", value = ""),
                            radioButtons("Randomize", label = "Select a type of sample randomization:", choices = c("None", "All", "Block"), selected = "None"),
                            conditionalPanel(
                              condition = "input.Randomize == 'Block'",
                              textAreaInput(
                                "conditions",
                                "Paste sample conditions (one per line, must match sample order):",
                                placeholder = "Control\nTreatment\nTreatment\nControl",
                                rows = 10
                              )
                            ),
                            actionButton("GenerateSequence", label = "Generate Sequence", class = "btn btn-success"),
                            textInput("FileName", label = "File Name", value = ""),
                            downloadButton("DownloadSequence", label = "Download Sequence", icon = icon("download"), class = "btn btn-primary")
                          ),
                          mainPanel(
                            DT::dataTableOutput("df")
                          )))

server <- function(input, output, session){
  
  sequence_generator = eventReactive(input$GenerateSequence,{
    
    seq = xcalibur(ID = input$uniqueID,
                   method = input$method,
                   numSamples = input$nSamples,
                   techReps = input$nTechReps,
                   numPools = input$nPools,
                   StartingVialPosition = input$StartingVialPosition,
                   Inj_Vol = input$Inj_Volume,
                   conditions = trimws(strsplit(input$conditions, "\n")[[1]]),
                   comment = input$comment,
                   randomize = input$Randomize)
  })
  
  
  observeEvent(input$method, {
    
    if (grepl("^Astral", input$method)) {
      
      updateSelectInput(
        session,
        "StartingVialPosition",
        selected = "R:B1"
      )
      
    } else {
      
      updateSelectInput(
        session,
        "StartingVialPosition",
        selected = "RB1"
      )
      
    }
  })
  
  DT = reactive({
    seqTable = sequence_generator()
  })
  
  # display seq
  output$df = DT::renderDataTable(DT())
  
  # download sequence
  output$DownloadSequence = downloadHandler(
    filename = function(){
      paste0(format(Sys.time(), '%Y%m%d'), "_", input$FileName, "_sequence.csv")
    },
    content = function(file){
      
      df <- sequence_generator()
      
      con <- file(file, open = "wt")
      
      # Row 1
      writeLines("Bracket Type=4", con)
      
      # Row 2 onward (column headers + data)
      write.table(
        df,
        file = con,
        sep = ",",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = "",
        append = TRUE
      )
      
      close(con)
    }
  )
  
}

# Run the app
shinyApp(ui, server)