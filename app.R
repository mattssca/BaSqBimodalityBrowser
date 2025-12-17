# Load required packages
required_packages <- c("shiny", "ggplot2", "patchwork", "dplyr", "DT")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Load pre-computed data and functions
load("data/app_data.RData")  # Contains all data including binomial_results
source("functions/bimodality_functions.R")

# Validate data loaded correctly
cat("Data loaded successfully!\n")
cat("Samples:", ncol(predicted_imvigor$data), "\n")
cat("Genes:", nrow(predicted_imvigor$data), "\n")
cat("BaSq samples:", nrow(basq_proliferations), "\n") 
cat("Bimodality results:", nrow(binomial_results), "genes\n")
cat("Proliferation genes:", length(proliferation_genes), "\n")
cat("Infiltration genes:", length(infiltration_genes), "\n")

# UI
ui <- fluidPage(
  titlePanel("BaSq Bimodality Browser"),
  
  tabsetPanel(
    # Tab 1: Gene Explorer
    tabPanel("Gene Explorer",
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 
                 # Gene input
                 h4("Gene Selection"),
                 textInput("gene_input", 
                           "Enter gene name:",
                           value = "RFC4",
                           placeholder = "e.g., MKI67, CD8A, TP53"),
                 
                 helpText("Enter a gene symbol (case sensitive)"),
                 
                 # Validation message
                 uiOutput("gene_validation"),
                 
                 # Gene summary
                 h4("Gene Summary"),
                 verbatimTextOutput("gene_summary"),
                 
                 # Download button
                 br(),
                 downloadButton("download_plot", "Download Plot", class = "btn-primary")
               ),
               
               mainPanel(
                 width = 9,
                 
                 # Plot output
                 conditionalPanel(
                   condition = "output.gene_exists",
                   plotOutput("combined_plot", height = "500px")
                 ),
                 
                 # Error message if gene not found
                 conditionalPanel(
                   condition = "!output.gene_exists",
                   div(
                     h3("Gene not found", style = "color: red;"),
                     p("Please check the gene name and try again."),
                     p("Make sure the gene name is spelled correctly and is present in the dataset."),
                     style = "text-align: center; margin-top: 100px;"
                   )
                 ),
                 
                 # Gene statistics table
                 conditionalPanel(
                   condition = "output.gene_exists",
                   br(),
                   h4("Bimodality Analysis Results"),
                   tableOutput("gene_stats")
                 )
               )
             )
    ),
    
    # Tab 2: Bimodality Results Table  
    tabPanel("Bimodality Results",
             fluidRow(
               column(12,
                      h3("Pre-computed Bimodality Analysis Results"),
                      p("Browse all genes and their bimodality metrics. Click on a gene to visualize it in the Gene Explorer tab."),
                      
                      # Summary stats
                      div(
                        style = "background-color: #f8f9fa; padding: 15px; margin-bottom: 20px; border-radius: 5px;",
                        h4("Dataset Summary"),
                        fluidRow(
                          column(3, strong("Total Genes:"), textOutput("total_genes", inline = TRUE)),
                          column(3, strong("BaSq Samples:"), textOutput("total_samples", inline = TRUE)),
                          column(3, strong("Proliferation Genes:"), textOutput("prolif_count", inline = TRUE)),
                          column(3, strong("Infiltration Genes:"), textOutput("infiltr_count", inline = TRUE))
                        )
                      ),
                      
                      # Filter options
                      fluidRow(
                        column(3,
                               selectInput("filter_blacklist",
                                           "Gene Filter:",
                                           choices = list(
                                             "All genes" = "all",
                                             "Non-blacklisted only" = "exclude_blacklist",
                                             "Proliferation genes only" = "proliferation_only",
                                             "Infiltration genes only" = "infiltration_only"
                                           ),
                                           selected = "exclude_blacklist")
                        ),
                        column(3,
                               numericInput("min_bimodality",
                                            "Min Bimodality:",
                                            value = 0,
                                            min = 0,
                                            max = 1,
                                            step = 0.01)
                        ),
                        column(3,
                               numericInput("min_separation",
                                            "Min Separation:",
                                            value = 0,
                                            min = 0,
                                            max = 1,
                                            step = 0.01)
                        ),
                        column(3,
                               numericInput("table_rows",
                                            "Show rows:",
                                            value = 100,
                                            min = 10,
                                            max = 1000,
                                            step = 10)
                        )
                      ),
                      
                      hr(),
                      
                      # Data table
                      DT::dataTableOutput("bimodality_table")
               )
             )
    ),
    
    # Tab 3: About
    tabPanel("About",
             div(
               style = "padding: 20px;",
               h3("About This App"),
               p("This app explores bimodal gene expression patterns in Basal-Squamous (BaSq) bladder cancer samples."),
               
               h4("Data Source"),
               p("Expression data from IMvigor210 clinical trial BaSq samples, classified using LundTaxR."),
               
               h4("Analysis"),
               tags$ul(
                 tags$li("Samples were clustered based on proliferation scores using K-means (k=2)"),
                 tags$li("Bimodality analysis performed on ~15,000 genes"),
                 tags$li("Genes ranked by combined bimodality and cluster separation metrics")
               ),
               
               h4("Signatures"),
               tags$ul(
                 tags$li(strong("Proliferation genes:"), "Cell cycle and DNA replication related genes"),
                 tags$li(strong("Infiltration genes:"), "Immune cell marker genes"),
                 tags$li(strong("Blacklisted:"), "Genes from proliferation + infiltration signatures")
               ),
               
               h4("Usage"),
               p("Use the Gene Explorer to visualize individual genes, or browse the Bimodality Results table to discover interesting candidates.")
             )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Summary statistics for display
  output$total_genes <- renderText({
    nrow(binomial_results)
  })
  
  output$total_samples <- renderText({
    nrow(basq_proliferations)  
  })
  
  output$prolif_count <- renderText({
    length(proliferation_genes)
  })
  
  output$infiltr_count <- renderText({
    length(infiltration_genes)
  })
  
  # Check if gene exists in data
  gene_exists <- reactive({
    req(input$gene_input)
    toupper(trimws(input$gene_input)) %in% toupper(rownames(predicted_imvigor$data))
  })
  
  # Get the correctly formatted gene name
  formatted_gene <- reactive({
    req(input$gene_input)
    user_input <- toupper(trimws(input$gene_input))
    available_genes <- rownames(predicted_imvigor$data)
    
    # Find exact match (case insensitive)
    match_idx <- which(toupper(available_genes) == user_input)
    if (length(match_idx) > 0) {
      return(available_genes[match_idx[1]])
    } else {
      return(NULL)
    }
  })
  
  # Output gene existence for conditional panels
  output$gene_exists <- reactive({
    gene_exists()
  })
  outputOptions(output, "gene_exists", suspendWhenHidden = FALSE)
  
  # Gene validation message
  output$gene_validation <- renderUI({
    req(input$gene_input)
    
    if (gene_exists()) {
      div(icon("check", style = "color: green;"), 
          span("Gene found!", style = "color: green;"))
    } else {
      div(icon("times", style = "color: red;"), 
          span("Gene not found in dataset", style = "color: red;"))
    }
  })
  
  # Gene summary
  output$gene_summary <- renderText({
    req(formatted_gene())
    
    gene_name <- formatted_gene()
    
    # Check signature membership
    in_proliferation <- gene_name %in% proliferation_genes
    in_infiltration <- gene_name %in% infiltration_genes
    
    # Get expression statistics
    expr_raw <- predicted_imvigor$data[gene_name, basq_proliferations$sample_id]
    expr_vals <- as.numeric(expr_raw)
    expr_vals <- expr_vals[!is.na(expr_vals)]
    
    # Create summary text
    signature_status <- "None"
    if (in_proliferation && in_infiltration) {
      signature_status <- "Both proliferation AND infiltration"
    } else if (in_proliferation) {
      signature_status <- "Proliferation signature"
    } else if (in_infiltration) {
      signature_status <- "Infiltration signature"
    }
    
    if (length(expr_vals) == 0 || !is.numeric(expr_vals)) {
      return(paste0(
        "Gene: ", gene_name, "\n",
        "Signature: ", signature_status, "\n",
        "Blacklisted: ", ifelse(gene_name %in% blacklisted_genes, "Yes", "No"), "\n\n",
        "Expression Statistics: No valid data available"
      ))
    }
    
    paste0(
      "Gene: ", gene_name, "\n",
      "Signature: ", signature_status, "\n", 
      "Blacklisted: ", ifelse(gene_name %in% blacklisted_genes, "Yes", "No"), "\n\n",
      "Expression Statistics:\n",
      "  Range: ", round(min(expr_vals), 2), " - ", round(max(expr_vals), 2), "\n",
      "  Mean: ", round(mean(expr_vals), 2), "\n",
      "  Median: ", round(median(expr_vals), 2), "\n",
      "  Samples: ", length(expr_vals)
    )
  })
  
  # Generate the combined plot
  output$combined_plot <- renderPlot({
    req(formatted_gene(), gene_exists())
    
    tryCatch({
      plot_bimodal_gene(this_gene = formatted_gene(),
                        this_expr = predicted_imvigor$data,
                        this_cluster = basq_proliferations)
    }, error = function(e) {
      ggplot() + 
        annotate("text", x = 0, y = 0, 
                 label = paste("Error plotting", formatted_gene(), "\n", e$message), 
                 size = 6) +
        theme_void()
    })
  })
  
  # Gene statistics table
  output$gene_stats <- renderTable({
    req(formatted_gene(), gene_exists())
    
    gene_name <- formatted_gene()
    
    # Get statistics from pre-computed binomial_results
    if (gene_name %in% binomial_results$gene) {
      gene_stats <- binomial_results[binomial_results$gene == gene_name, ]
      
      data.frame(
        Metric = c("Overall Bimodality", "Cluster Separation Quality", 
                   "Valley Position", "Cluster 1 in High Mode (%)", 
                   "Cluster 2 in High Mode (%)", "Cluster Explains Bimodality", 
                   "Combined Score"),
        Value = c(
          round(gene_stats$overall_bimodality, 4),
          round(gene_stats$cluster_separation_quality, 4),
          round(gene_stats$valley_position, 4),
          paste0(round(gene_stats$cluster1_in_high_mode * 100, 1), "%"),
          paste0(round(gene_stats$cluster2_in_high_mode * 100, 1), "%"),
          round(gene_stats$cluster_explains_bimodality, 4),
          round(gene_stats$combined_score, 4)
        )
      )
    } else {
      data.frame(
        Metric = "Analysis Status",
        Value = "Gene not included in bimodality analysis"
      )
    }
  }, striped = TRUE, hover = TRUE)
  
  # Filtered bimodality results for table
  filtered_bimodality <- reactive({
    data <- binomial_results
    
    # Add signature information
    data$signature <- case_when(
      data$gene %in% proliferation_genes & data$gene %in% infiltration_genes ~ "Both",
      data$gene %in% proliferation_genes ~ "Proliferation",
      data$gene %in% infiltration_genes ~ "Infiltration", 
      TRUE ~ "None"
    )
    
    # Apply gene filter
    if (input$filter_blacklist == "exclude_blacklist") {
      data <- data %>% filter(!gene %in% blacklisted_genes)
    } else if (input$filter_blacklist == "proliferation_only") {
      data <- data %>% filter(gene %in% proliferation_genes)
    } else if (input$filter_blacklist == "infiltration_only") {
      data <- data %>% filter(gene %in% infiltration_genes)
    }
    
    # Apply numeric filters
    data <- data %>%
      filter(
        overall_bimodality >= input$min_bimodality,
        cluster_separation_quality >= input$min_separation
      )
    
    # Format for display
    data %>%
      select(
        Gene = gene,
        Signature = signature,
        `Bimodality Score` = overall_bimodality,
        `Cluster Separation` = cluster_separation_quality,
        `Valley Position` = valley_position,
        `Cluster 1 High %` = cluster1_in_high_mode,
        `Cluster 2 High %` = cluster2_in_high_mode,
        `Explains Bimodality` = cluster_explains_bimodality,
        `Combined Score` = combined_score
      ) %>%
      mutate(
        `Bimodality Score` = round(`Bimodality Score`, 4),
        `Cluster Separation` = round(`Cluster Separation`, 4),
        `Valley Position` = round(`Valley Position`, 4),
        `Cluster 1 High %` = paste0(round(`Cluster 1 High %` * 100, 1), "%"),
        `Cluster 2 High %` = paste0(round(`Cluster 2 High %` * 100, 1), "%"),
        `Explains Bimodality` = round(`Explains Bimodality`, 4),
        `Combined Score` = round(`Combined Score`, 4)
      ) %>%
      arrange(desc(as.numeric(`Combined Score`)))
  })
  
  # Bimodality results table
  output$bimodality_table <- DT::renderDataTable({
    data <- filtered_bimodality()
    
    DT::datatable(
      data,
      options = list(
        pageLength = input$table_rows,
        scrollX = TRUE,
        scrollY = "600px",
        dom = 'frtip'
      ),
      selection = 'single',
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        'Signature',
        backgroundColor = DT::styleEqual(
          c('Proliferation', 'Infiltration', 'Both', 'None'),
          c('#ffcccc', '#ccccff', '#ffccff', '#ffffff')
        )
      )
  })
  
  # Update gene input when table row is selected
  observeEvent(input$bimodality_table_rows_selected, {
    if (length(input$bimodality_table_rows_selected) > 0) {
      selected_row <- input$bimodality_table_rows_selected[1]
      selected_gene <- filtered_bimodality()$Gene[selected_row]
      updateTextInput(session, "gene_input", value = selected_gene)
      updateTabsetPanel(session, "tabset", selected = "Gene Explorer")
    }
  })
  
  # Download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      req(formatted_gene())
      paste0(formatted_gene(), "_bimodal_plot.png")
    },
    content = function(file) {
      req(formatted_gene(), gene_exists())
      plot <- plot_bimodal_gene(this_gene = formatted_gene(),
                                this_expr = predicted_imvigor$data,
                                this_cluster = basq_proliferations)
      ggsave(file, plot, width = 12, height = 6, dpi = 300)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)