# Load packages explicitly
library(shiny)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DT)

# Keep the installation logic as backup for the server
required_packages <- c("shiny", "ggplot2", "patchwork", "dplyr", "DT")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    options(repos = c(CRAN = "https://cran.rstudio.com/"))
    install.packages(pkg, dependencies = TRUE)
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
cat("Genes that passed filters:", sum(binomial_results$passes_filters), "\n")
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
                 
                 # Plot type selector
                 h4("Plot Options"),
                 radioButtons("plot_type", 
                              "Choose plot type:",
                              choices = list(
                                "Density plots" = "density",
                                "Ranked expression" = "ranked"
                              ),
                              selected = "density"),
                 
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
                          column(2, strong("Total Genes:"), textOutput("total_genes", inline = TRUE)),
                          column(2, strong("Passed Filters:"), textOutput("passed_genes", inline = TRUE)),
                          column(2, strong("BaSq Samples:"), textOutput("total_samples", inline = TRUE)),
                          column(3, strong("Proliferation Genes:"), textOutput("prolif_count", inline = TRUE)),
                          column(3, strong("Infiltration Genes:"), textOutput("infiltr_count", inline = TRUE))
                        )
                      ),
                      
                      # Filter options
                      fluidRow(
                        column(2,
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
                        column(2,
                               selectInput("filter_status",
                                           "Filter Status:",
                                           choices = list(
                                             "All genes" = "all",
                                             "Passed filters only" = "passed_only",
                                             "Failed filters only" = "failed_only"
                                           ),
                                           selected = "passed_only")
                        ),
                        column(2,
                               numericInput("min_bimodality",
                                            "Min Bimodality:",
                                            value = 0,
                                            min = 0,
                                            max = 5,
                                            step = 0.1)
                        ),
                        column(2,
                               numericInput("min_separation",
                                            "Min Separation:",
                                            value = 0,
                                            min = 0,
                                            max = 1,
                                            step = 0.01)
                        ),
                        column(2,
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
               
               h4("Analysis Overview"),
               tags$ul(
                 tags$li("Samples were clustered based on proliferation scores using K-means (k=2)"),
                 tags$li("Bimodality analysis performed on ~15,000 genes"),
                 tags$li("Genes ranked by combined bimodality and cluster separation metrics"),
                 tags$li("Quality filters applied: minimum samples (≥20), bimodality threshold (≥0.1), balanced modes (≥10 each)"),
                 tags$li("Differential expression analysis performed using t-tests between clusters with FDR correction")
               ),
               
               h4("Visualization Options"),
               tags$ul(
                 tags$li(strong("Density plots:"), "Show expression distribution patterns and cluster separation"),
                 tags$li(strong("Ranked expression:"), "Display samples ordered by expression level to visualize bimodal pattern and cluster transitions")
               ),
               
               h4("Bimodality Analysis Variables"),
               div(
                 style = "background-color: #f8f9fa; padding: 15px; margin-bottom: 20px; border-radius: 5px;",
                 
                 h5(strong("Core Metrics")),
                 tags$dl(
                   tags$dt(strong("overall_bimodality")),
                   tags$dd("Strength of bimodal distribution across all samples. Calculated as valley depth × peak separation from density estimation. Higher values indicate clearer separation between two expression modes."),
                   
                   tags$dt(strong("cluster_separation_quality")),
                   tags$dd("How well the two proliferation clusters separate into high/low expression modes. Measured as absolute difference in proportions of each cluster in the high expression mode (0-1 scale)."),
                   
                   tags$dt(strong("valley_position")),
                   tags$dd("Expression threshold that separates low and high expression modes. Identified as the deepest valley between the two highest peaks in the density distribution."),
                   
                   tags$dt(strong("cluster1_in_high_mode")),
                   tags$dd("Proportion of cluster 1 samples with expression above the valley position (0-1 scale). Indicates what fraction of cluster 1 shows high expression."),
                   
                   tags$dt(strong("cluster2_in_high_mode")),
                   tags$dd("Proportion of cluster 2 samples with expression above the valley position (0-1 scale). Indicates what fraction of cluster 2 shows high expression."),
                   
                   tags$dt(strong("cluster_explains_bimodality")),
                   tags$dd("Score measuring how well cluster membership explains the bimodal pattern (0-1 scale). Calculated as separation quality × balance factor, where perfect separation yields higher scores."),
                   
                   tags$dt(strong("combined_score")),
                   tags$dd("Overall bimodality quality metric calculated as: overall_bimodality × cluster_explains_bimodality. Prioritizes genes that are both bimodal AND explained by cluster structure.")
                 ),
                 
                 h5(strong("Differential Expression")),
                 tags$dl(
                   tags$dt(strong("ttest_padj")),
                   tags$dd("FDR-adjusted p-value from two-sample t-test comparing expression between proliferation clusters. Values < 0.05 indicate significant differential expression after multiple testing correction.")
                 ),
                 
                 h5(strong("Quality Control")),
                 tags$dl(
                   tags$dt(strong("filter_status")),
                   tags$dd("Indicates which quality filter the gene passed or failed:"),
                   tags$ul(
                     tags$li(strong("passed_all_filters:"), "Gene meets all quality criteria"),
                     tags$li(strong("insufficient_samples:"), "Less than 20 valid expression values"),
                     tags$li(strong("zero_inflated:"), ">30% of values below 0.5 (too low expression)"),
                     tags$li(strong("weak_bimodality:"), "Bimodality score below 0.1 threshold"),
                     tags$li(strong("no_valley_detected:"), "Cannot identify clear valley between modes"),
                     tags$li(strong("unbalanced_modes:"), "Either expression mode has <10 samples")
                   ),
                   
                   tags$dt(strong("passes_filters")),
                   tags$dd("Boolean indicating whether the gene passed all quality filters (TRUE/FALSE)."),
                   
                   tags$dt(strong("quality")),
                   tags$dd("Categorical quality assessment based on combined_score: Excellent (top tier), Good (middle tier), or Fair (lower tier).")
                 ),
                 
                 h5(strong("Ranking Variables")),
                 tags$dl(
                   tags$dt(strong("rank_overall")),
                   tags$dd("Overall rank based on combined_score among all genes that passed filters. Lower numbers = higher quality bimodal genes."),
                   
                   tags$dt(strong("rank_bimodality")),
                   tags$dd("Rank based solely on overall_bimodality score among all analyzed genes. Measures pure bimodality strength regardless of cluster explanation."),
                   
                   tags$dt(strong("rank_separation")),
                   tags$dd("Rank based on cluster_separation_quality among all analyzed genes. Measures how well clusters separate into expression modes."),
                   
                   tags$dt(strong("percentile_overall")),
                   tags$dd("Percentile rank of the combined_score (0-100 scale). Higher percentiles indicate better bimodal quality relative to all analyzed genes.")
                 )
               ),
               
               h4("Gene Signatures"),
               tags$ul(
                 tags$li(strong("Proliferation genes:"), "Cell cycle and DNA replication related genes from established signatures"),
                 tags$li(strong("Infiltration genes:"), "Immune cell marker genes indicating tumor infiltration"),
                 tags$li(strong("Blacklisted:"), "Genes from proliferation + infiltration signatures (may confound analysis)")
               ),
               
               h4("Analysis Workflow"),
               tags$ol(
                 tags$li("Load BaSq samples and cluster by proliferation score (K-means, k=2)"),
                 tags$li("For each gene, estimate expression density using Sheather-Jones bandwidth"),
                 tags$li("Identify peaks and valleys in density to calculate bimodality strength"),
                 tags$li("Find optimal valley position separating high/low expression modes"),
                 tags$li("Calculate cluster separation quality based on mode assignments"),
                 tags$li("Perform two-sample t-tests between clusters with FDR correction"),
                 tags$li("Apply quality filters to ensure robust bimodal patterns"),
                 tags$li("Rank genes by combined bimodality and cluster explanation scores"),
                 tags$li("Assign quality categories and generate final rankings")
               ),
               
               h4("Usage Tips"),
               tags$ul(
                 tags$li("Use ", strong("Gene Explorer"), " to visualize individual genes and see bimodal patterns"),
                 tags$li("Browse ", strong("Bimodality Results"), " table to discover high-quality candidates"),
                 tags$li("Filter by ", strong("'Passed filters only'"), " to focus on robust bimodal genes"),
                 tags$li("Sort by ", strong("combined_score"), " to find genes where clusters best explain bimodality"),
                 tags$li("Check ", strong("filter_status"), " to understand why certain genes failed quality checks"),
                 tags$li("Look for genes with low ", strong("ttest_padj"), " values for significant differential expression"),
                 tags$li("Click table rows to automatically load genes in the visualizer"),
                 tags$li("Use ", strong("ranked expression plot"), " to clearly see bimodal transitions and cluster patterns")
               ),
               
               h4("Technical Notes"),
               tags$ul(
                 tags$li("Expression data is log2-transformed and normalized"),
                 tags$li("Density estimation uses Sheather-Jones bandwidth selection for robustness"),
                 tags$li("Bimodality detection requires at least 2 clear peaks with significant valley depth"),
                 tags$li("Cluster separation is measured independently of initial bimodality detection"),
                 tags$li("T-tests assume normal distribution; FDR correction controls false discovery rate"),
                 tags$li("Combined scoring emphasizes genes where proliferation clusters drive bimodal expression")
               )
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
  
  output$passed_genes <- renderText({
    sum(binomial_results$passes_filters)
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
    
    # Get filter status and t-test results if available
    filter_info <- ""
    ttest_info <- ""
    if (gene_name %in% binomial_results$gene) {
      gene_row <- binomial_results[binomial_results$gene == gene_name, ]
      filter_status <- gene_row$filter_status
      passes_filters <- gene_row$passes_filters
      ttest_padj <- gene_row$ttest_padj
      
      filter_info <- paste0("\nFilter Status: ", filter_status, 
                            "\nPassed All Filters: ", ifelse(passes_filters, "Yes", "No"))
      
      # Format t-test p-value
      if (!is.na(ttest_padj)) {
        if (ttest_padj < 0.001) {
          ttest_display <- "<0.001"
        } else {
          ttest_display <- round(ttest_padj, 4)
        }
        significance <- ifelse(ttest_padj < 0.05, "Yes", "No")
        ttest_info <- paste0("\nT-test FDR p-value: ", ttest_display,
                             "\nSignificant (FDR<0.05): ", significance)
      } else {
        ttest_info <- "\nT-test FDR p-value: Not available"
      }
    }
    
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
        "Blacklisted: ", ifelse(gene_name %in% blacklisted_genes, "Yes", "No"),
        filter_info,
        ttest_info, "\n\n",
        "Expression Statistics: No valid data available"
      ))
    }
    
    paste0(
      "Gene: ", gene_name, "\n",
      "Signature: ", signature_status, "\n", 
      "Blacklisted: ", ifelse(gene_name %in% blacklisted_genes, "Yes", "No"),
      filter_info,
      ttest_info, "\n\n",
      "Expression Statistics:\n",
      "  Range: ", round(min(expr_vals), 2), " - ", round(max(expr_vals), 2), "\n",
      "  Mean: ", round(mean(expr_vals), 2), "\n",
      "  Median: ", round(median(expr_vals), 2), "\n",
      "  Samples: ", length(expr_vals)
    )
  })
  
  # Generate the plot based on selected type
  output$combined_plot <- renderPlot({
    req(formatted_gene(), gene_exists())
    
    tryCatch({
      if (input$plot_type == "density") {
        # Original density plot
        plot_bimodal_gene(this_gene = formatted_gene(),
                          this_expr = predicted_imvigor$data,
                          this_cluster = basq_proliferations)
      } else {
        # New ranked expression plot
        plot_ranked_expression(this_gene = formatted_gene(),
                               this_expr = predicted_imvigor$data,
                               this_cluster = basq_proliferations)
      }
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
      
      # Format t-test p-value for display
      ttest_padj_display <- if (!is.na(gene_stats$ttest_padj)) {
        if (gene_stats$ttest_padj < 0.001) {
          "<0.001"
        } else {
          round(gene_stats$ttest_padj, 4)
        }
      } else {
        "Not available"
      }
      
      data.frame(
        Metric = c("Filter Status", "Passed All Filters", "Overall Bimodality", "Cluster Separation Quality", 
                   "Valley Position", "Cluster 1 in High Mode (%)", 
                   "Cluster 2 in High Mode (%)", "Cluster Explains Bimodality", 
                   "Combined Score", "T-test FDR p-value", "Rank Overall", "Quality"),
        Value = c(
          gene_stats$filter_status,
          ifelse(gene_stats$passes_filters, "Yes", "No"),
          round(gene_stats$overall_bimodality, 4),
          round(gene_stats$cluster_separation_quality, 4),
          ifelse(is.na(gene_stats$valley_position), "Not detected", round(gene_stats$valley_position, 4)),
          paste0(round(gene_stats$cluster1_in_high_mode * 100, 1), "%"),
          paste0(round(gene_stats$cluster2_in_high_mode * 100, 1), "%"),
          round(gene_stats$cluster_explains_bimodality, 4),
          round(gene_stats$combined_score, 4),
          ttest_padj_display,
          ifelse(gene_stats$passes_filters, gene_stats$rank_overall, "Not ranked"),
          gene_stats$quality
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
    
    # Apply filter status
    if (input$filter_status == "passed_only") {
      data <- data %>% filter(passes_filters == TRUE)
    } else if (input$filter_status == "failed_only") {
      data <- data %>% filter(passes_filters == FALSE)
    }
    
    # Apply numeric filters
    data <- data %>%
      filter(
        overall_bimodality >= input$min_bimodality,
        cluster_separation_quality >= input$min_separation
      )
    
    # Format for display - include ALL columns
    data %>%
      select(
        Gene = gene,
        Signature = signature,
        `Filter Status` = filter_status,
        `Passes Filters` = passes_filters,
        `Bimodality Score` = overall_bimodality,
        `Cluster Separation` = cluster_separation_quality,
        `Valley Position` = valley_position,
        `Cluster 1 High %` = cluster1_in_high_mode,
        `Cluster 2 High %` = cluster2_in_high_mode,
        `Explains Bimodality` = cluster_explains_bimodality,
        `Combined Score` = combined_score,
        `T-test FDR p-value` = ttest_padj,
        `Rank Overall` = rank_overall,
        `Rank Bimodality` = rank_bimodality,
        `Rank Separation` = rank_separation,
        `Percentile Overall` = percentile_overall,
        Quality = quality
      ) %>%
      mutate(
        `Bimodality Score` = round(`Bimodality Score`, 4),
        `Cluster Separation` = round(`Cluster Separation`, 4),
        `Valley Position` = ifelse(is.na(`Valley Position`), "Not detected", round(`Valley Position`, 4)),
        `Cluster 1 High %` = paste0(round(`Cluster 1 High %` * 100, 1), "%"),
        `Cluster 2 High %` = paste0(round(`Cluster 2 High %` * 100, 1), "%"),
        `Explains Bimodality` = round(`Explains Bimodality`, 4),
        `Combined Score` = round(`Combined Score`, 4),
        `T-test FDR p-value` = ifelse(is.na(`T-test FDR p-value`), "Not available", 
                                      ifelse(`T-test FDR p-value` < 0.001, "<0.001", round(`T-test FDR p-value`, 4))),
        `Percentile Overall` = round(`Percentile Overall`, 4),
        `Passes Filters` = ifelse(`Passes Filters`, "Yes", "No"),
        `Rank Overall` = ifelse(`Passes Filters` == "Yes", `Rank Overall`, "Not ranked")
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
        dom = 'frtip',
        columnDefs = list(
          list(width = '100px', targets = c(0, 1)),  # Gene and Signature columns
          list(width = '120px', targets = c(2)),     # Filter Status column
          list(width = '80px', targets = c(3)),      # Passes Filters column
          list(width = '100px', targets = c(11)),    # T-test FDR p-value column
          list(width = '80px', targets = c(12:16))   # Rank columns
        )
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
      ) %>%
      DT::formatStyle(
        'Filter Status',
        backgroundColor = DT::styleEqual(
          c('passed_all_filters', 'insufficient_samples', 'zero_inflated', 'weak_bimodality', 'no_valley_detected', 'unbalanced_modes'),
          c('#d4edda', '#f8d7da', '#ffeaa7', '#fdcb6e', '#e17055', '#fd79a8')
        )
      ) %>%
      DT::formatStyle(
        'Passes Filters',
        backgroundColor = DT::styleEqual(
          c('Yes', 'No'),
          c('#d4edda', '#f8d7da')
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
      plot_type_name <- ifelse(input$plot_type == "density", "density", "ranked")
      paste0(formatted_gene(), "_", plot_type_name, "_plot.png")
    },
    content = function(file) {
      req(formatted_gene(), gene_exists())
      
      if (input$plot_type == "density") {
        plot <- plot_bimodal_gene(this_gene = formatted_gene(),
                                  this_expr = predicted_imvigor$data,
                                  this_cluster = basq_proliferations)
      } else {
        plot <- plot_ranked_expression(this_gene = formatted_gene(),
                                       this_expr = predicted_imvigor$data,
                                       this_cluster = basq_proliferations)
      }
      
      ggsave(file, plot, width = 12, height = 6, dpi = 300)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)