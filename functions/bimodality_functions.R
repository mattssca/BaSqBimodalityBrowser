library(ggplot2)
library(dplyr)

# Function to plot gene expression distribution by clusters
plot_gene_distribution <- function(expr_matrix, cluster_data, gene_name, 
                                   plot_type = "density", add_stats = TRUE, 
                                   show_clusters = TRUE) {
  
  # Prepare data
  cluster_labels <- setNames(cluster_data$cluster_label, cluster_data$sample_id)
  common_samples <- intersect(colnames(expr_matrix), names(cluster_labels))
  
  # Check if gene exists
  if (!gene_name %in% rownames(expr_matrix)) {
    cat("Gene", gene_name, "not found in expression matrix.\n")
    cat("Available genes starting with", substr(gene_name, 1, 3), ":\n")
    matching_genes <- rownames(expr_matrix)[grepl(paste0("^", substr(gene_name, 1, 3)), 
                                                  rownames(expr_matrix), ignore.case = TRUE)]
    print(head(matching_genes, 10))
    return(NULL)
  }
  
  # Extract gene expression
  gene_expr <- expr_matrix[gene_name, common_samples]
  clusters_subset <- cluster_labels[common_samples]
  
  # Create plotting dataframe
  plot_data <- data.frame(
    sample_id = common_samples,
    expression = as.numeric(gene_expr),
    cluster = clusters_subset,
    proliferation_score = cluster_data$proliferation_score[
      match(common_samples, cluster_data$sample_id)
    ]
  )
  
  # Calculate statistics
  if (show_clusters) {
    stats_summary <- plot_data %>%
      group_by(cluster) %>%
      summarise(
        n = n(),
        mean_expr = mean(expression, na.rm = TRUE),
        median_expr = median(expression, na.rm = TRUE),
        sd_expr = sd(expression, na.rm = TRUE),
        min_expr = min(expression, na.rm = TRUE),
        max_expr = max(expression, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Statistical test
    t_test <- t.test(expression ~ cluster, data = plot_data)
    wilcox_test <- wilcox.test(expression ~ cluster, data = plot_data)
    
    # From your binomial results, get binomial test info if available
    binomial_info <- NULL
    if (exists("binomial_results_fast")) {
      binomial_row <- binomial_results_fast[binomial_results_fast$gene == gene_name, ]
      if (nrow(binomial_row) > 0) {
        binomial_info <- binomial_row[1, ]
      }
    }
  } else {
    # Overall statistics without clustering
    stats_summary <- data.frame(
      n = nrow(plot_data),
      mean_expr = mean(plot_data$expression, na.rm = TRUE),
      median_expr = median(plot_data$expression, na.rm = TRUE),
      sd_expr = sd(plot_data$expression, na.rm = TRUE),
      min_expr = min(plot_data$expression, na.rm = TRUE),
      max_expr = max(plot_data$expression, na.rm = TRUE)
    )
  }
  
  # Create the plot based on type
  if (plot_type == "density") {
    
    if (show_clusters) {
      p <- ggplot(plot_data, aes(x = expression)) +
        
        # Density curves
        geom_density(aes(fill = cluster, color = cluster), alpha = 0.6, size = 1.2) +
        
        # Rug plot
        geom_rug(aes(color = cluster), alpha = 0.7, size = 0.8) +
        
        # Vertical lines for means
        geom_vline(data = stats_summary, 
                   aes(xintercept = mean_expr, color = cluster),
                   linetype = "dashed", size = 1.2) +
        
        scale_fill_manual(values = c("cluster1" = "#F87B1B", "cluster2" = "#134686")) +
        scale_color_manual(values = c("cluster1" = "#9E4705", "cluster2" = "#09213E")) +
        
        labs(
          title = paste("Expression Distribution:", gene_name),
          subtitle = paste0("t-test p-value: ", format.pval(t_test$p.value, digits = 3),
                            " | Wilcoxon p-value: ", format.pval(wilcox_test$p.value, digits = 3)),
          x = "Log2(TPM + 1)",
          y = "Density",
          fill = "Cluster",
          color = "Cluster"
        )
      
    } else {
      # Plot without cluster information
      p <- ggplot(plot_data, aes(x = expression)) +
        
        # Single density curve
        geom_density(fill = "lightgray", color = "black", alpha = 0.6, size = 1.2) +
        
        # Rug plot
        geom_rug(alpha = 0.7, size = 0.8, color = "black") +
        
        # Vertical line for mean
        geom_vline(xintercept = stats_summary$mean_expr,
                   linetype = "dashed", size = 1.2, color = "red") +
        
        labs(
          title = paste("Expression Distribution:", gene_name),
          subtitle = paste0("Mean: ", round(stats_summary$mean_expr, 3), 
                            " | SD: ", round(stats_summary$sd_expr, 3)),
          x = "Log2(TPM + 1)",
          y = "Density"
        )
    }
    
    # Overall median (for both cases)
    p <- p + 
      geom_vline(xintercept = median(plot_data$expression), 
                 linetype = "dotted", color = "black", size = 1) +
      
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "top"
      )
    
  } else if (plot_type == "boxplot") {
    
    if (show_clusters) {
      p <- ggplot(plot_data, aes(x = cluster, y = expression, fill = cluster)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
        
        scale_fill_manual(values = c("cluster1" = "#f46517", "cluster2" = "#113473")) +
        
        labs(
          title = paste("Expression Distribution:", gene_name),
          subtitle = paste0("t-test p-value: ", format.pval(t_test$p.value, digits = 3)),
          x = "Cluster", 
          y = "Log2(TPM + 1)",
          fill = "Cluster"
        )
    } else {
      p <- ggplot(plot_data, aes(x = "", y = expression)) +
        geom_boxplot(fill = "lightgray", alpha = 0.7, outlier.alpha = 0.6) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
        
        labs(
          title = paste("Expression Distribution:", gene_name),
          subtitle = paste0("Median: ", round(stats_summary$median_expr, 3)),
          x = "", 
          y = "Log2(TPM + 1)"
        )
    }
    
    p <- p + 
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "none"
      )
    
  } else if (plot_type == "violin") {
    
    if (show_clusters) {
      p <- ggplot(plot_data, aes(x = cluster, y = expression, fill = cluster)) +
        geom_violin(alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        
        scale_fill_manual(values = c("cluster1" = "#f46517", "cluster2" = "#113473")) +
        
        labs(
          title = paste("Expression Distribution:", gene_name),
          subtitle = paste0("Mean difference: ", 
                            round(diff(stats_summary$mean_expr), 3)),
          x = "Cluster",
          y = "Log2(TPM + 1)",
          fill = "Cluster"
        )
    } else {
      p <- ggplot(plot_data, aes(x = "", y = expression)) +
        geom_violin(fill = "lightgray", alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8) +
        geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
        
        labs(
          title = paste("Expression Distribution:", gene_name),
          subtitle = paste0("Range: [", round(stats_summary$min_expr, 2), ", ",
                            round(stats_summary$max_expr, 2), "]"),
          x = "",
          y = "Log2(TPM + 1)"
        )
    }
    
    p <- p + 
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "none"
      )
  }
  
  # Add statistical annotations if requested
  if (add_stats && show_clusters) {
    
    # Calculate positions for annotations
    y_max <- max(density(plot_data$expression)$y) * 1.1
    x_min <- min(plot_data$expression)
    x_max <- max(plot_data$expression)
    
    if (plot_type == "density") {
      # Add binomial test info in top left (BLACK TEXT)
      if (!is.null(binomial_info)) {
        p <- p +
          annotate("text",
                   x = x_min + (x_max - x_min) * 0.05, y = y_max * 0.95,
                   label = paste0("Binomial Test:\n",
                                  "Prop C1: ", round(binomial_info$cluster1_prop, 3), "\n",
                                  "Prop C2: ", round(binomial_info$cluster2_prop, 3), "\n",
                                  "p-adj: ", format.pval(binomial_info$p_adjusted, digits = 2)),
                   hjust = 0, vjust = 1,
                   color = "black", size = 3, fontface = "bold",
                   bbox = list(boxstyle = "round,pad=0.3", facecolor = "white", alpha = 0.8))
      }
      
      # Add cluster statistics in top right
      p <- p + 
        annotate("text", 
                 x = x_max - (x_max - x_min) * 0.05, y = y_max * 0.95,
                 label = paste0("Cluster 1 (n=", stats_summary$n[1], ")\n",
                                "Mean: ", round(stats_summary$mean_expr[1], 2), "\n",
                                "SD: ", round(stats_summary$sd_expr[1], 2)),
                 hjust = 1, vjust = 1, 
                 color = "#9E4705", size = 3.5, fontface = "bold",
                 bbox = list(boxstyle = "round,pad=0.3", facecolor = "white", alpha = 0.8)) +
        
        annotate("text", 
                 x = x_max - (x_max - x_min) * 0.05, y = y_max * 0.65,
                 label = paste0("Cluster 2 (n=", stats_summary$n[2], ")\n",
                                "Mean: ", round(stats_summary$mean_expr[2], 2), "\n", 
                                "SD: ", round(stats_summary$sd_expr[2], 2)),
                 hjust = 1, vjust = 1,
                 color = "#09213E", size = 3.5, fontface = "bold",
                 bbox = list(boxstyle = "round,pad=0.3", facecolor = "white", alpha = 0.8))
    }
  }
  
  # Print summary statistics
  if (add_stats) {
    cat("\n=== GENE EXPRESSION SUMMARY:", gene_name, "===\n")
    
    if (show_clusters) {
      print(stats_summary)
      
      cat("\nStatistical Tests:\n")
      cat("t-test p-value:", format.pval(t_test$p.value), "\n")
      cat("Wilcoxon p-value:", format.pval(wilcox_test$p.value), "\n")
      cat("Mean difference:", round(diff(stats_summary$mean_expr), 3), "\n")
      cat("Effect size (Cohen's d):", round(abs(diff(stats_summary$mean_expr)) / 
                                              sqrt(mean(stats_summary$sd_expr^2)), 3), "\n")
      
      if (!is.null(binomial_info)) {
        cat("\nBinomial Test Results:\n")
        cat("Cluster 1 high proportion:", round(binomial_info$cluster1_prop, 3), "\n")
        cat("Cluster 2 high proportion:", round(binomial_info$cluster2_prop, 3), "\n") 
        cat("Effect size:", round(binomial_info$effect_size, 3), "\n")
        cat("Adjusted p-value:", format.pval(binomial_info$p_adjusted), "\n")
      }
    } else {
      print(stats_summary)
      cat("\nOverall distribution statistics without clustering shown above.\n")
    }
  }
  
  return(p)
}

plot_bimodal_gene = function(this_gene = NULL, this_expr = NULL, this_cluster = NULL){
  p1 = plot_gene_distribution(expr_matrix = this_expr, cluster_data = this_cluster, gene_name = this_gene, plot_type = "density", show_clusters = FALSE)
  p2 = plot_gene_distribution(expr_matrix = this_expr, cluster_data = this_cluster, gene_name = this_gene, plot_type = "density", show_clusters = TRUE)
  combined_plot <- p1 + p2
  
  return(combined_plot)
}
