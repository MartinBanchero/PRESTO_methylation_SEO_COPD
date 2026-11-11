# Box plots
boxplots.top.cpgs <- function(top.sites, values, type, fdr, pheno.color.var){
par(mfrow = c(5, 3), cex.lab=1.2, cex.axis=1.2, cex.main=3)
  for (cpg in top.sites) {
  boxplot(
    values[cpg, ] ~ pheno.color.var,
    col = c("darkred", "cadetblue"),
    ylab = paste0("methylation (", type, "-values)"),
    xlab = "Sample type",
    main = "" 
    #ylim = c(0, 1)
  )
  title(main = cpg, cex.main = 1.6)
 }
}
### Box plot for top genes
boxplots.top.genes <- function(genes, values, type, pheno.color.var, id.mapping.df, 
                               ensembl.id.col, symbol.col) {
  
  # Create mapping from symbol to ENSEMBL ID
  mapping <- id.mapping.df[, c(ensembl.id.col, symbol.col)]
  colnames(mapping) <- c("ensembl", "symbol")
  
  # Find ENSEMBL IDs for requested genes
  ensembl_ids <- c()
  gene_labels <- c()
  
  for (gene in genes) {
    # Check if gene is already an ENSEMBL ID in the matrix
    if (gene %in% rownames(values)) {
      ensembl_ids <- c(ensembl_ids, gene)
      gene_labels <- c(gene_labels, gene)
    } 
    # Check if gene is a symbol that can be mapped
    else {
      matched <- mapping[mapping$symbol == gene & mapping$ensembl %in% rownames(values), ]
      if (nrow(matched) > 0) {
        ensembl_ids <- c(ensembl_ids, matched$ensembl[1])
        gene_labels <- c(gene_labels, gene)
      } else {
        warning(paste("Gene not found or not mappable:", gene))
      }
    }
  }
  
  if (length(ensembl_ids) == 0) {
    stop("None of the genes could be found or mapped to ENSEMBL IDs in the matrix")
  }
  
  # Set up plot layout
  n_genes <- length(ensembl_ids)
  n_cols <- min(3, n_genes)
  n_rows <- ceiling(n_genes / n_cols)
  
  par(mfrow = c(n_rows, n_cols), 
      cex.lab = 1.8, 
      cex.axis = 1.6, 
      cex.main = 2.2,
      mar = c(6, 6, 5, 2) + 0.1)
  
  for (i in 1:length(ensembl_ids)) {
    expr_values <- as.numeric(values[ensembl_ids[i], ])
    
    boxplot(
      expr_values ~ pheno.color.var,
      col = c("darkred", "cadetblue"),
      ylab = paste0("Gene Expression (", type, ")"),
      xlab = "Sample type",
      main = paste0("GEX for gene: ", gene_labels[i]),
      cex.lab = 2.0,      # Explicitly set label size
      cex.axis = 1.8,     # Explicitly set axis size
      cex.main = 2.4      # Explicitly set title size
    )
  }
}
# Swarm plots
swarm.nocolor.plot <- function(top.sites, values, m.or.beta ,  pheno.color.var){
par(mfrow = c(5, 3), cex.lab=1.2, cex.axis=1.2, cex.main=3)
for (cpg in top.sites) {
  plotCpg(
    values[cpg, ,drop = FALSE],
    type = "categorical",
    cpg = cpg,
    measure = m.or.beta,
    pheno = pheno.color.var,
    ylab = paste0("methylation (", m.or.beta, "-values)"),
  )
 }
}


grid.swarm.nocolor.plot <- function(top.sites, values, pheno.color.var) {
  par(mfrow = c(5, 3), cex.lab = 1.2, cex.axis = 1.2, cex.main = 3)
  
  # Ensure categories align with m.values
  if (!all(names(pheno.color.var) %in% base::colnames(values))) {
    stop("Sample names in pheno.color.var and m.values do not match.")
  }
  categories <- pheno.color.var[base::colnames(values)]
  
  for (cpg in top_sites) {
    # Check if CpG exists in m.values
    if (!cpg %in% rownames(values)) {
      warning(paste("CpG site", cpg, "not found in m.values"))
      next
    }
    
    # Extract data for the current CpG
    data <- as.vector(values[cpg, ])
    
    # Debugging: Check lengths of data and categories
    if (length(data) != length(categories)) {
      warning(paste("Data and categories length mismatch for CpG:", cpg))
      next
    }
    
    # Create an empty plot
    plot(1, type = "n", xlab = "Category", ylab = paste0("methylation (", type, "-values)"), 
         xlim = c(0.5, length(unique(categories)) + 0.5), 
         ylim = range(data, na.rm = TRUE), main = cpg)
    
    # Add swarm plot
    beeswarm::beeswarm(data ~ categories, add = TRUE, pch = 16, col = "blue")
    
    # Add boxplot
    boxplot(data ~ categories, 
            add = TRUE,
            outline = FALSE,
            col = "gray",
            border = "black",
            xaxt = "n")
  }
}

library(ggplot2)
library(gridExtra)

grid.swarm.color.plot<- function(top_sites, values, pheno.data, pheno.color.var, filename){
  # Ensure colnames of m.values are correctly assigned
  base::colnames(values) <- pheno.data$Sample_Name
  
  # Initialize an empty list to store individual plots
  plot_list <- list()
  
  # Loop through each CpG site and create a plot, adding it to plot_list
  for (i in seq_along(top_sites)) {
    cpg <- top_sites[i]
    
    # Prepare data frame for plotting
    data_to_plot <- data.frame(
      SampleID = base::colnames(values),
      Mvalues = as.numeric(values[cpg, , drop = FALSE]),
      SEO_COPD_Y_or_N = pheno.color.var
    )
    
    # Generate the plot for this CpG site
    p <- ggplot(data_to_plot, aes(x = SEO_COPD_Y_or_N, y = Mvalues, color = SampleID)) +
      geom_point(size = 5) +
      geom_text(aes(label = SampleID), vjust = -0.5, hjust = 0.5, size = 3) +  # Add labels
      labs(
        x = "COPD Status",
        y = "Mvalues",
        title = paste("CpG:", cpg)
      ) +
      theme_bw() +
      scale_color_discrete(name = "Sample ID")
    
    # Add the plot to the list
    plot_list[[i]] <- p
    
    # Stop if we have reached 9 plots (3x3 grid)
    if (i == 9) break
  }
  
  # Arrange the plots in a 3x3 grid
  combined_plot <- grid.arrange(grobs = plot_list, ncol = 5, nrow = 5)
  
  # Save the combined plot
  ggsave(
    filename = filename,
    plot = combined_plot,
    device = "png",
    dpi = 300,
    width = 20,    # Width for the combined plot
    height = 20    # Height for the combined plot
   )
}

# Cell proportions
boxplot.cell.proportions<- function(df, x.cells, y.condition){
plt <- ggplot(
  df_long,
  aes(factor(x.cells), props, fill = factor(y.condition))
) + 
  geom_boxplot(outlier.size = 1, outlier.shape = 21, color = "black", alpha = 0.7) +  # Adjusted boxplot styling
  labs(
    title = "Box Plot of Cell Type Proportions",
    x = "Cell Types",
    y = "Proportion",
    fill = "COPD Status"
  ) +
  scale_fill_brewer(palette = "Set2") +  # Use a color-blind-friendly palette
  theme_bw(base_size = 16) +  # Set base font size for readability
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center and bold title
    axis.title.y = element_text(size = 18, face = "bold"),             # Larger, bold y-axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),      # Rotate x-axis text for clarity
    axis.text.y = element_text(size = 12),                             # Larger y-axis text
    legend.title = element_text(size = 14, face = "bold"),             # Larger, bold legend title
    legend.text = element_text(size = 12),                             # Larger legend text
    legend.position = "top"                                            # Position legend on top for balance
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))  # Center legend title


# Save the combined plot
ggsave(
  filename = filename,
  plot = plt,
  device = "png",
  dpi = 300,
  width = 6,    # Width for the combined plot
  height = 6    # Height for the combined plot
)

}

###############################
# Swarm plot with group means #
###############################
library(ggplot2)
library(ggbeeswarm)

swarm.with.group.means <- function(top.sites, values, pheno.color.var) {
  # Create a grid layout for multiple plots
  plots <- list()
  
  for (cpg in top.sites) {
    # Create data frame for current CpG
    df <- data.frame(
      Methylation = values[cpg, ],
      Group = pheno.color.var
    )
    
    # Calculate group means
    means <- aggregate(Methylation ~ Group, data = df, mean, na.rm = TRUE)
    
    # Create plot
    p <- ggplot(df, aes(x = Group, y = Methylation, color = Group)) +
      geom_beeswarm(cex = 3) +  # Swarm plot
      stat_summary(fun=mean, geom="crossbar",color="grey", size=0.5) +
      stat_summary(fun=mean, geom="point", color="black", size=2) +
      labs(title = cpg, y = "Methylation Value (betas)") +
      #theme_minimal() +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 12, face = "bold"))
    plots[[cpg]] <- p
  }
  
  # Arrange plots in grid
  do.call(grid.arrange, c(plots, ncol = 3, nrow = 5))
}

