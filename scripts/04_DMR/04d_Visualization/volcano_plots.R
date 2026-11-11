# Modify form Rolph Vedder's code
#################################################################
# Helpers                                                       #
#################################################################

select.columns.in.order <- function(dataframe, columns) {
  dataframe[, columns]
}

drop.columns.if.all.same.value <- function(dataframe) {
  for (name in colnames(dataframe)) {
    is.all.same <- (dataframe[, name] %>% unique() %>% length()) <= 1
    if (is.all.same) {
      dataframe <- dataframe %>%
        dplyr::select(
          -tidyselect::one_of(name)
        )
    }
  }
  dataframe
}

normalize.ignore.na <- function(df) {
  for (row in 1:nrow(df)) {
    this.row <- as.numeric(df[row,])
    this.row.mean <- mean(this.row, na.rm=TRUE)
    df[row, !is.na(this.row)] <- df[row, !is.na(this.row)] - this.row.mean
  }
  df
}

select.rows.in.order <- function(dataframe, rows) {
  dataframe[rows,]
}

#################################################################
# Cutoffs etd                                                   #
#################################################################

significance.cutoff <- 0.05 # FDR
fc.cutoff <- 0.5

significance.filter <- (function() {
  sco <- significance.cutoff
  fco <- fc.cutoff
  
  function(df) {
    df %>%
      dplyr::mutate(
        fc.include = dplyr::case_when(
          fco == 0  ~ TRUE,
          TRUE ~ abs(logFC) > fco
        )
      ) %>%
      dplyr::filter(
        adj.P.Val <= sco &
          fc.include
      ) %>%
      dplyr::select(-fc.include)
  }
})()

volcanoplot.mutate <- function(df) {
  df %>%
    dplyr::mutate(
      Legend.raw = dplyr::case_when(
        fc.cutoff == 0 ~ TRUE,
        TRUE ~ abs(logFC) > fc.cutoff
      ),
      Legend = dplyr::case_when(
        adj.P.Val > significance.cutoff ~ "Not.Significant",
        logFC < 0 & Legend.raw ~ "Hypomethylated",
        logFC > 0 & Legend.raw ~ "Hypermethylated",
        TRUE ~ "Not.Significant"
      )
    ) %>%
    dplyr::select(-Legend.raw)
}

volcanoplot.fc.filter <- function(df) {
  if (fc.cutoff > 0) {
    df %>%
      dplyr::filter(
        logFC >= fc.cutoff
      )
  } else {
    df
  }
}

volcanoplot.vertical.fc.lines <- function(plt) {
  if (fc.cutoff > 0) {
    plt +
      geom_vline(
        xintercept = fc.cutoff,
        colour = "#990000",
        linetype = "dashed"
      ) +
      geom_vline(
        xintercept = (fc.cutoff)*(-1),
        colour = "#990000",
        linetype = "dashed"
      )
  } else {
    plt
  }
}

###############
# Volcano Plots
###############
#c.result.dir <- file.path(script.result.dir, "volcano plots")
#dir.create(c.result.dir, recursive=T)

make.volcano.plot <- function(limma.result.to.plot) {
   fdr.threshold <- max(
    limma.result.to.plot$P.Value[limma.result.to.plot$adj.P.Val < significance.cutoff],
    na.rm = TRUE
  )
  y.bar.intercept1 <- -log10(fdr.threshold)
  y.bar.intercept2 <- -log10(significance.cutoff)
  p <- (limma.result.to.plot %>%
          ggplot2::ggplot(
            mapping = ggplot2::aes(
              x = logFC,
              y = -log10(P.Value)
            )
          ) +
          ggplot2::geom_point(
            mapping = ggplot2::aes(color = Legend)
          ) +
          #ggplot2::geom_hline(
          #    yintercept = y.bar.intercept2,
          #    colour = "black",
          #    linetype = "dashed"  
          #) +
          ggplot2::geom_hline(
            yintercept = y.bar.intercept1,
            colour = "#990000",
            linetype = "dashed"
    )) %>%
    volcanoplot.vertical.fc.lines() +
    ggrepel::geom_text_repel(
      data = limma.result.to.plot[1:10,],
      mapping = ggplot2::aes(label = CpG.id),
      size = 2,# size of top10 significant
      box.padding = unit(0.1, "lines"),
      point.padding = unit(0.1, "lines"),
      max.overlaps = Inf)+
    ggplot2::ggtitle("Volcano plot") +
    ggplot2::scale_color_manual(
      values = c(
        Not.Significant = "grey",
        Hypomethylated = "blue",
        Hypermethylated = "red"
      )
    ) +
    ggplot2::theme_bw(base_size = 20) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white"),  # Ensure background is white
      axis.title.x = ggplot2::element_text(size = 24, face = "bold", margin = margin(t = 10)),
      axis.title.y = ggplot2::element_text(size = 24, face = "bold", margin = margin(r = 10)),
      axis.text.x = ggplot2::element_text(size = 18, color = "black"),
      axis.text.y = ggplot2::element_text(size = 18, color = "black"),
      panel.grid = element_blank(),
      legend.text = ggplot2::element_text(size = 24),    # Legend text size
      legend.title = ggplot2::element_text(size = 14, face = "bold"),  # Legend title size
      plot.margin = margin(4, 4, 4, 4, "cm")  # Increased bottom margin
    ) +
    xlab("log2FC") +
    #ggprism::theme_prism() +
    #ggprism::scale_colour_prism() +
    ggprism::scale_shape_prism()
  
  return(p)
}