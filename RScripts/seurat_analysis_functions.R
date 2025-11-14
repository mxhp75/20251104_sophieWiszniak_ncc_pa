
# Function to add metadata to a Seurat object ---------------------------------------------------------------------

add_metadata <- function(object, metadata){
  object[[]] <-
    object[[]] %>%
    as_tibble(rownames = "barcode") %>%
    left_join(y = metadata, by = "barcode") %>%
    column_to_rownames(var = "barcode")
  
  return(object)
}



# Convert human genes to mouse --------------------------------------------

# Use this for converting Seurat's list of cell cycle genes to mouse genes
# https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/

# Note that this isn't based on biology; the genes were derived from human experiments
# but the developers still recommend this solution:
# https://github.com/satijalab/seurat/issues/2493#issuecomment-575702274

human_genes_to_mouse <- function(human_genes){
  hom_data <- homologene::homologene(human_genes, inTax = 9606, outTax = 10090)
  result <- setNames(hom_data$`10090`, hom_data$`9606`)
  return(result)
}

# Check
# s.genes <- human_genes_to_mouse(Seurat::cc.genes.updated.2019$s.genes)
# g2m.genes <- human_genes_to_mouse(Seurat::cc.genes.updated.2019$g2m.genes)
# 
# (s.genes %in% (map(seurat_list, rownames) %>% unlist())) %>% table()
# (g2m.genes %in% (map(seurat_list, rownames) %>% unlist())) %>% table()
# All are found in at least one Seurat object



# Function to create QC metrics then plot them --------------------------------------------------------------------

perform_qc <- function(object, n_dims = 30, clustering_resolution = 0.8, grouping_variable = NULL, dim_reduction = FALSE, plot_out_dir){
  
  ##########
  ### Checks
  ##########
  
  assertthat::assert_that(class(object) == "Seurat",
                          msg = "Object must be a Seurat object")
  
  assertthat::assert_that(rlang::is_integerish(n_dims),
                          msg = "n_dims must be a whole number")
  
  assertthat::assert_that(is.numeric(clustering_resolution),
                          msg = "Clustering resolution must be numeric")
  
  if(!is.null(grouping_variable)){
    assertthat::assert_that(assertthat::is.string(grouping_variable),
                            grouping_variable %in% colnames(object[[]]),
                            msg = "Grouping variable must be a string, and be present in the Seurat object metadata") 
  }
  
  assertthat::assert_that(is.logical(dim_reduction),
                          msg = "dim_reduction must be TRUE or FALSE (No quotation marks)")

  
  
  #########
  ### Setup
  #########
  
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.most_abundant", "complexity"#,
                  # "log10complexity"
                  )
  cell_cycle_metrics <- c("S.Score", "G2M.Score")
  all_qc_metrics <- c(qc_metrics, cell_cycle_metrics)
  
  s.genes <- human_genes_to_mouse(Seurat::cc.genes.updated.2019$s.genes)
  g2m.genes <- human_genes_to_mouse(Seurat::cc.genes.updated.2019$g2m.genes)
  
  # Generated using:
  # ~/localwork/facility/projects/250408SophieWiszniak-mouse_heart_Mib1_KD_scRNA/02-download_mouse_ribosomal_genes.R
  ribo_genes <- c("Rpl10", "Rpl10a", "Rpl11", "Rpl12", "Rpl13", "Rpl13a", "Rpl14", 
                  "Rpl15", "Rpl17", "Rpl18", "Rpl18a", "Rpl19", "Rpl21", "Rpl22", 
                  "Rpl23", "Rpl23a", "Rpl24", "Rpl26", "Rpl27", "Rpl27a", "Rpl28", 
                  "Rpl29", "Rpl3", "Rpl30", "Rpl31", "Rpl32", "Rpl34", "Rpl35", 
                  "Rpl35a", "Rpl36", "Rpl37", "Rpl37a", "Rpl38", "Rpl39", "Rpl4", 
                  "Rpl40", "Rpl41", "Rpl44", "Rpl5", "Rpl6", "Rpl7", "Rpl7a", "Rpl8", 
                  "Rpl9", "Rplp0", "Rplp1", "Rplp2", "Rps10", "Rps11", "Rps12", 
                  "Rps13", "Rps14", "Rps15", "Rps15a", "Rps16", "Rps17", "Rps18", 
                  "Rps19", "Rps2", "Rps20", "Rps21", "Rps23", "Rps24", "Rps25", 
                  "Rps26", "Rps27", "Rps27a", "Rps28", "Rps29", "Rps3", "Fau", 
                  "Rps3a", "Rps4x", "Rps5", "Rps6", "Rps7", "Rps8", "Rps9", "Rpsa"
  )
  
  # Generated using:
  # ~/localwork/facility/projects/250408SophieWiszniak-mouse_heart_Mib1_KD_scRNA/02-download_mouse_haemoglobin_genes.R
  hb_genes <- c("Hba-a1", "Hba-a2", "Hba-x", "Hbb-b1", "Hbb-b2", "Hbb-bh1", 
                "Hbb-bh2", "Hbb-bs", "Hbb-bt", "Hbb-y", "Hbq1a", "Hbq1b")
  
  if(!is.null(grouping_variable)){
    group_var_sym <- rlang::sym(grouping_variable)
  }
  
  n_dims <- as.integer(n_dims)
  
  basic_plot_dir <- file.path(plot_out_dir, "qc_plots-basic")
  clustered_plot_dir <- file.path(plot_out_dir, "qc_plot-clustered")
  dir.create(path = plot_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = basic_plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = clustered_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  
  ########################  
  ### Calculate QC metrics
  ########################
  
  # Percent mitochondrial
  message("Calculating percentage of mitochondrial reads")
  object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "^mt-")
  
  # Percent ribosomal
  message("Calculating percentage of ribosomal reads")
  object[["percent.ribo"]] <- 
    PercentageFeatureSet(object,
                         features = ribo_genes[which(ribo_genes %in% Features(object))])
  
  # Percent haemoglobin
  message("Calculating percentage of heomoglobin reads")
  object[["percent.hb"]] <-
    PercentageFeatureSet(object,
                         features = hb_genes[which(hb_genes %in% Features(object))])
  
  # Percent most abundant gene
  message("Determining most abundant gene per cell")
  max_expression <-
    LayerData(object = object, layer = "counts") %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "barcode", values_to = "counts") %>%
    group_by(barcode) %>%
    slice_max(counts, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::rename(max_counts = counts, # I have assumed Nick wants dplyr
           max_gene = gene)
  
  object <- add_metadata(object = object, metadata = max_expression)
  
  object[[]] <-
    object[[]] %>%
    mutate(percent.most_abundant = max_counts / nCount_RNA * 100 )
  
  
  # Transcriptome complexity
  # I won't model the expected complexity (nFeature_RNA ~ nCount_RNA) as it isn't linear and it is skewed by outliers
  # It's close to linear once logged and with outliers removed, but outlier removal (eg. < 100 genes) isn't something I want to implement automatically
  # I'll just make violin scatter plots
  message("Calculating transcriptome complexity")
  object[[]] <-
    object[[]] %>%
    mutate(complexity = nFeature_RNA / nCount_RNA,
           log10complexity = log10(complexity + 1))
  # ^ I've calculated this differently to https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html#Plotting_complexity
  
  
  
  ####################
  ### Initial QC plots
  ####################
  
  # If a grouping variable is provided, just plot the number of cells in each group
  if(!is.null(grouping_variable)){
    p <-
      object[[]] %>%
      group_by(!!group_var_sym) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      ggplot(aes(x = !!group_var_sym, y = n)) +
      geom_col() +
      # labs(subtitle = "Top 15 most abundant genes by frequency") +
      theme_bw() +
      guides(x = guide_axis(angle = 65)) +
      theme(axis.title.x = element_blank())
    print(p)
    ggsave(filename = file.path(basic_plot_dir, "01-Bar-cells_per_group.png"),
           plot = p)
  }
  
  # Bar plot of max gene (how many cells for each max gene)
  message("Producing plot of top most abundant genes")
  p <-
    object[[]] %>%
    group_by(max_gene) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    slice_max(n, n = 20) %>%
    ggplot(aes(x = reorder(max_gene, dplyr::desc(n)), y = n)) +
    geom_col() +
    labs(subtitle = "Top 15 most abundant genes by frequency") +
    theme_bw() +
    guides(x = guide_axis(angle = 65)) +
    theme(axis.title.x = element_blank())
  print(p)
  ggsave(filename = file.path(basic_plot_dir, "01-Bar-most_abundant_gene.png"),
         plot = p)
  
  
  # Violin plot: QC metrics - ungrouped 
  # message("Producing ungrouped QC plots")
  # p <- 
  #   object[[]] %>%
  #   as_tibble(rownames = "barcode") %>%
  #   select(barcode, all_of(qc_metrics)) %>%
  #   pivot_longer(-barcode, names_to = "metric", values_to = "value") %>%
  #   ggplot(aes(x = metric, y = value)) +
  #   gghalves::geom_half_violin(side = "r",
  #                              alpha = 0.5,
  #                              fill = "dodgerblue",
  #                              draw_quantiles = 0.5,
  #                              scale = "width") +
  #   gghalves::geom_half_point(side = "l", alpha = 0.3) +
  #   facet_wrap(~ metric, scales = "free") +
  #   theme_bw()
  # # Used suppressMessages to hide "Picking joint bandwidth of 0.012" for every faceted plot
  # suppressMessages(print(p))
  # suppressMessages(
  #   ggsave(filename = file.path(basic_plot_dir, "02-Violin-ungrouped_QC_metrics.png"),
  #          plot = p,
  #          width = 11, height = 9)
  # )
  # Grok tell me there is a compatability issue with ggplot 4.0.0 and gghalves 0.1.4
  # this is a known issue and the recommended fix is to move to ggdist (the new and better gghalves)
  message("Producing ungrouped QC plots")
  
  p <- object[[]] %>%
    as_tibble(rownames = "barcode") %>%
    select(barcode, all_of(qc_metrics)) %>%
    pivot_longer(-barcode, names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = metric, y = value)) +
    
    # Right half: density + median line
    stat_halfeye(
      side = "right",
      fill = "dodgerblue",
      alpha = 0.5,
      .width = 0.5,           # median line
      justification = -0.2,
      scale = 0.9,
      linewidth = 0.4,
      color = "black"
    ) +
    
    # Left half: dots
    stat_dots(
      side = "left",
      alpha = 0.3,
      justification = 1.25,
      quantiles = NULL,       # prevents binwidth error
      dotsize = 0.7
    ) +
    
    facet_wrap(~ metric, scales = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  print(p)
  ggsave(
    filename = file.path(basic_plot_dir, "02-Violin-ungrouped_QC_metrics.png"),
    plot = p,
    width = 11, height = 9, dpi = 300
  )
  
  
  
  # Violin plot: QC metrics - grouped by user-defined variable
  # if (!is.null(grouping_variable)) {
  #   message("Producing QC plots grouped by ", grouping_variable)
  #   
  #   p <- 
  #     object[[]] %>%
  #     as_tibble(rownames = "barcode") %>%
  #     select(barcode, !!group_var_sym, all_of(qc_metrics)) %>%
  #     pivot_longer(-c(barcode, !!group_var_sym), names_to = "metric", values_to = "values") %>%
  #     ggplot(aes(x = !!group_var_sym, y = values)) +
  #     gghalves::geom_half_violin(side = "r",
  #                                alpha = 0.5,
  #                                fill = "dodgerblue",
  #                                draw_quantiles = 0.5,
  #                                scale = "width") +
  #     gghalves::geom_half_point(side = "l", alpha = 0.3) +
  #     theme_bw() + facet_wrap(~ metric, scales = "free_y") +
  #     guides(x = guide_axis(angle = 65))
  #   
  #   suppressMessages(print(p))
  #   ggsave(filename = file.path(basic_plot_dir, paste0("03-Violin-QC_metric_by_", grouping_variable, ".png")),
  #          plot = p,
  #          width = 13, height = 9)
  # } else {
  #   message("No grouping variable provided. Skipping grouped plot.")
  # }
  if (!is.null(grouping_variable)) {
    message("Producing QC plots grouped by ", grouping_variable)
    
    group_var_sym <- sym(grouping_variable)
    
    p <- object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, !!group_var_sym, all_of(qc_metrics)) %>%
      pivot_longer(
        cols = -c(barcode, !!group_var_sym),
        names_to = "metric",
        values_to = "values"
      ) %>%
      ggplot(aes(x = !!group_var_sym, y = values)) +
      
      # Right half: density + median
      stat_halfeye(
        side = "right",
        fill = "dodgerblue",
        alpha = 0.5,
        .width = 0.5,
        justification = -0.2,
        scale = 0.9,
        linewidth = 0.4,
        color = "black"
      ) +
      
      # Left half: dots (fixed for categorical x)
      stat_dots(
        side = "left",
        alpha = 0.3,
        justification = 1.25,
        quantiles = NULL,       # THIS IS THE KEY FIX
        dotsize = 0.7
      ) +
      
      facet_wrap(~ metric, scales = "free_y") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 65, hjust = 1),
        panel.grid.major.x = element_blank()
      ) +
      labs(x = grouping_variable, y = "Value")
    
    print(p)
    ggsave(
      filename = file.path(basic_plot_dir, paste0("03-Violin-QC_metric_by_", grouping_variable, ".png")),
      plot = p,
      width = 13, height = 9, dpi = 300
    )
    
  } else {
    message("No grouping variable provided. Skipping grouped plot.")
  }
    
 
  ####################################################
  ### Normalisation, dimensional reduction, clustering
  ####################################################
  
  if(dim_reduction){
    # Run a non-optimised data normalisation, scaling, dimension reduction pipeline
    message("Running normalisation and dimensional reduction algorithms")
    message("- Clustering resolution: ", clustering_resolution)
    message("- Using ", n_dims, " dimensions for neighbourhood and UMAP analyses")
    object <-
      object %>%
      NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>%
      ScaleData(features = rownames(object), verbose = FALSE) %>%
      RunPCA(verbose = FALSE) %>% # Defaults to features = VariableFeatures(object)), but better not to specify as these are not yet stored in that object
      FindNeighbors(dims = 1:n_dims, verbose = FALSE) %>%
      FindClusters(resolution = clustering_resolution, verbose = FALSE) %>%
      RunUMAP(dims = 1:n_dims, verbose = FALSE)
    
    n_clusters <- object[[]] %>% pull(seurat_clusters) %>% unique() %>% length()
    message("- Found ", n_clusters, " clusters")
    
    
    
    ####################
    ### Score cell cycle
    ####################
    message("Calculating cell cycle scores")
    object <-
      CellCycleScoring(object = object,
                       s.features = s.genes,
                       g2m.features = g2m.genes)
    
    
    
    #######################
    ### Dim reduction plots
    #######################

    # PCA plot by user-specified group
    if(!is.null(grouping_variable)){
      message("Plotting ", grouping_variable, " onto PCA")
      p <-
        object %>%
        DimPlot(reduction = "pca", group.by = grouping_variable)
      print(p)
      ggsave(filename = file.path(clustered_plot_dir, paste0("01-PCA-", grouping_variable, ".png")),
             plot = p)
    }

    # PCA plot showing QC metrics
    message("Plotting QC metrics onto PCA plot")
    p <-
      object %>%
      FeaturePlot(reduction = "pca",
                  features = all_qc_metrics) %>%
      lapply(function(x) {x + scale_colour_distiller(palette = "Spectral", direction = -1)}) %>%
      patchwork::wrap_plots(nrow = 3, guides = "keep")
    
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "02-PCA-QC_metrics.png"),
           plot = p, scale = 2)

    # PCA plot showing clusters
    message("Plotting clusters onto PCA")
    p <-
      object %>%
      DimPlot(reduction = "pca", label = TRUE, label.box = TRUE, repel = TRUE)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "03-PCA-clusters.png"),
           plot = p)


    # UMAP by user-specified group
    if(!is.null(grouping_variable)){
      message("Plotting ", grouping_variable, " onto UMAP")

      plots <-
        lapply(object[[]] %>% pull(group_var_sym) %>% unique(), function(x){
          p <-
            object@reductions$umap@cell.embeddings %>%
            as_tibble(rownames = "barcode") %>%
            left_join(y = object[[]] %>% as_tibble(rownames = "barcode"),
                      by = "barcode") %>%
            ggplot(aes(x = umap_1,
                       y = umap_2)) +
            # Plot the background points
            geom_point(data = . %>% filter(!!group_var_sym != x),
                       size = 0.2, colour = "grey", alpha = 0.3) +
            # Plot the points of interest
            geom_point(data = . %>% filter(!!group_var_sym == x),
                       size = 0.2, colour = "dodgerblue", alpha = 0.6) +
            labs(subtitle = x) +
            theme_bw() +
            theme(plot.subtitle = element_text(colour = "dodgerblue", face = "bold"))
        }) %>%
        patchwork::wrap_plots()

      print(plots)
      ggsave(filename = file.path(clustered_plot_dir, paste0("04-UMAP-", grouping_variable, ".png")),
             plot = plots,
             width = 15, height = 9)
    }


    # UMAP showing clusters
    message("Plotting clusters onto UMAP")
    p <-
      object %>%
      DimPlot(reduction = "umap", label = TRUE, label.box = TRUE, repel = TRUE)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "05-UMAP-clusters.png"),
           plot = p)


    # UMAP showing QC metrics
    message("Plotting QC metrics onto UMAP")
    p <-
      object %>%
      FeaturePlot(reduction = "umap",
                  features = all_qc_metrics) %>%
      lapply(function(x) {x + scale_colour_distiller(palette = "Spectral", direction = -1)}) %>%
      patchwork::wrap_plots(nrow = 3, guides = "keep")
    
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "06-UMAP-QC_metrics.png"),
           plot = p,
           # scale = 2,
           width = 15, height = 11)


  #   # Violin plot: QC metrics by cluster
  #   message("Producing QC plots grouped by cluster")
  #   p <-
  #     object[[]] %>%
  #     as_tibble(rownames = "barcode") %>%
  #     select(barcode, seurat_clusters, all_of(all_qc_metrics)) %>%
  #     pivot_longer(-c(barcode, seurat_clusters), names_to = "metric", values_to = "values") %>%
  #     ggplot(aes(x = seurat_clusters, y = values)) +
  #     gghalves::geom_half_violin(side = "r",
  #                                alpha = 0.5,
  #                                fill = "dodgerblue",
  #                                draw_quantiles = 0.5,
  #                                scale = "width") +
  #     gghalves::geom_half_point(side = "l", alpha = 0.3) +
  #     theme_bw() + facet_wrap(~ metric, scales = "free_y") +
  #     guides(x = guide_axis(angle = 65))
  #   suppressMessages(print(p))
  #   suppressMessages(
  #     ggsave(filename = file.path(clustered_plot_dir, "07-Violin-QC_metrics_by_cluster.png"),
  #            plot = p,
  #            width = 15, height = 9)
  #   )
    
    # Violin plot: QC metrics by cluster
    message("Producing QC plots grouped by cluster")
    
    p <- object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, seurat_clusters, all_of(all_qc_metrics)) %>%
      pivot_longer(
        cols = -c(barcode, seurat_clusters),
        names_to = "metric",
        values_to = "values"
      ) %>%
      ggplot(aes(x = seurat_clusters, y = values)) +
      
      # Right half: smooth density + median line
      stat_halfeye(
        side = "right",
        fill = "dodgerblue",
        alpha = 0.5,
        .width = 0.5,           # draws the median (same as draw_quantiles = 0.5)
        justification = -0.2,
        scale = 0.9,
        linewidth = 0.4,
        color = "black"
      ) +
      
      # Left half: dots — fixed for categorical x (clusters)
      stat_dots(
        side = "left",
        alpha = 0.3,
        justification = 1.25,
        quantiles = NULL,       # ← THIS PREVENTS THE binwidth / xscale ERROR
        dotsize = 0.7
      ) +
      
      facet_wrap(~ metric, scales = "free_y") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 65, hjust = 1),
        panel.grid.major.x = element_blank()
      ) +
      labs(x = "Seurat Cluster", y = "Value")
    
    # Print and save cleanly
    print(p)
    ggsave(
      filename = file.path(clustered_plot_dir, "07-Violin-QC_metrics_by_cluster.png"),
      plot = p,
      width = 15, height = 9, dpi = 300
    )
  }
  
  return(object)
}




# Function to create QC metrics for FILTERED data, then plot them -------------------------------------------------

filtered_qc <- function(object, clustering_resolution, grouping_variable = NULL, plot_out_dir){
  
  ##########
  ### Checks
  ##########
  
  assertthat::assert_that(class(object) == "Seurat",
                          msg = "Object must be a Seurat object")
  
  assertthat::assert_that(is.numeric(clustering_resolution),
                          msg = "Clustering resolution must be numeric")
  
  if(!is.null(grouping_variable)){
    assertthat::assert_that(assertthat::is.string(grouping_variable),
                            grouping_variable %in% colnames(object[[]]),
                            msg = "Grouping variable must be a string, and be present in the Seurat object metadata") 
  }
  
  assertthat::assert_that(!is.null(object@reductions),
                          msg = "Must run RunUMAP() prior to running this function")
  
  assertthat::assert_that(max(grepl(pattern = "cluster", x = object[[]] %>% colnames())) == 1,
                          msg = "Seurat object must contain cluster information prior to running this function")
  
  
  
  #########
  ### Setup
  #########
  
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.most_abundant", "complexity"#,
                  # "log10complexity"
  )
  cell_cycle_metrics <- c("S.Score", "G2M.Score")
  all_qc_metrics <- c(qc_metrics, cell_cycle_metrics)
  
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  
  
  if(!is.null(grouping_variable)){
    group_var_sym <- rlang::sym(grouping_variable)
  }
  
  
  clustered_plot_dir <- file.path(plot_out_dir, "qc_plot-clustered")
  dir.create(path = plot_out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path = clustered_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  
  
  ####################
  ### Initial QC plots
  ####################
  
  # If a grouping variable is provided, just plot the number of cells in each group
  if(!is.null(grouping_variable)){
    p <-
      object[[]] %>%
      group_by(!!group_var_sym) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      ggplot(aes(x = !!group_var_sym, y = n)) +
      geom_col() +
      theme_bw() +
      guides(x = guide_axis(angle = 65)) +
      theme(axis.title.x = element_blank())
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "01-Bar-cells_per_group.png"),
           plot = p)
  }
  
  
  # Violin plot: QC metrics - grouped by user-defined variable
  if (!is.null(grouping_variable)) {
    message("Producing QC plots grouped by ", grouping_variable)
    
    p <- 
      object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, !!group_var_sym, all_of(all_qc_metrics)) %>%
      pivot_longer(-c(barcode, !!group_var_sym), names_to = "metric", values_to = "values") %>%
      ggplot(aes(x = !!group_var_sym, y = values)) +
      gghalves::geom_half_violin(side = "r",
                                 alpha = 0.5,
                                 fill = "dodgerblue",
                                 draw_quantiles = 0.5,
                                 scale = "width") +
      gghalves::geom_half_point(side = "l", alpha = 0.3) +
      theme_bw() + facet_wrap(~ metric, scales = "free_y") +
      guides(x = guide_axis(angle = 65))
    
    suppressMessages(print(p))
    ggsave(filename = file.path(clustered_plot_dir, paste0("02-Violin-QC_metric_by_", grouping_variable, ".png")),
           plot = p,
           width = 13, height = 9)
  } else {
    message("No grouping variable provided. Skipping grouped plot.")
  }
  
  
  
  #######################
  ### Dim reduction plots
  #######################
  
  # PCA plot by user-specified group
  if(!is.null(grouping_variable)){
    message("Plotting ", grouping_variable, " onto PCA")
    p <-
      object %>%
      DimPlot(reduction = "pca", group.by = grouping_variable)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, paste0("03-PCA-", grouping_variable, ".png")),
           plot = p)
  }
  
  # PCA plot showing QC metrics
  message("Plotting QC metrics onto PCA plot")
  p <-
    object %>%
    FeaturePlot(reduction = "pca",
                features = all_qc_metrics)
  print(p)
  ggsave(filename = file.path(clustered_plot_dir, "04-PCA-QC_metrics.png"),
         plot = p, scale = 2)
  
  # PCA plot showing clusters
  message("Plotting clusters onto PCA")
  p <-
    object %>%
    DimPlot(reduction = "pca", label = TRUE, label.box = TRUE, repel = TRUE)
  print(p)
  ggsave(filename = file.path(clustered_plot_dir, "05-PCA-clusters.png"),
         plot = p)
  
  
  # UMAP by user-specified group
  if(!is.null(grouping_variable)){
    message("Plotting ", grouping_variable, " onto UMAP")
    
    plots <-
      lapply(object[[]] %>% pull(group_var_sym) %>% unique(), function(x){
        p <-
          object@reductions$umap@cell.embeddings %>%
          as_tibble(rownames = "barcode") %>%
          left_join(y = object[[]] %>% as_tibble(rownames = "barcode"),
                    by = "barcode") %>%
          mutate(assignment = ifelse(!!group_var_sym == x, !!group_var_sym, "other")) %>%
          ggplot(aes(x = umap_1,
                     y = umap_2,
                     colour = assignment,
                     alpha = assignment)) +
          geom_point(size = 0.2) +
          scale_alpha_manual(values = c(0.4, 0.9)) +
          scale_colour_manual(values = c("grey", "darkorange2")) +
          labs(subtitle = x) +
          theme_bw() +
          theme(plot.subtitle = element_text(colour = "darkorange2", face = "bold")) +
          NoLegend()
      }) %>%
      patchwork::wrap_plots()
    
    print(plots)
    ggsave(filename = file.path(clustered_plot_dir, paste0("06-UMAP-", grouping_variable, ".png")),
           plot = plots,
           width = 15, height = 9)
    
    # UMAP showing clusters
    message("Plotting clusters onto UMAP")
    p <-
      object %>%
      DimPlot(reduction = "umap", label = TRUE, label.box = TRUE, repel = TRUE)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "07-UMAP-clusters.png"),
           plot = p)
    
    # UMAP
    message("Plotting QC metrics onto UMAP")
    p <-
      object %>%
      FeaturePlot(reduction = "umap",
                  features = all_qc_metrics)
    print(p)
    ggsave(filename = file.path(clustered_plot_dir, "08-UMAP-QC_metrics.png"),
           plot = p,
           # scale = 2,
           width = 15, height = 11)
    
    
    # Violin plot: QC metrics - grouped by cluster
    message("Producing QC plots grouped by cluster")
    p <-
      object[[]] %>%
      as_tibble(rownames = "barcode") %>%
      select(barcode, seurat_clusters, all_of(all_qc_metrics)) %>%
      pivot_longer(-c(barcode, seurat_clusters), names_to = "metric", values_to = "values") %>%
      ggplot(aes(x = seurat_clusters, y = values)) +
      gghalves::geom_half_violin(side = "r",
                                 alpha = 0.5,
                                 fill = "dodgerblue",
                                 draw_quantiles = 0.5,
                                 scale = "width") +
      gghalves::geom_half_point(side = "l", alpha = 0.3) +
      theme_bw() + facet_wrap(~ metric, scales = "free_y") +
      guides(x = guide_axis(angle = 65))
    suppressMessages(print(p))
    suppressMessages(
      ggsave(filename = file.path(clustered_plot_dir, "09-Violin-QC_metrics_by_cluster.png"),
             plot = p,
             width = 15, height = 9)
    )
  }
  
  # return(object)
  invisible(NULL)
}



# Violin/scatter plot -----------------------------------------------------

violin_points <- function(object, grouping_variable, variables_to_plot, fill = "dodgerblue", alpha = 0.5) {

  assertthat::assert_that(class(object) == "Seurat",
                          msg = "Object must be a Seurat object")
  
  if(!is.null(grouping_variable)){
    assertthat::assert_that(assertthat::is.string(grouping_variable),
                            grouping_variable %in% colnames(object[[]]),
                            msg = "Grouping variable must be a string, and be present in the Seurat object metadata")
    
    group_var_sym <- rlang::sym(grouping_variable)
  }
  
  assertthat::assert_that(is.character(variables_to_plot),
                          all(variables_to_plot %in% colnames(object[[]])),
                          msg = "variables_to_plot must be a character vector, and all variables must be present in the Seurat object metadata")
  
  assertthat::assert_that(substr(x = fill, start = 1, stop = 1) == "#" | fill %in% colours(),
                          msg = "Please supply valid colour from `colours()` or in hexadecimal format") 
  
  assertthat::assert_that(assertthat::is.number(alpha) & alpha <= 1,
                          msg = "Please supply a numeric value less than 1 for the `alpha` parameter")

  # object[[]] %>%
  #   as_tibble(rownames = "barcode") %>%
  #   select(barcode, !!group_var_sym, all_of(variables_to_plot)) %>%
  #   pivot_longer(-c(barcode, !!group_var_sym), names_to = "metric", values_to = "values") %>%
  #   ggplot(aes(x = !!group_var_sym, y = values)) +
  #   gghalves::geom_half_violin(side = "r",
  #                              alpha = alpha,
  #                              fill = fill,
  #                              draw_quantiles = 0.5,
  #                              scale = "width") +
  #   gghalves::geom_half_point(side = "l", alpha = 0.3) +
  #   theme_bw() +
  #   facet_wrap(~ metric, scales = "free_y") +
  #   guides(x = guide_axis(angle = 65))
  
  object[[]] %>%
    as_tibble(rownames = "barcode") %>%
    select(barcode, !!group_var_sym, all_of(variables_to_plot)) %>%
    pivot_longer(
      cols = -c(barcode, !!group_var_sym),
      names_to = "metric",
      values_to = "values"
    ) %>%
    ggplot(aes(x = !!group_var_sym, y = values)) +
    
    # Right half: smooth density + median line
    stat_halfeye(
      side = "right",
      fill = fill,
      alpha = alpha,
      .width = 0.5,                    # draws median (same as draw_quantiles = 0.5)
      justification = -0.2,
      scale = 0.9,
      linewidth = 0.4,
      color = "black"
    ) +
    
    # Left half: dots — FIXED for categorical x-axis
    stat_dots(
      side = "left",
      alpha = 0.3,
      justification = 1.25,
      quantiles = NULL,                # THIS PREVENTS THE binwidth/xscale ERROR
      dotsize = 0.7
    ) +
    
    facet_wrap(~ metric, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 65, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    labs(x = as_label(group_var_sym))
}



# Dimensional reduction by user-specified group ---------------------------

# # UMAP by user-specified group
# if(!is.null(grouping_variable)){
#   message("Plotting ", grouping_variable, " onto UMAP")
#   
#   plots <-
#     lapply(object[[]] %>% pull(group_var_sym) %>% unique(), function(x){
#       p <-
#         object@reductions$umap@cell.embeddings %>%
#         as_tibble(rownames = "barcode") %>%
#         left_join(y = object[[]] %>% as_tibble(rownames = "barcode"),
#                   by = "barcode") %>%
#         ggplot(aes(x = umap_1,
#                    y = umap_2)) +
#         # Plot the background points
#         geom_point(data = . %>% filter(!!group_var_sym != x),
#                    size = 0.2, colour = "grey", alpha = 0.3) +
#         # Plot the points of interest
#         geom_point(data = . %>% filter(!!group_var_sym == x),
#                    size = 0.2, colour = "dodgerblue", alpha = 0.6) +
#         labs(subtitle = x) +
#         theme_bw() +
#         theme(plot.subtitle = element_text(colour = "dodgerblue", face = "bold"))
#     }) %>%
#     patchwork::wrap_plots()
#   
#   print(plots)
#   ggsave(filename = file.path(clustered_plot_dir, paste0("04-UMAP-", grouping_variable, ".png")),
#          plot = plots,
#          width = 15, height = 9)
# }
