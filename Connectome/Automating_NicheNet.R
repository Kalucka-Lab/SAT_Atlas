######################################################################
## Title        : Automating NicheNet analysis
##
## Description  :
##
## Author       : "Mohamed Hassan"
##
## Date         : "`r format(Sys.time(),  '%d %B %Y')`"
######################################################################
##
## Loading required libraries:########################################
pacman::p_load(tidyverse, Seurat, patchwork, remotes, nichenetr, cowplot, ggpubr)

# Load data
seuratObj <- readRDS("R_Projects/scAtlas/Analysis/BigObj/Objects/Atlas_IntegratedObj_with_all_sublabels_transferred.RDS")

# Remove low-quality clusters
clusters2keep <- grep(pattern = "Low-quality", x = (seuratObj$all_sub_labels %>% factor() %>% levels()), value = TRUE, invert = TRUE)
seuratObj <- subset(seuratObj, all_sub_labels %in% clusters2keep)

# Load NicheNet networks
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds")) %>% distinct(from, to)
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# Define Scenarios
scenarios <- list(
  "Scenario_1" = list(sender = c("Activated Vascular Endothelial cells"),
                      receivers = c("Adipocyte_3", "Adipocyte_5",
                                    "Adipocyte-like Endothelial cells",
                                    "ASPCs_4", "ASPCs_5",
                                    "FAP-like Endothelial cells",
                                    # "Lymphatic Endothelial cells",
                                    "Lymphoid-like Endothelial cells",
                                    # "Myeloid-like Endothelial cells 2",
                                    )),
  "Scenario_2" = list(sender = c("FAP-like Endothelial cells"),
                      receivers = c("Adipocyte_3",
                                    "Adipocyte-like Endothelial cells",
                                    "ASPCs_4", "ASPCs_5", "ASPCs_6",
                                    "Lymphoid-like Endothelial cells",
                                    "Lymphatic Endothelial cells",
                                    "Myeloid-like Endothelial cells 2",
                                    "Activated Vascular Endothelial cells")),
  "Scenario_3" = list(sender = c("Lymphoid-like Endothelial cells"),
                      receivers = c("Adipocyte_3",
                                    "Adipocyte-like Endothelial cells",
                                    "ASPCs_4", "ASPCs_5",
                                    "FAP-like Endothelial cells",
                                    "Lymphatic Endothelial cells",
                                    "Myeloid-like Endothelial cells 2"))
)

gc()

# Define conditions
condition_oi <- "Diabetic Obese"
condition_reference <- "Lean"

# Initialize an empty list to store plots
list_of_combined_plots <- list()


# Loop through scenarios and process each receiver separately
for (scenario in names(scenarios)) {
  sender_celltypes <- scenarios[[scenario]]$sender
  receivers <- scenarios[[scenario]]$receivers
  
  for (receiver in receivers) {
    print(paste("Processing:", scenario, "Sender:", sender_celltypes, "Receiver:", receiver))
    
    # Create unique output directory
    output_dir <- paste0("R_Projects/scAtlas/Analysis/BigObj/NicheNet/Loop/Results/", scenario, "/", sender_celltypes, "/", receiver, "/")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Define expressed genes
    Idents(seuratObj) <- seuratObj$all_sub_labels
    expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)
    expressed_genes_sender <- get_expressed_genes(sender_celltypes, seuratObj, pct = 0.05)
    
    # Define potential ligands
    all_receptors <- unique(lr_network$to)
    expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
    potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
    potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)
    
    # Differential expression analysis
    
    # Subset the receiver cell type
    seurat_obj_receiver <- subset(seuratObj, idents = receiver)
    
    # **Skip if the receiver has too few cells**
    if (ncol(seurat_obj_receiver) < 10) { 
      print(paste("Skipping:", receiver, "- Too few cells present."))
      next
    }
    
    # Differential expression
    DE_table_receiver <- FindMarkers(object = seurat_obj_receiver,
                                     ident.1 = condition_oi, ident.2 = condition_reference,
                                     group.by = "Condition",
                                     latent.vars = c("gender", "Study_chemistry"),
                                     test.use = "LR",
                                     min.pct = 0.05) %>% rownames_to_column("gene")
    
    # Define gene set of interest
    geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    # **Skip if the gene set is empty**
    if (length(geneset_oi) == 0) {
      print(paste("Skipping:", receiver, "- No significant DE genes found."))
      next
    }
    
    # **Skip if all genes have the same response (all upregulated or all downregulated)**
    if (length(unique(sign(DE_table_receiver$avg_log2FC))) == 1) {
      print(paste("Skipping:", receiver, "- All genes have the same response (no variation)."))
      next
    }
    
    
    
    # Define background genes
    background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    length(background_expressed_genes)
    length(geneset_oi)
    
    # Perform NicheNet ligand activity analysis
    ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                   background_expressed_genes = background_expressed_genes,
                                                   ligand_target_matrix = ligand_target_matrix,
                                                   potential_ligands = potential_ligands)
    ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
    ligand_activities
    
    # Save ligand activities table
    write.csv(ligand_activities, paste0(output_dir, "ligand_activities.csv"), row.names = FALSE)
    
    # Save ligand activity histogram
    p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
      geom_histogram(color="black", fill="darkorange")  + 
      geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
                 color="red", linetype="dashed", size=1) + 
      labs(x="Ligand Activity (PCC)", y="# Ligands") +
      theme_classic()
    
    ggsave(paste0(output_dir, "ligand_activity_histogram.png"), p_hist_lig_activity, width = 6, height = 4)
    
    best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% 
      arrange(-aupr_corrected) %>% pull(test_ligand)
    
    vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
      column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
    
    (make_heatmap_ggplot(vis_ligand_aupr,
                         "Prioritized ligands", "Ligand activity", 
                         legend_title = "AUPR", color = "darkorange") + 
        theme(axis.text.x.top = element_blank()))  
    
    # Infer target genes and receptors of top-ranked ligands
    
    active_ligand_target_links_df <- best_upstream_ligands %>%
      lapply(get_weighted_ligand_target_links,
             geneset = geneset_oi,
             ligand_target_matrix = ligand_target_matrix,
             n = 100) %>%
      bind_rows() %>% drop_na()
    
    nrow(active_ligand_target_links_df)
    head(active_ligand_target_links_df)
    
    active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                      ligand_target_matrix = ligand_target_matrix,
                                                                      cutoff = 0.33)
    nrow(active_ligand_target_links)
    head(active_ligand_target_links)
    
    order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
    order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
    
    vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
    
    make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                        color = "purple", legend_title = "Regulatory potential") +
      scale_fill_gradient2(low = "whitesmoke",  high = "purple")
    
    
    
    # Receptors of top-ranked ligands
    # Receptor plot
    ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
      best_upstream_ligands, expressed_receptors,
      lr_network, weighted_networks$lr_sig) 
    
    vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
      ligand_receptor_links_df,
      best_upstream_ligands,
      order_hclust = "both") 
    
    
    p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                             y_name = "Ligands", x_name = "Receptors",  
                                             color = "mediumvioletred", legend_title = "Prior interaction potential")
    
    p_ligand_receptor
    
    # Sender-focused approach
    
    ligand_activities_all <- ligand_activities
    best_upstream_ligands_all <- best_upstream_ligands
    
    ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
    best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
    
    
    ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>% column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
    vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 
    
    p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                         "Prioritized ligands", "Ligand activity", 
                                         legend_title = "AUPR", color = "darkorange") + 
      theme(axis.text.x.top = element_blank())
    
    p_ligand_aupr
    
    
    
    
    # Target gene plot
    active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix,  n = 100) %>% bind_rows() %>% drop_na()
    
    active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix,cutoff = 0.33) 
    
    order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
    order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
    
    vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
    p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                           color = "purple", legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
    
    p_ligand_target
    
    
    # Receptor plot
    ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig)
    
    vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df, best_upstream_ligands, order_hclust = "both")
    
    p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                             y_name = "Ligands", x_name = "Receptors",  
                                             color = "mediumvioletred", legend_title = "Prior interaction potential")
    
    p_ligand_receptor
    
    
    
    best_upstream_ligands_all %in% rownames(seuratObj) %>% table()
    
    
    # Dotplot of sender-focused approach
    p_dotplot <- DotPlot(subset(seuratObj, all_sub_labels %in% sender_celltypes),
                         features = rev(best_upstream_ligands), cols = "RdYlBu") + 
      coord_flip() +
      scale_y_discrete(position = "right")
    
    p_dotplot
    
    celltype_order <- levels(Idents(seuratObj)) 
    
    DE_table_top_ligands <- lapply(
      celltype_order[celltype_order %in% sender_celltypes],
      get_lfc_celltype, 
      seurat_obj = seuratObj,
      condition_colname = "Condition",
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      celltype_col = "all_sub_labels",
      min.pct = 0, logfc.threshold = 0,
      features = best_upstream_ligands 
    ) 
    
    DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
      column_to_rownames("gene") 
    
    vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE]) 
    
    p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                            "Prioritized ligands", "LFC in Sender",
                                            low_color = "midnightblue", mid_color = "white",
                                            mid = median(vis_ligand_lfc), high_color = "red",
                                            legend_title = "LFC")
    
    p_lfc
    
    
    (make_line_plot(ligand_activities = ligand_activities_all,
                    potential_ligands = potential_ligands_focused) +
        theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))
    
    # Generate summary visualization
    figures_without_legend <- cowplot::plot_grid(
      p_ligand_aupr + theme(legend.position = "none"),
      p_dotplot + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0)),
      p_lfc + theme(legend.position = "none"),
      p_ligand_target + theme(legend.position = "none"),
      align = "hv",
      nrow = 1)
    
    legends <- cowplot::plot_grid(
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
      nrow = 1)
    
    combined_plot <- cowplot::plot_grid(figures_without_legend, 
                                        legends, rel_heights = c(10,5), 
                                        nrow = 2, align = "hv")
    
    # Define sender and receiver names dynamically
    scenario_name <- paste0(gsub(" ", "_", sender_celltypes), 
                            "_to_", gsub(" ", "_", receiver))
    
    # Store the plot in the list with a meaningful name
    list_of_combined_plots[[scenario_name]] <- combined_plot
    
    ggsave(paste0(output_dir, "summary_plot.png"), 
           combined_plot, width = 30, 
           height = 10)
    ggsave(paste0(output_dir, "summary_plot.pdf"), combined_plot, 
           width = 30, height = 12, 
           device = cairo_pdf)
    pdf(paste0(output_dir, "summary_combined_plot.pdf"), 
        width = 30, height = 12)
    print(combined_plot)
    dev.off()
    print(paste("Finished:", scenario, "Sender:", sender_celltypes, "Receiver:", receiver))
    gc()
  }
}
