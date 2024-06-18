
Vascular_IntegratedObj <- readRDS("R_Projects/scAtlas/Analysis/Vascular/Objects/Vascular_IntegratedObj_label_transferred_from_SubECs_updated_names.RDS")
all_clusters <- Vascular_IntegratedObj$Vascular_labels %>% factor() %>% levels()
clusters2remove <- all_clusters[8:10]

## Creating a function by negation
"%out%" <- Negate("%in%")

Vascular_IntegratedObj <- subset(Vascular_IntegratedObj, technology == "Single Nucleus")
Vascular_IntegratedObj <- subset(Vascular_IntegratedObj, Vascular_labels %out% clusters2remove)


Vascular_IntegratedObj$Vascular_labels <- droplevels(Vascular_IntegratedObj$Vascular_labels)
Vascular_IntegratedObj$barcode <- rownames(Vascular_IntegratedObj@meta.data) 


Endo <- c("PECAM1", "CDH5", "VWF")
Mesenchymal <- c("ZEB2", "S100A4", "AIFM2", "COL1A1","COL3A1", 
                 "TWIST1", "ACTA2", "FAP", "SNAI1", "SNAI2", 
                 "SMAD2", "SMAD3", "SMAD4", "KDM4B", "MAPK8", "MAPK14",
                 "YAP1", "TGFB1","TGFB2","TGFB3","TGFBR1")


exp_mat <- as.matrix(Vascular_IntegratedObj[["RNA"]]@counts[c(columns_to_process),])
meta <- Vascular_IntegratedObj@meta.data %>% select(Vascular_labels)
meta <- dplyr::bind_cols(meta, exp_mat %>% t() %>% as.data.frame())
meta <- mutate(meta, barcode = rownames(meta))

columns_to_process <- c(Endo, Mesenchymal)


# Apply zscore to each column and create corresponding zscore columns
for (col in columns_to_process) {
  zscore_col <- paste0("Z_", col)
  meta[[zscore_col]] <- scale(meta[[col]])
}

# Calculate means for EMT and vascular markers
meta$mean_endo <- rowMeans(meta[Endo])
meta$mean_mesenchymal <- rowMeans(meta[Mesenchymal])
meta$diff <- meta$mean_mesenchymal - meta$mean_endo

main_df <- Vascular_IntegratedObj@meta.data

# meta <- meta[, c(26,51:53)]
meta <- meta %>% select(barcode, mean_endo, mean_mesenchymal, diff)

merged_df <- merge(main_df, meta, by= "barcode")
merged_df <- column_to_rownames(merged_df, var = "barcode")
Vascular_IntegratedObj@meta.data <- merged_df

pdf("R_Projects/scAtlas/Analysis/Vascular/Analysis/EndMT_score/ScatterPlot_Split_by_cluster.pdf", 
    width = 10, height = 10)
scCustomize::Split_FeatureScatter(Vascular_IntegratedObj,
                                  feature1 = "mean_endo",
                                  feature2 = "mean_mesenchymal",
                                  split.by = "Vascular_labels",
                                  group.by = "Vascular_labels") & NoLegend()
dev.off()



all_scatters <- scCustomize::Split_FeatureScatter(Vascular_IntegratedObj,
                                  feature1 = "mean_endo",
                                  feature2 = "mean_mesenchymal",
                                  split.by = "Vascular_labels",
                                  group.by = "Vascular_labels") & NoLegend()




pdf("SubECs_scatter.pdf", width = 10, height = 10);
for (i in all_scatters){
  print(i)
}
dev.off()


pdf("R_Projects/scAtlas/Analysis/Vascular/Analysis/EndMT_score/ScatterPlot_grouped_by_cluster_updated_colors_NoAxes.pdf", 
    width = 10, height = 10)
FeatureScatter(Vascular_IntegratedObj, feature1 = "mean_endo", feature2 = "mean_mesenchymal",
               group.by = "Vascular_labels", cols = colors2use) + NoLegend() + NoAxes()
dev.off()

ggplot(Vascular_IntegratedObj@meta.data, 
       aes(x = Vascular_labels, y = diff, fill = Vascular_labels)) + 
  geom_boxplot() + 
  theme_classic() + 
  NoLegend()



Vascular_IntegratedObj@meta.data$diff2 <- Vascular_IntegratedObj$mean_mesenchymal - Vascular_IntegratedObj$mean_endo


Vascular_IntegratedObj@meta.data %>% 
  select(Vascular_labels, mean_endo, mean_mesenchymal) %>% 
  as.data.frame() %>% 
  writexl::write_xlsx("R_Projects/scAtlas/Analysis/Vascular/Analysis/EndMT_score/MeanEndo_MeanMesenchymal_scores.xlsx")




ggplot(Vascular_IntegratedObj@meta.data, 
       aes(x=Vascular_labels,
           y=diff2_norm,
           fill=Vascular_labels)) +
  geom_violin(trim=FALSE,scale="width") +
  geom_boxplot(width=0.1, fill="white",outlier.shape = NA) +
  scale_fill_manual(values = colors2use) + 
  theme(axis.text.x = element_text(size = 0))


############################################################################
############################################################################
#########  Pseudotime cells with mean Endo and Mesenchymal scores ##########
############################################################################
############################################################################
  

Vascular_IntegratedObj <- readRDS("R_Projects/scAtlas/Analysis/SubECs/Analysis/Pseudotime/Scenario5/SeuratPseudotime.RDS")
Vascular_IntegratedObj$barcode <- rownames(Vascular_IntegratedObj@meta.data) 


Endo <- c("PECAM1", "CDH5", "VWF")
Mesenchymal <- c("ZEB2", "S100A4", "AIFM2", "COL1A1","COL3A1", 
                 "TWIST1", "ACTA2", "FAP", "SNAI1", "SNAI2", 
                 "SMAD2", "SMAD3", "SMAD4", "KDM4B", "MAPK8", "MAPK14",
                 "YAP1", "TGFB1","TGFB2","TGFB3","TGFBR1")


### For automation
# Columns to apply zscore and calculate means
# emt_columns <- c("PLIN1", "NEGR1", "PTPRC")
# vascular_columns <- c("CDH5", "PECAM1", "VWF")

columns_to_process <- c(Endo, Mesenchymal)



exp_mat <- as.matrix(Vascular_IntegratedObj[["RNA"]]@counts[c(columns_to_process),])
meta <- Vascular_IntegratedObj@meta.data %>% select(transferred_labels_from_Act, monocle3_pseudotime)
meta <- dplyr::bind_cols(meta, exp_mat %>% t() %>% as.data.frame())
meta <- mutate(meta, barcode = rownames(meta))

# Apply zscore to each column and create corresponding zscore columns
for (col in columns_to_process) {
  zscore_col <- paste0("Z_", col)
  meta[[zscore_col]] <- scale(meta[[col]])
}

# Calculate means for EMT and vascular markers
meta$mean_endo <- rowMeans(meta[Endo])
meta$mean_mesenchymal <- rowMeans(meta[Mesenchymal])
meta$diff <- meta$mean_mesenchymal - meta$mean_endo

main_df <- Vascular_IntegratedObj@meta.data

# meta <- meta[, c(26,51:53)]
meta <- meta %>% select(barcode, monocle3_pseudotime, mean_endo, mean_mesenchymal, diff)
# merged_df <- merge(main_df, meta, by= "barcode")
merged_df <- column_to_rownames(merged_df, var = "barcode")
Vascular_IntegratedObj@meta.data <- merged_df




ggplot(Vascular_IntegratedObj@meta.data, 
       aes(x=monocle3_pseudotime,
           y=diff, 
           group = transferred_labels_from_Act,
           fill = transferred_labels_from_Act)) +
  geom_boxplot() + theme_classic()
# scale_fill_manual(values = colors2use) + 
theme(axis.text.x = element_text(size = 0))


