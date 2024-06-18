## a) take top 100 GO terms based on q-value
## b) summarize GO terms so that they merge
## c) choose a name for the merged GO term
## d) plot common GO terms on a radar plot
## e) plot non common GO terms on a dot plot

library(GSEABase)
library(Seurat)
library(AUCell)
library(tidyverse)
library(readxl)
load("R_Projects/scAtlas/Analysis/Vascular/Analysis/DEA/LR/GO_terms/DEGs_Up/DEGs_Up_Vascular.Rdata")

##############
## subECs
##############

subECs_HO_vs_L <- df_Enrichment_Vascular_Clusters_DEG_HO_vs_L$`SubEndothelial cells`
subECs_DO_vs_L <- df_Enrichment_Vascular_Clusters_DEG_DO_vs_L$`SubEndothelial cells`
# subECs_DO_vs_HO <- df_Enrichment_Vascular_Clusters_DEG_DO_vs_HO$subECs

## arrange the df by q value and extract top 100 GO terms
subECs_HO_vs_L <- subECs_HO_vs_L %>% arrange(qvalue) %>% head(100)
subECs_DO_vs_L <- subECs_DO_vs_L %>% arrange(qvalue) %>% head(100)
# subECs_DO_vs_HO <- subECs_DO_vs_HO %>% arrange(qvalue) %>% head(100)


## Finding common and unique GO terms
# HO_vs_L_intersect_DO_vs_L <- intersect(subECs_HO_vs_L$ID, subECs_DO_vs_L$ID)
# HO_vs_L_intersect_DO_vs_HO <- intersect(subECs_HO_vs_L$ID, subECs_DO_vs_HO$ID)
# DO_vs_L_intersect_DO_vs_HO <- intersect(subECs_DO_vs_L$ID, subECs_DO_vs_HO$ID)

common_GO <- intersect(subECs_HO_vs_L$ID, subECs_DO_vs_L$ID)


####################
### Common GO terms
####################
# 1. HO vs L
common_GO_subECs_HO_vs_L <- subECs_HO_vs_L[(subECs_HO_vs_L$ID %in% common_GO),]

# 2. DO vs L
# common_GO_subECs_DO_vs_L <- subECs_DO_vs_L[(subECs_DO_vs_L$ID %in% common_GO),]

# 3. DO vs HO  
# common_GO_subECs_DO_vs_HO <- subECs_DO_vs_HO[(subECs_DO_vs_HO$ID %in% common_GO),]

###################################################################################################
###################################################################################################


## Summarising the common GO terms
common_GO_subECs_HO_vs_L_summarised <- 
  common_GO_subECs_HO_vs_L %>% 
  group_by(geneID) %>% 
  summarise(Description_updated = paste(Description, collapse = ', ')) %>% 
  dplyr::select(Description_updated, geneID) 



writexl::write_xlsx(common_GO_subECs_HO_vs_L_summarised, 
                    "subECs_common_GO_UP.xlsx")
# write.csv(uniqueGO_subECs_HO_vs_L_summarised, "subECs_uniqueGO_HO_vs_L.csv")
# write.csv(uniqueGO_subECs_DO_vs_L_summarised, "subECs_uniqueGO_DO_vs_L.csv")
# write.csv(uniqueGO_subECs_DO_vs_HO_summarised, "subECs_uniqueGO_DO_vs_HO.csv")

commonGO_subECs <- readxl::read_excel("Downloads/subECs_common_GO_UP.xlsx", 
                                         sheet = 1)
colnames(commonGO_subECs) <- commonGO_subECs[1,]
commonGO_subECs <- commonGO_subECs[-1,-1]
commonGO_subECs <- commonGO_subECs %>% filter(Filter  == "Yes")

gene_lists <- strsplit(commonGO_subECs$geneID, "/")
names(gene_lists) <- commonGO_subECs$Description_updated


###################################################################################################
###################################################################################################
## Reading the vascular object

Vascular_IntegratedObj <- readRDS("R_Projects/scAtlas/Analysis/Vascular/Objects/Vascular_IntegratedObj_label_transferred_from_subECs.RDS")

##########################################
## Subsetting only subECs cluster
##########################################
Vascular_IntegratedObj <- subset(Vascular_IntegratedObj, 
                                 Vascular_labels == "SubEndothelial cells")

## Subsets
Vascular_lean <- subset(Vascular_IntegratedObj, Condition == "Lean")
Vascular_obese <- subset(Vascular_IntegratedObj, Condition == "Obese")
Vascular_D_obese <- subset(Vascular_IntegratedObj, Condition == "Diabetic Obese")


## extract expression matrix
Vascular_lean_matrix <- GetAssayData(Vascular_lean, assay = "RNA", slot = "count")
Vascular_obese_matrix <- GetAssayData(Vascular_obese, assay = "RNA", slot = "count")
Vascular_D_obese_matrix <- GetAssayData(Vascular_D_obese, assay = "RNA", slot = "count")


library(Matrix)
Vascular_lean_matrix <- as(Vascular_lean_matrix, "dgCMatrix")
Vascular_obese_matrix <- as(Vascular_obese_matrix, "dgCMatrix")
Vascular_D_obese_matrix <- as(Vascular_D_obese_matrix, "dgCMatrix")



###################################################################################################
###################################################################################################

#################################
############### Lean ############ 
#################################


# Define a list of gene sets with corresponding names
# Create an empty data frame to store results
cells_AUC_lean <- data.frame()

# Iterate through the gene sets
for (i in names(gene_lists)) {
  # Extract the current gene set and name
  currentGeneSet <- gene_lists[[i]]
  currentSetName <- as.character(i)  # You can customize the name as needed
  
  # Create a GeneSet object with the current name
  geneSetObj <- GeneSet(currentGeneSet, setName = currentSetName)
  
  # Calculating AUC score for the current GeneSet in the D_obese expression matrix
  cells_AUC <- AUCell_run(Vascular_lean_matrix, geneSetObj)
  
  # Build rankings
  cells_rankings <- AUCell_buildRankings(Vascular_lean_matrix, plotStats = TRUE)
  
  # Calculate AUC
  cells_AUC <- AUCell_calcAUC(geneSetObj, cells_rankings)
  
  # Append the AUC values to the result data frame
  cells_AUC_lean <- rbind(cells_AUC_lean, cells_AUC@assays@data@listData[["AUC"]])
}

###################################################################################################
###################################################################################################

###################################
###############  Obese ############ 
###################################


# Define a list of gene sets with corresponding names
# Create an empty data frame to store results
cells_AUC_obese <- data.frame()

# Iterate through the gene sets
for (i in names(gene_lists)) {
  # Extract the current gene set and name
  currentGeneSet <- gene_lists[[i]]
  currentSetName <- as.character(i)  # You can customize the name as needed
  
  # Create a GeneSet object with the current name
  geneSetObj <- GeneSet(currentGeneSet, setName = currentSetName)
  
  # Calculating AUC score for the current GeneSet in the D_obese expression matrix
  cells_AUC <- AUCell_run(Vascular_obese_matrix, geneSetObj)
  
  # Build rankings
  cells_rankings <- AUCell_buildRankings(Vascular_obese_matrix, plotStats = TRUE)
  
  # Calculate AUC
  cells_AUC <- AUCell_calcAUC(geneSetObj, cells_rankings)
  
  # Append the AUC values to the result data frame
  cells_AUC_obese <- rbind(cells_AUC_obese, cells_AUC@assays@data@listData[["AUC"]])
}

###################################################################################################
###################################################################################################

###########################################
############### Diabetic obese ############ 
###########################################


# Define a list of gene sets with corresponding names
# Create an empty data frame to store results
cells_AUC_D_obese <- data.frame()

# Iterate through the gene sets
for (i in names(gene_lists)) {
  # Extract the current gene set and name
  currentGeneSet <- gene_lists[[i]]
  currentSetName <- as.character(i)  # You can customize the name as needed
  
  # Create a GeneSet object with the current name
  geneSetObj <- GeneSet(currentGeneSet, setName = currentSetName)
  
  # Calculating AUC score for the current GeneSet in the D_obese expression matrix
  cells_AUC <- AUCell_run(Vascular_D_obese_matrix, geneSetObj)
  
  # Build rankings
  cells_rankings <- AUCell_buildRankings(Vascular_D_obese_matrix, plotStats = TRUE)
  
  # Calculate AUC
  cells_AUC <- AUCell_calcAUC(geneSetObj, cells_rankings)
  
  # Append the AUC values to the result data frame
  cells_AUC_D_obese <- rbind(cells_AUC_D_obese, cells_AUC@assays@data@listData[["AUC"]])
}


###################################################################################################
###################################################################################################



## Final dataframe to plot the radar chart

AUC_list <- list(cells_AUC_lean, cells_AUC_obese, cells_AUC_D_obese)
names(AUC_list) <- c("Lean", "Obese", "D_obese")



mean_lean <- as.numeric(apply(cells_AUC_lean, 1, mean))
mean_obese <- as.numeric(apply(cells_AUC_obese, 1, mean))
mean_D_obese <- as.numeric(apply(cells_AUC_D_obese, 1, mean))

final_radar_df <- cbind(rownames(cells_AUC_lean), mean_lean, mean_obese, mean_D_obese) %>% as.data.frame()
colnames(final_radar_df)[1] <- "GO"
# rownames(final_radar_df) <- NULL
t_df <- t(final_radar_df)
t_df <- as.data.frame(t_df)
colnames(t_df) <- t_df[1,]
t_df <- t_df[-1,]

numeric_df <- apply(t_df, 2, as.numeric) %>% as.data.frame()
rownames(numeric_df) <- c("lean", "obese", "D_obese")

radar_subECs_numeric_df_UP <- numeric_df

Condition_colors <- c("#9FCECA","#7C4700","#191D54")

pdf("R_Projects/scAtlas/Analysis/Vascular/Analysis/DEA/LR/GO_terms/DEGs_Up/RadarPlots/subECs_common_UP_GOs_RadarPlot.pdf", 
    width = 20, height = 16)
fmsb::radarchart(radar_subECs_numeric_df_UP, maxmin = F, 
                 pcol = Condition_colors,
                 plty = 1, plwd = 1, pty = 16, 
                 pfcol = Condition_colors, cglwd = 1,
                 pdensity = 60, cglty = 2, title = "subECs", 
                 seg = 4, centerzero = T, palcex = 4, vlcex = 1)
dev.off()


list_radar_numeric_df_UP <- list("radar_subECs_numeric_df_UP" = radar_subECs_numeric_df_UP,
                                 "radar_pericytes_numeric_df_UP" = radar_pericytes_numeric_df_UP,
                                 "radar_venous1_numeric_df_UP" = radar_venous1_numeric_df_UP,
                                 "radar_venous2_numeric_df_UP" = radar_venous2_numeric_df_UP,
                                 "radar_VSMCs_numeric_df_UP" = radar_VSMCs_numeric_df_UP)

pdf("R_Projects/scAtlas/Analysis/Vascular/Analysis/DEA/LR/GO_terms/DEGs_Up/RadarChart_subECs_commonGO_letters.pdf", 
    width = 20, height = 16)
fmsb::radarchart(radar_subECs_numeric_df_letters, maxmin = F, 
                 pcol = Condition_colors, 
                 plty = 1, plwd = 1, pty = 16, 
                 pfcol = Condition_colors, 
                 pdensity = 60, cglty = 2, title = "subECs ECs", 
                 seg = 4, centerzero = T, palcex = 4, vlcex = 1)
dev.off()


