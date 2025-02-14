
library(tidyverse)
library(Seurat)
library(dittoSeq)

genelist <- read.table("R_Projects/scAtlas/Analysis/Vascular/Analysis/Druggable_genes/Druggable_genes_selected") %>% 
  pull(V1) %>% 
  unique() %>% 
  as.character()


Atlas_IntegratedObj <- readRDS("R_Projects/scAtlas/Analysis/BigObj/Objects/Atlas_IntegratedObj_with_all_labels_transferred.RDS")
Vascular_IntegratedObj <- readRDS("R_Projects/scAtlas/Analysis/Vascular/Objects/Vascular_IntegratedObj_label_transferred_from_SubECs_updated_names.RDS")

## Transferring labels
Atlas_IntegratedObj$transferred_lab_from_Vascular <- as.character(Atlas_IntegratedObj$final_annotation)
subset_cells <- rownames(Vascular_IntegratedObj@meta.data)
Atlas_IntegratedObj$transferred_lab_from_Vascular[rownames(Atlas_IntegratedObj@meta.data) %in% subset_cells] <- as.character(Vascular_IntegratedObj$Vascular_labels)


## Calculating AUC
Idents(Atlas_IntegratedObj) <- Atlas_IntegratedObj$transferred_lab_from_Vascular
vascular_labels_AUC <- FindAllMarkers(Atlas_IntegratedObj, 
                                      features = genelist, 
                                      test.use = "roc", 
                                      return.thresh = F)

mygenes <- readxl::read_xlsx("Downloads/genes_all_labels_AUC.xlsx", sheet = 1) %>% 
  select(myAUC,cluster, gene)

mygenes$cluster <- factor(mygenes$cluster)

cluster_gene_list <- split.data.frame(mygenes, mygenes$cluster)


## subsetting cells
arterial <- subset(Vascular_IntegratedObj, Vascular_labels == "Arterial Endothelial cells") %>% colnames()
venous1_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "Venous Endothelial cells 1") %>% colnames()
venous2_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "Venous Endothelial cells 2") %>% colnames()
venous3_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "Venous Endothelial cells 3") %>% colnames()
cap1_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "Capillary Endothelial cells 1") %>% colnames()
cap2_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "Capillary Endothelial cells 2") %>% colnames()
SubECs_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "SubEndothelial cells") %>% colnames()
lymphatic_cells <- subset(Vascular_IntegratedObj, Vascular_labels == "Lymphatic Endothelial cells") %>% colnames()
pericytes <- subset(Vascular_IntegratedObj, Vascular_labels == "Pericytes") %>% colnames()
vsmcs <- subset(Vascular_IntegratedObj, Vascular_labels == "VSMCs") %>% colnames()


############################
######## Arterial ##########
############################
pdf("RidgePlot_PCSK5.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "PCSK5", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = arterial) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_VEGFC.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "VEGFC", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = arterial) + NoLegend() + ylab("") + xlab("")
dev.off()

############################
######## Capillary 1 #######
############################
pdf("RidgePlot_FABP4_cap1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "FABP4", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap1_cells) + NoLegend() + ylab("") + xlab("")
dev.off()


############################
######## Capillary 2 #######
############################
pdf("RidgePlot_FABP4_cap2.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "FABP4", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_HLA-B.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "HLA-B", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_CD74.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "CD74", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_AQP1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "AQP1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_HSPB1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "HSPB1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_BCAM.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "BCAM", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = cap2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()


############################
######## Pericytes #########
############################
pdf("RidgePlot_POSTN.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "POSTN", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = pericytes) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_GUCY1A2.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "GUCY1A2", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = pericytes) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_NOTCH3.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "NOTCH3", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = pericytes) + NoLegend() + ylab("") + xlab("")
dev.off()



############################
######## Venous 2 #######
############################
pdf("RidgePlot_ADAMTS1_ven2.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "ADAMTS1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous2_cells) + NoLegend() + ylab("") + xlab("")
dev.off()


############################
######## Venous 3 ##########
############################
pdf("RidgePlot_BCAM_noDiabetic.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "BCAM", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_B2M_noDiabetic.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "B2M", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_HSP90AB1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "HSP90AB1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_GAPDH.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "GAPDH", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_ADAMTS1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "ADAMTS1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_CDKN1A.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "CDKN1A", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_IFITM1_noDiabetic.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "IFITM1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_DNAJA1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "DNAJA1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()

pdf("RidgePlot_HSPA1A.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "HSPA1A", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = venous3_cells) + NoLegend() + ylab("") + xlab("")
dev.off()



###############
#### VSMCS ####
###############
pdf("RidgePlot_PRKG1.pdf", width = 9, height = 3)
dittoRidgePlot(Vascular_IntegratedObj, var = "PRKG1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = vsmcs) + NoLegend() + ylab("") + xlab("")
dev.off()


dittoRidgePlot(Vascular_IntegratedObj, var = "PRKG1", 
               group.by = "Vascular_labels", 
               split.by = "Condition", 
               cells.use = vsmcs)


pdf("selected_52_genes_RidgePlots_VSMCs.pdf", width = 8, height = 8)
for (gene in mygenes){
  myplot <- dittoSeq::dittoRidgePlot(Vascular_IntegratedObj, var = gene, 
                                     group.by = "Vascular_labels", 
                                     split.by = "Condition", 
                                     cells.use = vsmcs)
  print(myplot)
}
dev.off()










adipo_like_genes <- cluster_gene_list$`Adipocyte-like Endothelial cells`$gene
lymphoid_like_genes <- cluster_gene_list$`Lymphoid-like Endothelial cells`$gene
mesenchymal_like_genes <- cluster_gene_list$`Mesenchymal-like Endothelial cells`$gene
myeloid1_like_genes <- cluster_gene_list$`Myeloid-like Endothelial cells 1`$gene
myeloid2_like_genes <- cluster_gene_list$`Myeloid-like Endothelial cells 2`$gene

adipo_like_cells <- subset(SubECsECs_Integrated, SubECsECs_labels == "Adipocyte-like Endothelial cells") %>% colnames()
mesenchymal_like_cells <- subset(SubECsECs_Integrated, SubECsECs_labels == "Mesenchymal-like Endothelial cells") %>% colnames()
lymphoid_like_cells <- subset(SubECsECs_Integrated, SubECsECs_labels == "Lymphoid-like Endothelial cells") %>% colnames()
myeloid1_like_cells <- subset(SubECsECs_Integrated, SubECsECs_labels == "Myeloid-like Endothelial cells 1") %>% colnames()
myeloid2_like_cells <- subset(SubECsECs_Integrated, SubECsECs_labels == "Myeloid-like Endothelial cells 2") %>% colnames()

##################################
######## Adipocyte-like ##########
##################################
pdf("RidgePlot_Adipocyte_like_genes.pdf", width = 9, height = 3)
for (gene in adipo_like_genes){
  myplot <- dittoSeq::dittoRidgePlot(SubECsECs_Integrated, var = gene, 
                                     group.by = "SubECsECs_labels", 
                                     split.by = "Condition", 
                                     cells.use = adipo_like_cells)
  print(myplot)
}
dev.off()


##################################
######## Lymphoid-like ##########
##################################
pdf("RidgePlot_lymphoid_like_genes.pdf", width = 9, height = 3)
for (gene in lymphoid_like_genes){
  myplot <- dittoSeq::dittoRidgePlot(SubECsECs_Integrated, var = gene, 
                                     group.by = "SubECsECs_labels", 
                                     split.by = "Condition", 
                                     cells.use = lymphoid_like_cells)
  print(myplot)
}
dev.off()





##################################
######## Mesenchymal-like ##########
##################################
pdf("RidgePlot_mesenchymal_like_genes.pdf", width = 9, height = 3)
for (gene in mesenchymal_like_genes){
  myplot <- dittoSeq::dittoRidgePlot(SubECsECs_Integrated, var = gene, 
                                     group.by = "SubECsECs_labels", 
                                     split.by = "Condition", 
                                     cells.use = mesenchymal_like_cells)
  print(myplot)
}
dev.off()





##################################
######## Myeloid1-like ##########
##################################
pdf("RidgePlot_Myeloid1_like_genes.pdf", width = 9, height = 3)
for (gene in myeloid1_like_genes){
  myplot <- dittoSeq::dittoRidgePlot(SubECsECs_Integrated, var = gene, 
                                     group.by = "SubECsECs_labels", 
                                     split.by = "Condition", 
                                     cells.use = myeloid1_like_cells)
  print(myplot)
}
dev.off()





##################################
######## Myeloid2-like ##########
##################################
pdf("RidgePlot_myeloid2_like_genes.pdf", width = 9, height = 3)
for (gene in myeloid2_like_genes){
  myplot <- dittoSeq::dittoRidgePlot(SubECsECs_Integrated, var = gene, 
                                     group.by = "SubECsECs_labels", 
                                     split.by = "Condition", 
                                     cells.use = myeloid2_like_cells)
  print(myplot)
}
dev.off()




