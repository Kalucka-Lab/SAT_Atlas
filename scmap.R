## River plot human vs pig

library(Seurat)

Human_vascular <- readRDS("R_Projects/scAtlas/Analysis/Vascular/Objects/Vascular_IntegratedObj_label_transferred.RDS")
Vascular_Pig <- readRDS("R_Projects/scAtlas/Pig/Vascular_Pig_annotated_updated.RDS")

query_object <- Vascular_Pig
ref_object <- Human_vascular


library(nichenetr)
library(dplyr)
ref_object <- as.SingleCellExperiment(ref_object, assay = "RNA")
query_object <- as.SingleCellExperiment(query_object, assay = "RNA")
#rownames(query_object) <- toupper(rownames(query_object))
# rownames(ref_object) <- toupper(rownames(ref_object))


shared_genes<-Reduce(intersect, lapply(list(rownames(query_object),rownames(ref_object)),
                                       FUN = function(x){toupper(rownames(x))}))


library(SingleCellExperiment)
library(scmap)
rowData(query_object)$feature_symbol <- rownames(query_object)
# remove features with duplicated names
query_object <- query_object[!duplicated(rownames(query_object)), ]
rowData(ref_object)$feature_symbol <- rownames(ref_object)
# remove features with duplicated names
ref_object <- ref_object[!duplicated(rownames(ref_object)), ]



library_list_combined<-list("query"=query_object,
                            "ref_1"=ref_object)


library(scmap)
library_list_combined[[1]]<-selectFeatures(library_list_combined[[1]],suppress_plot = F,n_features = 5000)
library_list_combined[[1]] <- indexCluster(library_list_combined[[1]],cluster_col = "Vascular_labels")
library_list_combined[[2]]<-selectFeatures(library_list_combined[[2]],suppress_plot = F,n_features = 5000)
library_list_combined[[2]] <- indexCluster(library_list_combined[[2]],cluster_col = "Vascular_labels")


#Select out the cluster index of each reference dataset
scmap_cluster_ref_dataset_cluster_index<-lapply(library_list_combined[c("query")],FUN=function(x){metadata(x)$scmap_cluster_index})
#Run scmapCluster, with threshold = 0.5
#change the threshold accordigly
scmapCluster_results_human_pig <- scmapCluster(
  projection = library_list_combined[["ref_1"]],
  index_list = scmap_cluster_ref_dataset_cluster_index,
  threshold = 0.1
)


scmap_colors_human_pig <- c("navy",
                            "#cd919e", 
                            "#cd3333", 
                            "darkgoldenrod", 
                            "#ffa54f", 
                            "#8b636c", 
                            "#cd853f",
                            "#6959cd", 
                            "forestgreen", 
                            "#473c8b")

plot(
  getSankey(
    colData(ref_object)$Vascular_labels, 
    scmapCluster_results_human_pig$scmap_cluster_labs,
    colors = scmap_colors_human_pig,
    plot_height = 400
  )
)
webshot2::webshot(url = "R_Projects/scAtlas/scMAP/Sankey_Human_pig_colored.html",
                  file = "/home/lucamannino/R_Projects/scAtlas/scMAP/Sankey_human_pig.pdf")

###############################################################################################
###############################################################################################


###################
## Mouse and human
###################
Human_vascular <- readRDS("R_Projects/scAtlas/Analysis/Vascular/Objects/Vascular_IntegratedObj_label_transferred.RDS")
Mouse_Vascular_Integrated <- readRDS("R_Projects/scAtlas/Mouse/Emont_only/Emont_Vascular_Integrated_annotated.RDS")

query_object <- Mouse_Vascular_Integrated
ref_object <- Human_vascular


library(nichenetr)
library(dplyr)
ref_object <- as.SingleCellExperiment(ref_object, assay = "RNA")
query_object <- as.SingleCellExperiment(query_object, assay = "RNA")
rownames(query_object) <- toupper(rownames(query_object))
# rownames(ref_object) <- toupper(rownames(ref_object))


shared_genes<-Reduce(intersect, lapply(list(rownames(query_object),rownames(ref_object)),
                                       FUN = function(x){toupper(rownames(x))}))


library(SingleCellExperiment)
library(scmap)
rowData(query_object)$feature_symbol <- rownames(query_object)
# remove features with duplicated names
query_object <- query_object[!duplicated(rownames(query_object)), ]
rowData(ref_object)$feature_symbol <- rownames(ref_object)
# remove features with duplicated names
ref_object <- ref_object[!duplicated(rownames(ref_object)), ]



library_list_combined<-list("query"=query_object,
                            "ref_1"=ref_object)


library(scmap)
library_list_combined[[1]]<-selectFeatures(library_list_combined[[1]],suppress_plot = F,n_features = 5000)
library_list_combined[[1]] <- indexCluster(library_list_combined[[1]],cluster_col = "vascular_labels")
library_list_combined[[2]]<-selectFeatures(library_list_combined[[2]],suppress_plot = F,n_features = 5000)
library_list_combined[[2]] <- indexCluster(library_list_combined[[2]],cluster_col = "Vascular_labels")


#Select out the cluster index of each reference dataset
scmap_cluster_ref_dataset_cluster_index<-lapply(library_list_combined[c("query")],FUN=function(x){metadata(x)$scmap_cluster_index})
#Run scmapCluster, with threshold = 0.5
#change the threshold accordigly
scmapCluster_results_human_mouse <- scmapCluster(
  projection = library_list_combined[["ref_1"]],
  index_list = scmap_cluster_ref_dataset_cluster_index,
  threshold = 0.1
)

scmap_colors_human_mouse <- c("#8b6508",
                              "navy",
                              "#eea9b8",
                              "#cd3333",
                              "#ffa54f",
                              "#cd853f",
                              "#8b636c",
                              "#0000cd",
                              "#00008b",
                              "forestgreen")

plot(
  getSankey(
    colData(ref_object)$Vascular_labels, 
    scmapCluster_results_human_mouse$scmap_cluster_labs,
    plot_height = 400,
    colors = scmap_colors_human_mouse
  )
)

webshot2::webshot(url = "R_Projects/scAtlas/scMAP/Sankey_Human_Emont_mouse_colored.html",
                  file = "/home/lucamannino/R_Projects/scAtlas/scMAP/Sankey_human_mouse.pdf")

save.image("R_Projects/scAtlas/scMAP_for_paper.Rdata")
