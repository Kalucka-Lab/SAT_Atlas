## Enrichment Lean vs Diabetic Obese
```{r}
Enrichment_Vascular_Clusters_DEG_DO_vs_L <- lapply(Vascular_Clusters_DEG_DO_vs_L, function(x) {
 enrichGO(
  gene= x$gene,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "ALL",
  pAdjustMethod = "BH",
  #universe = universe,
  qvalueCutoff = 0.05,
  pvalueCutoff = 0.05,
  readable=TRUE
)
})
```


## Enrichment Obese vs Diabetic Obese
```{r}
Enrichment_Vascular_Clusters_DEG_DO_vs_HO <- lapply(Vascular_Clusters_DEG_DO_vs_HO, function(x) {
 enrichGO(
  gene= x$gene,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "ALL",
  pAdjustMethod = "BH",
  #universe = universe,
  qvalueCutoff = 0.05,
  pvalueCutoff = 0.05,
  readable=TRUE
)
})
```


```{r, eval=FALSE}
df_Enrichment_Vascular_Clusters_DEG_DO_vs_HO <- lapply(Enrichment_Vascular_Clusters_DEG_DO_vs_HO, function(x){
  x <- as.data.frame(x)
})

df_Enrichment_Vascular_Clusters_DEG_DO_vs_L <- lapply(Enrichment_Vascular_Clusters_DEG_DO_vs_L, function(x){
  x <- as.data.frame(x)
})

df_Enrichment_Vascular_Clusters_DEG_HO_vs_L <- lapply(Enrichment_Vascular_Clusters_DEG_HO_vs_L, function(x){
  x <- as.data.frame(x)
})
```

# Dot plots
## Diabetic Obese vs Obese
```{r, fig.height=12, fig.width=12}
pdf("GO_terms/DO_vs_HO/DO_vs_HO_Dotplots.pdf", width = 12, height = 12)
lapply(Enrichment_Vascular_Clusters_DEG_DO_vs_HO, function(df) {
  if (nrow(df) > 0) {
    dotplot(df, showCategory = 15, split = "ONTOLOGY") + facet_grid(~ONTOLOGY) + theme(axis.text.y = element_text(size = 7))
     } else {
       NULL
     }
})
dev.off()
```


## Diabetic Obese vs Lean
```{r, fig.height=12, fig.width=12}
pdf("GO_terms/DO_vs_L/DO_vs_L_Dotplots.pdf", width = 12, height = 12)
lapply(Enrichment_Vascular_Clusters_DEG_DO_vs_L, function(df) {
  if (nrow(df) > 0) {
    dotplot(df, showCategory = 15, split = "ONTOLOGY") + facet_grid(~ONTOLOGY) + theme(axis.text.y = element_text(size = 7))
     } else {
       NULL
     }
})
dev.off()
```




## Lean vs Obese
```{r, fig.height=12, fig.width=12}
pdf("GO_terms/HO_vs_L/HO_vs_L_Dotplots.pdf", width = 12, height = 12)
lapply(Enrichment_Vascular_Clusters_DEG_HO_vs_L, function(df) {
  if (nrow(df) > 0) {
    dotplot(df, showCategory = 15, split = "ONTOLOGY") + facet_grid(~ONTOLOGY) + theme(axis.text.y = element_text(size = 7))
     } else {
       NULL
     }
})
dev.off()
```


