#####################################################################
########### New strategy 2.0, Lasso tool and all ####################
#####################################################################


library(Seurat)
library(semla)
library(patchwork)
library(dplyr)
library(ggplot2)


infoTable <- tibble::tibble(
  samples   = c("data/raw/220019_GEX_PEX/filtered_feature_bc_matrix.h5",
                "data/raw/230066_GEX_PEX/filtered_feature_bc_matrix.h5"),
  imgs      = c("data/raw/220019_GEX_PEX/spatial/tissue_lowres_image.png",
                "data/raw/230066_GEX_PEX/spatial/tissue_lowres_image.png"),
  spotfiles = c("data/raw/220019_GEX_PEX/spatial/tissue_positions.csv",
                "data/raw/230066_GEX_PEX/spatial/tissue_positions.csv"),
  json      = c("data/raw/220019_GEX_PEX/spatial/scalefactors_json.json",
                "data/raw/230066_GEX_PEX/spatial/scalefactors_json.json"),
  sample_id = c("myo", "ctl")
)



install.packages(c('shinyBS', 'beakr', 'colourpicker')) #needed for lasso tool to work
library("shinyBS", "beakr")
library("colourpicker")


# Load the full object
se_all <- ReadVisiumData(infoTable) |> LoadImages()

# Subset for each sample *without* splitting
# Subset the 'myo' sample
se_myo <- SubsetSTData(se_all, expression = sample_id == "myo")

# Subset the 'ctl' sample
se_ctl <- SubsetSTData(se_all, expression = sample_id == "ctl")


# For MYO
se_myo$annotation_myo <- rep("", ncol(se_myo))
se_myo <- FeatureViewer(se_myo)
se_myo_core <- subset(se_myo, subset = annotation_myo == "core")

# For CTL
se_ctl$annotation_ctl <- rep("", ncol(se_ctl))
se_ctl <- FeatureViewer(se_ctl)


#saving
saveRDS(se_myo, "data/se_myo_annotated.rds")
saveRDS(se_myo_core, "data/se_myo_core.rds")
saveRDS(se_ctl, "data/se_ctl_annotated.rds")

# Inspect annotations (ensure lasso tool saved to 'selected_region')
table(se_myo$annotation_myo)
table(se_ctl$annotation_ctl)

# Subset to retain only spots labeled as "core"
se_myo_filtered <- subset(se_myo, subset = annotation_myo == "core")
se_ctl_filtered <- subset(se_ctl, subset = annotation_ctl == "core")


# Merge the filtered control and myocardial samples
se_clean <- merge(
  x = se_ctl_filtered,
  y = se_myo_filtered,
  add.cell.ids = c("ctl", "myo")
)


# Step 7: Standard Seurat processing
# Normalize and scale MYO sample
# Normalize and scale each sample independently
se_myo_filtered <- NormalizeData(se_myo_filtered)
se_myo_filtered <- FindVariableFeatures(se_myo_filtered)
se_myo_filtered <- ScaleData(se_myo_filtered)

se_ctl_filtered <- NormalizeData(se_ctl_filtered)
se_ctl_filtered <- FindVariableFeatures(se_ctl_filtered)
se_ctl_filtered <- ScaleData(se_ctl_filtered)

# Merge preprocessed samples
se_clean <- merge(
  x = se_ctl_filtered,
  y = se_myo_filtered,
  add.cell.ids = c("ctl", "myo")
)

# Set assay
DefaultAssay(se_clean) <- "Spatial"

# Safest – re-run ScaleData after merge
# This ensures `scale.data` has proper cell names post-merge
se_clean <- ScaleData(se_clean)

# Now continue as usual
se_clean <- RunPCA(se_clean)
se_clean <- RunUMAP(se_clean, dims = 1:30)
se_clean <- FindNeighbors(se_clean, dims = 1:30)
se_clean <- FindClusters(se_clean, resolution = 0.8)


#Subset myo from merged object
# Use Seurat's subset() instead
se_myo <- subset(se_clean, subset = sample_id == "myo")

Idents(se_myo) <- "seurat_clusters"

#Identify marker genes per cluster
myo_markers <- FindAllMarkers(
  se_myo,
  assay = "Spatial",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  return.thresh = 0.05
)

# View top genes per cluster
top_markers_per_cluster <- myo_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))

write.csv(top_markers_per_cluster, file = "top_myo_markers_per_cluster.csv", row.names = FALSE)


#Evidence supporting immune infiltration in Cluster 3,4:
#From your top genes in cluster 3,4, several are strongly associated with myeloid immune cells, particularly:
  
#S100A9: a canonical marker of neutrophils and monocytes; involved in inflammation.

#CD36: expressed on macrophages and monocytes, involved in lipid uptake and innate immunity.

#CTSB (Cathepsin B): enriched in macrophages, especially in inflammatory or phagocytic states.

#CST3 (Cystatin C): secreted by many cells, but upregulated in infiltrating immune cells during inflammation.

#PSAP (Prosaposin): highly expressed in macrophages, involved in lysosomal function.
#Note: Absence of classic lymphoid markers (e.g., CD3D for T cells, MS4A1 for B cells).

#Comparing against ctl sample



# Step: Subset Cluster 3 from both MYO and CTL safely
se_myo_cluster3 <- subset(se_clean, subset = sample_id == "myo" & seurat_clusters == 3)
se_ctl_cluster3 <- subset(se_clean, subset = sample_id == "ctl" & seurat_clusters == 3)

# Step: Merge both for DE analysis
se_cluster3_combined <- merge(
  se_myo_cluster3,
  y = se_ctl_cluster3,
  add.cell.ids = c("myo", "ctl")
)

# Step: Set identity class to sample_id
Idents(se_cluster3_combined) <- "sample_id"

# Step: Define top immune-related cluster 3 genes
cluster3_top_genes <- c("S100A9", "CST3", "CD36", "PSAP", "CTSB")

# Step: Differential expression (MYO vs CTL in cluster 3)
se_cluster3_combined <- JoinLayers(se_cluster3_combined)


cluster3_deg <- FindMarkers(
  se_cluster3_combined,
  ident.1 = "myo",
  ident.2 = "ctl",
  features = cluster3_top_genes,
  assay = "Spatial",
  logfc.threshold = 0,
  min.pct = 0.05
)

# Step: Save and view results
write.csv(cluster3_deg, "cluster3_myo_vs_ctl_DE.csv")
head(cluster3_deg)

##           p_val  avg_log2FC pct.1 pct.2 p_val_adj
# S100A9 0.03258523 10.2184432 0.105 0.000         1
# CD36   0.47204635 -0.7135316 0.047 0.075         1
# CTSB   0.98850345  0.3098380 0.074 0.075         1

#same can be done for cluster 4



#Go enrichment of cluster markers
#################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages if not already installed
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

library(clusterProfiler)
library(org.Hs.eg.db)

# Filter marker genes from your data
markers_c3 <- myo_markers %>% filter(cluster == 3, p_val_adj < 0.05)
markers_c4 <- myo_markers %>% filter(cluster == 4, p_val_adj < 0.05)

# Run GO analysis
ego_c3 <- enrichGO(gene = markers_c3$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", readable = TRUE)
ego_c4 <- enrichGO(gene = markers_c4$gene, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", readable = TRUE)

# View top terms
barplot(ego_c3, showCategory = 10)
barplot(ego_c4, showCategory = 10)










