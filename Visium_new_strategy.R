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
#se_myo_core <- subset(se_myo, subset = annotation_myo == "core")

# For CTL
se_ctl$annotation_ctl <- rep("", ncol(se_ctl))
se_ctl <- FeatureViewer(se_ctl)


#saving
saveRDS(se_myo, "data/se_myo_annotated.rds")
saveRDS(se_myo_core, "data/se_myo_core.rds")
saveRDS(se_ctl, "data/se_ctl_annotated.rds")


# Load annotated Seurat objects from disk
se_myo <- readRDS("data/se_myo_annotated.rds")
se_ctl <- readRDS("data/se_ctl_annotated.rds")


# Inspect annotations (ensure lasso tool saved to 'selected_region')
table(se_myo$annotation_myo)
table(se_ctl$annotation_ctl)

# Subset to retain only spots labeled as "core"
se_myo_filtered <- subset(se_myo, subset = annotation_myo == "core")
se_ctl_filtered <- subset(se_ctl, subset = annotation_ctl == "core")


# Merge the filtered control and myocardial samples
se_clean2 <- merge(
  x = se_ctl_filtered,
  y = se_myo_filtered,
  add.cell.ids = c("ctl", "myo")
)
se_clean3 <- NormalizeData(se_clean2)
#se_clean3 <- FindVariableFeatures(se_clean3) #this was being used before

#making FindVariableFeatures less stringent and using some other strategy inside
se_clean3 <- FindVariableFeatures(
  se_clean3,
  selection.method = "vst",  # default method, reliable
  nfeatures = 4000           # default is 2000, increase as needed
)

se_clean3 <- FindVariableFeatures(
  se_clean3,
  selection.method = "dispersion",  # useful for sparse/high-dispersion markers
  mean.cutoff = c(0.05, 8),
  dispersion.cutoff = c(0.5, Inf)
)





se_clean3 <- ScaleData(se_clean3, split.by = "sample_id")
se_clean4 <- JoinLayers(se_clean3)


# Now continue as usual
se_clean <- RunPCA(se_clean4)
se_clean <- RunUMAP(se_clean, dims = 1:50, min.dist = 0.1, spread = 0.1)
se_clean <- FindNeighbors(se_clean, dims = 1:30)
se_clean <- FindClusters(se_clean, resolution = 2)
DimPlot(se_clean, split.by = "sample_id")



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
  return.thresh = 0.25
)

# View top genes per cluster
top_markers_per_cluster <- myo_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))

write.csv(top_markers_per_cluster, file = "top_myo_markers_per_cluster.csv", row.names = FALSE)


# Use Seurat's subset() instead
Idents(se_clean) <- "seurat_clusters"

#Identify marker genes per cluster
markers <- FindAllMarkers(
  se_clean,
  assay = "Spatial",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  return.thresh = 0.25
)

# View top genes per cluster
top_markers_per_cluster <- markers %>%
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


####Using Ab data


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



# Important: ab.capture = TRUE
se_all_ab <- ReadVisiumData(infoTable) |> LoadImages()


# Subset to myo/ctl
se_myo_ab <- SubsetSTData(se_all_ab, expression = sample_id == "myo")
se_ctl_ab <- SubsetSTData(se_all_ab, expression = sample_id == "ctl")



# Transfer annotations from previously saved lasso object
se_myo_ab$annotation_myo <- readRDS("data/se_myo_annotated.rds")$annotation_myo
se_ctl_ab$annotation_ctl <- readRDS("data/se_ctl_annotated.rds")$annotation_ctl



# Subset "core"
se_myo_filtered <- subset(se_myo_ab, subset = annotation_myo == "core")
se_ctl_filtered <- subset(se_ctl_ab, subset = annotation_ctl == "core")

# Merge
se_clean2 <- merge(
  x = se_ctl_filtered,
  y = se_myo_filtered,
  add.cell.ids = c("ctl", "myo")
)


# Set the correct assay
DefaultAssay(se_clean2) <- "AbCapture"

# Normalize using CLR
se_clean2 <- NormalizeData(se_clean2, normalization.method = "CLR", margin = 2)

# Find variable protein features
se_clean2 <- FindVariableFeatures(se_clean2)

# Scale data separately by sample
se_clean2 <- ScaleData(se_clean2, split.by = "sample_id")

# Run dimensionality reduction and clustering
se_clean2 <- RunPCA(se_clean2, assay = "AbCapture", features = VariableFeatures(se_clean2))
se_clean2 <- RunUMAP(se_clean2, dims = 1:20)
se_clean2 <- FindNeighbors(se_clean2, dims = 1:20)
se_clean2 <- FindClusters(se_clean2, resolution = 0.5)

# Visualize
DimPlot(se_clean2, split.by = "sample_id", label = TRUE)


rownames(se_clean2[["AbCapture"]])  # or use head() to preview




# Subset MYO only
se_myo_ab <- subset(se_clean2, subset = sample_id == "myo")
DefaultAssay(se_myo_ab) <- "AbCapture"

# Set cluster identities
Idents(se_myo_ab) <- "seurat_clusters"

# Find cluster-specific Ab markers (AbCapture only)
myo_ab_markers <- FindAllMarkers(
  se_myo_ab,
  assay = "AbCapture",
  only.pos = TRUE,
  logfc.threshold = 0.1,
  min.pct = 0.01
)

# Save or view top markers
write.csv(myo_ab_markers, "myo_ab_markers_per_cluster.csv", row.names = FALSE)


#checking cluster 3 since bigger in myo
se_myo_ab <- subset(se_clean2, subset = sample_id == "myo")
Idents(se_myo_ab) <- "seurat_clusters"

cluster3_myo_ab_markers <- FindMarkers(
  se_myo_ab,
  ident.1 = "3",
  assay = "AbCapture",
  only.pos = TRUE,
  logfc.threshold = 0.1,
  min.pct = 0.05
)


#compare cluster 3 (myo vs ctl) directly

se_c3_myo <- subset(se_clean2, subset = seurat_clusters == 3 & sample_id == "myo")
se_c3_ctl <- subset(se_clean2, subset = seurat_clusters == 3 & sample_id == "ctl")

se_c3_combined <- merge(se_c3_ctl, se_c3_myo, add.cell.ids = c("ctl", "myo"))
Idents(se_c3_combined) <- "sample_id"

de_c3_ab <- FindMarkers(
  se_c3_combined,
  ident.1 = "myo",
  ident.2 = "ctl",
  assay = "AbCapture",
  logfc.threshold = 0.25,
  min.pct = 0.05
)
write.csv(de_c3_ab, "cluster3_myo_vs_ctl_ab.csv")


#compare cluster 1 since did not find cytokine difference in cluster 3:
# Subset cluster 1 (myo vs ctl)
se_c1_myo <- subset(se_clean2, subset = seurat_clusters == 1 & sample_id == "myo")
se_c1_ctl <- subset(se_clean2, subset = seurat_clusters == 1 & sample_id == "ctl")

# Merge and set identities
se_c1_combined <- merge(se_c1_ctl, se_c1_myo, add.cell.ids = c("ctl", "myo"))
Idents(se_c1_combined) <- "sample_id"

# Find differential antibodies for cluster 1
de_c1_ab <- FindMarkers(
  se_c1_combined,
  ident.1 = "myo",
  ident.2 = "ctl",
  assay = "AbCapture",
  logfc.threshold = 0.25,
  min.pct = 0.05
)
write.csv(de_c1_ab, "cluster1_myo_vs_ctl_ab.csv")

#checking cluster 2 for cytokine differences now:
# Subset cluster 2 (myo vs ctl)
se_c2_myo <- subset(se_clean2, subset = seurat_clusters == 2 & sample_id == "myo")
se_c2_ctl <- subset(se_clean2, subset = seurat_clusters == 2 & sample_id == "ctl")

# Merge and set identities
se_c2_combined <- merge(se_c2_ctl, se_c2_myo, add.cell.ids = c("ctl", "myo"))
Idents(se_c2_combined) <- "sample_id"

# Find differential antibodies for cluster 2
de_c2_ab <- FindMarkers(
  se_c2_combined,
  ident.1 = "myo",
  ident.2 = "ctl",
  assay = "AbCapture",
  logfc.threshold = 0.25,
  min.pct = 0.05
)
write.csv(de_c2_ab, "cluster2_myo_vs_ctl_ab.csv")



#Identify your top marker per cluster
#From previous FindAllMarkers() output

top_markers <- myo_ab_markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

#plot these top_markers locations in the myo sample on raw histology image

# Step 1: Start from spatially valid object
se_all_ab <- ReadVisiumData(infoTable) |> LoadImages()

# Step 2: Subset MYO with spatial metadata intact
se_myo_ab <- SubsetSTData(se_all_ab, expression = sample_id == "myo")

# Step 3: Set AbCapture assay and normalize
DefaultAssay(se_myo_ab) <- "AbCapture"
se_myo_ab <- NormalizeData(se_myo_ab, normalization.method = "CLR", margin = 2)

# Step 4: Visualize top AbCapture markers (use actual top marker list)
top_marker_list <- c("CD68.1", "ITGAM.1", "PECAM1.1", "PTPRC.2", 
                     "ACTA2.1", "mouse-IgG2bk", "CEACAM8.1", 
                     "PCNA.1", "KRT5.1")

for (marker in top_marker_list) {
  p <- MapFeatures(
    object = se_myo_ab,
    features = marker,
    image_use = "raw",
    pt_size = 2,
    override_plot_dims = TRUE,
    title = paste("Spatial expression of", marker)
  )
  
  ggsave(
    filename = paste0("spatial_expression_", marker, ".png"),
    plot = p,
    width = 6, height = 6, dpi = 300
  )
}

