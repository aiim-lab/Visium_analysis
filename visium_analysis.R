install.packages("Matrix")
install.packages("Seurat")
install.packages("patchwork")
install.packages("dplyr")
install.packages("ggplot2")


library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

install.packages("BiocManager")  # If not already installed
BiocManager::install("hdf5r")

#Load data
seurat_obj <- Load10X_Spatial(
  data.dir = "C:/Users/manpr/OneDrive/Documents/Mcgill/manpreet_data/Visium/Visium/data/raw/230066_GEX_PEX",
  assay = "Spatial",
  slice = "230066_GEX_PEX"
)
####################################################
#checking if loading was fine
seurat_obj  # Should print Seurat object with dimensions
Assays(seurat_obj)  # Should include "Spatial"
DefaultAssay(seurat_obj)


# Check image slot (Visium image)
seurat_obj@images

# Plot tissue outline and spot positions
SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial")
##########################################################
#Calculate mitochondrial gene percentage (optional but recommended)
# Identify mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize key QC metrics
VlnPlot(seurat_obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
#RESULT: most spots are not dominated by mitochondrial transcripts

#scatterplot
# Scatterplot 1: nCount_Spatial vs percent.mt
FeatureScatter(seurat_obj, feature1 = "nCount_Spatial", feature2 = "percent.mt")

# Scatterplot 2: nCount_Spatial vs nFeature_Spatial
FeatureScatter(seurat_obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

############################################################################
#Apply QC Filters
#Let’s now filter out:
  
#Spots with <200 genes (likely empty/low-quality).

#Spots with >5,000 genes (likely doublets/artifacts).

#Spots with >10% mitochondrial gene expression.


seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_Spatial > 200 & 
                       nFeature_Spatial < 5000 & 
                       percent.mt < 10)

# View how many spots remain
seurat_obj
#RESULTS
#Spots before filtering: 1,417

#Spots after filtering: 1,190

#Removed 227 low-quality spots (≈16%).

################################################################################
# Normalize the data using LogNormalize
seurat_obj <- NormalizeData(seurat_obj, assay = "Spatial", normalization.method = "LogNormalize")

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

################################################################################
#Scale the data (Gene Expression layer)

seurat_obj <- ScaleData(seurat_obj, assay = "Spatial", features = VariableFeatures(seurat_obj))

################################################################################
# Run PCA on the variable genes
VariableFeatures(seurat_obj, layer = "counts.Gene Expression")[1:5]  # just show first 5

seurat_obj <- RunPCA(
  seurat_obj,
  assay = "Spatial",
  features = VariableFeatures(seurat_obj, layer = "counts.Gene Expression")
)

################################################################################
#Visualize PCA results
# Elbow plot to identify significant PCs
ElbowPlot(seurat_obj, ndims = 30)

# Optional PCA plot (2D)
DimPlot(seurat_obj, reduction = "pca", group.by = NULL)

################################################################################
#Run UMAP and Clustering
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, reduction = "pca")

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

################################################################################
#Visualize Clusters
# UMAP plot colored by clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP: Spatial Clusters")

# Spatial plot to see cluster locations on the tissue
SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3) + ggtitle("Spatial Clusters")

################################################################################
# Identify marker genes across spatial clusters
# Join layers for proper differential expression testing
seurat_obj <- JoinLayers(seurat_obj, assay = "Spatial")


markers <- FindAllMarkers(
  seurat_obj,
  assay = "Spatial",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# View top markers per cluster
markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)


# Show top 3 markers per cluster by log fold change
markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

################################################################################
#score immune cell types using the 31-plex protein panel 

# Set the path to the features.tsv.gz file for your sample
features_path <- "C:/Users/manpr/OneDrive/Documents/Mcgill/manpreet_data/Visium/Visium/data/raw/230066_GEX_PEX/filtered_feature_bc_matrix/features.tsv"

# Read features file (assumes tab-separated with no header)
features_df <- read.table(gzfile(features_path), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(features_df) <- c("ID", "Name", "Modality")  # third column = modality

# Check modality distribution
table(features_df$Modality)

# Pull full count matrix
full_counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")

# Match features
rna_idx <- which(features_df$Modality == "Gene Expression")
adt_idx <- which(features_df$Modality == "Antibody Capture")

rna_counts <- full_counts[rna_idx, ]
adt_counts <- full_counts[adt_idx, ]

# Add separate assays
seurat_obj[["RNA"]] <- CreateAssayObject(counts = rna_counts)
seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)
DefaultAssay(seurat_obj) <- "RNA"

# Check success
Assays(seurat_obj)



# Check dimensions and some rownames (protein markers)
ADT_counts <- GetAssayData(seurat_obj, assay = "ADT", layer = "counts")
dim(ADT_counts)                       # Should be ~35 x 1190
rownames(ADT_counts)[1:10]           # First 10 markers

# Check a few values for selected markers
ADT_counts[c("CD3E", "CD4", "CD8A"), 1:5]  # Expression for first 5 spots

# See all ADT markers
rownames(ADT_counts)

# Try matching names with suffix ".1"
ADT_counts[c("CD3E.1", "CD4.1", "CD8A.1"), 1:5]


################################################################################
#Normalizing
seurat_obj <- NormalizeData(seurat_obj, assay = "ADT", normalization.method = "CLR")

################################################################################
DefaultAssay(seurat_obj) <- "ADT"


Matrix::colSums(GetAssayData(seurat_obj, slot = "counts"))[1:10]
Matrix::rowSums(GetAssayData(seurat_obj, slot = "counts"))  # For overall marker expression


library(Matrix)
mtx_path <- "C:/Users/manpr/OneDrive/Documents/Mcgill/manpreet_data/Visium/Visium/data/raw/230066_GEX_PEX/filtered_feature_bc_matrix/matrix.mtx.gz"
mtx <- readMM(mtx_path)
dim(mtx)  # Should match total features × spots (e.g., 18120 x 1417)

# Check last 35 rows for Antibody Capture data
adt_block <- mtx[(nrow(mtx)-34):nrow(mtx), ]
summary(adt_block)  # Should show non-zero entries if data exists






################################################################################
######## Chris code ############################################################


# library ----------------------------------------------------------------
install.packages("viridis")
install.packages("remotes")

install.packages("remotes")
remotes::install_github("ludvigla/semla")

# Load it
library(semla)

install.packages("Seurat")


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(semla)
library(viridis)



data_root_directory <- list.dirs("data/raw", full.names = TRUE, recursive = FALSE)

samples <- Sys.glob(paths = file.path(data_root_directory, 
                                      "filtered_feature_bc_matrix.h5"))

imgs <- Sys.glob(paths = file.path(data_root_directory, 
                                   "spatial", "tissue_lowres_image.png"))

spotfiles <- Sys.glob(paths = file.path(data_root_directory, 
                                        "spatial", "tissue_positions.csv"))

json <- Sys.glob(paths = file.path(data_root_directory, 
                                   "spatial", "scalefactors_json.json"))

infoTable <- tibble(samples, imgs, spotfiles, json, # Add required columns
                    sample_id = c("myo", "ctl")) # Add additional column

se <- semla::ReadVisiumData(infoTable)
saveRDS(se, file = "data/results/se_before_QC_MS.rds")

se <- readRDS(file = "data/results/se_before_QC_MS.rds")




se <- semla::LoadImages(se)
semla::ImagePlot(se)





# Plot with semla
MapFeatures(se, features = "nFeature_Spatial", 
            image_use = "raw",
            override_plot_dims = TRUE,
            pt_size = 1.5,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral")) + 
  ThemeLegendRight()

cols <- c("myo" = "#E41A1C", "ctl" = "#377EB8")  # Customize if needed


p <- MapFeaturesSummary(se, 
                        features = "nFeature_Spatial", 
                        subplot_type = "violin",
                        colors = cols)
# Visualize QC metrics as a violin plot
plot1 <- VlnPlot(se, features = "nCount_Spatial", pt.size = 0.1, log = TRUE) + NoLegend()
plot2 <- SpatialFeaturePlot(se, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


cols <- viridis::rocket(11, direction = -1)

# Plot with Seurat
MapFeatures(se, slot = "counts",
            features = "VIM", 
            colors = cols)
p




MapFeaturesSummary(se, features = "nFeature_Spatial", subplot_type = "histogram")
# Filter by number of unique genes
se <- SubsetSTData(se, expression = nFeature_Spatial > 40)
se
MapFeaturesSummary(se[,se@meta.data$sample_id=="myo"], features = "nFeature_Spatial", subplot_type = "histogram")

tt <- SubsetSTData(se, expression = sample_id == "ctl", spots = NULL, features = NULL, idents = NULL)




se1 <- se |> 
  NormalizeData() |>
  ScaleData() |> 
  FindVariableFeatures() |> 
  RunPCA()





MapFeatures(se1, 
            features = "PC_4", 
            center_zero = TRUE, 
            section_number = 1, 
            pt_size = 2,
            colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") |> rev())




se1 <- FindNeighbors(se1, reduction = "pca", dims = 1:30)
se1 <- FindClusters(se1, verbose = FALSE)
se1 <- RunUMAP(se1, reduction = "pca", dims = 1:30)





p1 <- DimPlot(se1, reduction = "umap", label = TRUE, split.by = "sample_id")
p2 <- SpatialDimPlot(se1, label = TRUE, label.size = 3)
p1 + p2








se <- LoadImages(se, verbose = FALSE)

MapLabels(se1, 
          column_name = "seurat_clusters", 
          image_use = "raw", 
          override_plot_dims = TRUE) +
  plot_layout(guides = "collect") &
  guides(fill = guide_legend(override.aes = list(size = 3), 
                             ncol = 2)) &
  theme(legend.position = "right")




se2 <- se1

# Normalize ADT data,
DefaultAssay(se2) <- "AbCapture"
se2 <- NormalizeData(se2, normalization.method = "CLR", margin = 2)
DefaultAssay(cbmc) <- "RNA"

# Note that the following command is an alternative but returns the same result
s2 <- NormalizeData(se2, normalization.method = "CLR", margin = 2, assay = "AbCapture")

# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2
```


```{r}
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(se2, "spatial_CD14", cols = c("lightgrey", "darkgreen")) + ggtitle("CD68 RNA")
p2 <- FeaturePlot(se2, "abcapture_CD14.1") + ggtitle("CD68 protein")
p1 | p2
```


```{r}
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(se2, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 100) %>%
  ungroup() -> top100
DoHeatmap(se2, features = top100$gene) + NoLegend()
```

########################################################
####### Own code: inspired by Chris code ##############

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(semla)

# Define actual folder names and sample IDs
sample_dirs <- c("data/raw/220019_GEX_PEX", "data/raw/230066_GEX_PEX")
sample_ids  <- c("myo", "ctl")

# Load each sample and assign sample ID
sample_list <- lapply(seq_along(sample_dirs), function(i) {
  obj <- Load10X_Spatial(data.dir = sample_dirs[i])
  obj$sample_id <- sample_ids[i]
  obj
})

# Merge into a single Seurat object
se <- merge(sample_list[[1]], y = sample_list[[2]], add.cell.ids = sample_ids, project = "visium_myo_ctl")

table(se$sample_id)

#### QC ####
## Violin plots (gene/UMI distribution across samples): P1

VlnPlot(se, features = c("nFeature_Spatial", "nCount_Spatial"), 
        group.by = "sample_id", pt.size = 0.1, log = TRUE) +
  plot_layout(ncol = 2)


##Spatial feature plot- visualize tissue

library(patchwork)

# Subset by sample
se_ctl <- subset(se, subset = sample_id == "ctl")
se_myo <- subset(se, subset = sample_id == "myo")

# Plot separately
p_ctl <- SpatialFeaturePlot(se_ctl, features = "nCount_Spatial") + ggtitle("ctl")
p_myo <- SpatialFeaturePlot(se_myo, features = "nCount_Spatial") + ggtitle("myo")

# Combine side-by-side
p_ctl + p_myo   #p2


# Remove spots with fewer than 100 genes OR fewer than 500 UMIs
# First, split the original unfiltered object
se_ctl_raw <- subset(se, subset = sample_id == "ctl")
se_myo_raw <- subset(se, subset = sample_id == "myo")

# Apply different thresholds
#se_ctl_filt <- subset(se_ctl_raw, subset = nFeature_Spatial > 100 & nCount_Spatial > 500)
#se_myo_filt <- subset(se_myo_raw, subset = nFeature_Spatial > 40 & nCount_Spatial > 300)

# Apply different thresholds for ctl and stricter for myo based on visual inspection
se_ctl_filt <- subset(se_ctl_raw, subset = nFeature_Spatial > 100 & nCount_Spatial > 500)
se_myo_filt <- subset(se_myo_raw, subset = nFeature_Spatial > 40 & nCount_Spatial > 400)  # increased from 300


# Combine back into one object
se_filtered <- merge(se_ctl_filt, y = se_myo_filt, add.cell.ids = c("ctl", "myo"), project = "visium_filtered")


# Check final filtered spot counts
table(se_filtered$sample_id)

# Visual inspection
p_ctl <- SpatialFeaturePlot(subset(se_filtered, sample_id == "ctl"), features = "nFeature_Spatial") + ggtitle("ctl")
p_myo <- SpatialFeaturePlot(subset(se_filtered, sample_id == "myo"), features = "nFeature_Spatial") + ggtitle("myo")
p_ctl + p_myo




# Using semla here for fun for normalization, pca etc

infoTable <- tibble(
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

# Read data with semla (returns Seurat object)
se <- semla::ReadVisiumData(infoTable)

se <- semla::LoadImages(se)


se <- NormalizeData(se)
se <- FindVariableFeatures(se)
se <- ScaleData(se)
se <- RunPCA(se)

# After PCA and before plotting
se <- RunUMAP(se, reduction = "pca", dims = 1:30)


# Clustering and dimensionality reduction
se <- FindNeighbors(se, reduction = "pca", dims = 1:30)
se <- FindClusters(se, resolution = 0.8, verbose = FALSE)
se <- RunUMAP(se, reduction = "pca", dims = 1:30)

# Now plotting
p1 <- DimPlot(se, reduction = "umap", label = TRUE, split.by = "sample_id") #saved as p4 locally


#Map clusters on histology (after clustering + umap)

table(se$seurat_clusters)
# Create a numeric version of seurat_clusters
se$seurat_clusters_numeric <- as.numeric(as.character(se$seurat_clusters))

# Now try MapFeatures on the numeric version
MapFeatures(se, features = "seurat_clusters_numeric")

## Normalize Antibody captyre data

# Set to AbCapture (protein) assay
DefaultAssay(se) <- "AbCapture"

# Normalize the protein data
se <- NormalizeData(se, normalization.method = "CLR", margin = 2)

# Visualize feature names if unsure what proteins are present
rownames(se[["AbCapture"]])[1:10]



# Replace "CD163" with your gene of interest, and "CD163.1" with matching protein
gene <- "CD163"
protein <- "CD163.1"

DefaultAssay(se) <- "Spatial"
p_rna <- FeaturePlot(se, gene, cols = c("lightgrey", "darkgreen")) + ggtitle(paste(gene, "RNA"))

DefaultAssay(se) <- "AbCapture"
p_protein <- FeaturePlot(se, protein, cols = c("lightgrey", "darkred")) + ggtitle(paste(gene, "Protein"))

p_rna | p_protein



# Continue from semla pipeline after clustering and UMAP embedding
# Goal: Visualize clusters on tissue images

# Convert cluster labels to numeric (for MapFeatures)
se$seurat_clusters_numeric <- as.numeric(as.character(se$seurat_clusters))

# Visualize clusters over histology
p9 <- MapFeatures(
  se,
  features = "seurat_clusters_numeric",
  image_use = "raw",
  override_plot_dims = TRUE,
  pt_size = 2,
  colors = RColorBrewer::brewer.pal(n = 9, name = "Set1")
) +
  ggtitle("Clusters mapped to tissue") +
  theme(legend.position = "right")

p9


#Find marker genes for each cluster

# Set Spatial assay as default
DefaultAssay(se) <- "Spatial"

# Identify marker genes for each cluster (positive markers only)
cluster_markers <- FindAllMarkers(se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top 5 markers per cluster
library(dplyr)
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

head(top_markers, 15)


#Visualize top markers on UMAP

# Pick top genes to visualize
top_genes <- unique(top_markers$gene)[1:6]  # or manually choose representative ones

# Plot on UMAP
FeaturePlot(se, features = top_genes, reduction = "umap", ncol = 3)

# Plot on histology (one by one or loop through)
MapFeatures(se, features = "PDK4", image_use = "raw", pt_size = 1.5)  # example gene

#Map marker expression to histology

library(viridis)

# Visualize spatial expression of representative markers from each cluster
MapFeatures(se,
            features = c("PDK4"),
            image_use = "raw",
            override_plot_dims = TRUE,
            pt_size = 1.5,
            colors = viridis::magma(100)) +
  plot_layout(ncol = 1)


#########################################################################
##### seeing top expressed genes #############

# Set Spatial assay
DefaultAssay(se) <- "Spatial"

# Split by sample
se_myo <- subset(se, subset = sample_id == "myo")
se_ctl <- subset(se, subset = sample_id == "ctl")

# Compute average expression per gene
avg_myo <- rowMeans(as.matrix(GetAssayData(se_myo, slot = "data")))
avg_ctl <- rowMeans(as.matrix(GetAssayData(se_ctl, slot = "data")))

# Get top 20 genes for each
top_myo <- sort(avg_myo, decreasing = TRUE)[1:500]
top_ctl <- sort(avg_ctl, decreasing = TRUE)[1:500]

# View gene names and expression
top_myo_genes <- names(top_myo)
top_ctl_genes <- names(top_ctl)

# Print unique and overlapping genes
cat("Top 20 MYO genes:\n"); print(top_myo_genes)
cat("\nTop 20 CTL genes:\n"); print(top_ctl_genes)

cat("\nCommon top genes:\n"); print(intersect(top_myo_genes, top_ctl_genes))
cat("\nUnique to MYO:\n"); print(setdiff(top_myo_genes, top_ctl_genes))
cat("\nUnique to CTL:\n"); print(setdiff(top_ctl_genes, top_myo_genes))


#### checking of the significant genes are statistically differentailly expressed

Idents(se) <- "sample_id"
de_myo_vs_ctl <- FindMarkers(se, ident.1 = "myo", ident.2 = "ctl", min.pct = 0.25, logfc.threshold = 0.25)
head(de_myo_vs_ctl)


######################################################################

##### here the goal is to perform proper immune profiling of the 2 sample ###
#1) identify immune-relevant protein markers from your AbCapture panel,
#(2) visualize them in MYO vs CTL side-by-side on tissue using MapFeatures().

##############################################################################

# Get all available protein markers in the AbCapture assay
all_proteins <- rownames(se[["AbCapture"]])
cat("All available proteins:\n")
print(all_proteins)


# Load required packages
library(semla)
library(ggplot2)
library(patchwork)
library(viridis)

# Set sample IDs
samples <- c("myo", "ctl")

# Define marker groups
rna_tcell_genes <- c("CD3E")#, "CD4", "CD8A")
protein_tcell_genes <- c("CD3E.1")#, "CD4.1", "CD8A.1")

# Function to generate plots
plot_marker_group <- function(se, markers, assay_name, sample_id, out_prefix) {
  DefaultAssay(se) <- assay_name
  
  plots <- lapply(markers, function(gene) {
    p <- MapFeatures(subset(se, sample_id == sample_id),
                     features = gene,
                     image_use = "raw",
                     pt_size = 1.5,
                     override_plot_dims = TRUE,
                     colors = viridis::magma(100)) +
      ggtitle(paste(sample_id, "-", gene))
    return(p)
  })
  
  # Combine and save
  combined <- wrap_plots(plots, ncol = 1)
  ggsave(paste0("plots/", out_prefix, "_", sample_id, ".png"),
         combined, width = 5, height = 12, dpi = 300)
  
  return(combined)
}

# --- T Cell RNA plots ---
plot_marker_group(se, rna_tcell_genes, "Spatial", "myo", "tcell_rna")
plot_marker_group(se, rna_tcell_genes, "Spatial", "ctl", "tcell_rna")

# --- T Cell Protein plots ---
plot_marker_group(se, protein_tcell_genes, "AbCapture", "myo", "mouse-IgG2a")
plot_marker_group(se, protein_tcell_genes, "AbCapture", "ctl", "tcell_protein")


VlnPlot(se, features = c("CD3E.1", "CD4.1", "CD8A.1"), group.by = "seurat_clusters")

########################################################################
######### Violin plots to check difference in immune expression profiles


# Load required library
library(ggplot2)

# Ensure correct assay and identities
DefaultAssay(se) <- "AbCapture"
Idents(se) <- "sample_id"

# Define protein markers for each immune group
tcell_proteins <- c("CD3E.1", "CD4.1", "CD8A.1")
bcell_proteins <- c("CD19.1", "MS4A1.1")
macrophage_proteins <- c("CD68.1", "CD163.1")

# Combine all immune markers
all_immune_proteins <- c(tcell_proteins, bcell_proteins, macrophage_proteins)

# Plot violin plots
v_immune <- VlnPlot(
  object = se,
  features = all_immune_proteins,
  pt.size = 0.1,
  group.by = "sample_id",
  combine = TRUE,
  stack = TRUE,
  flip = TRUE,
  ncol = 1
) + ggtitle("Immune Protein Expression (AbCapture): MYO vs CTL")

# Save the plot
ggsave("violin_t_b_macrophages.png", plot = v_immune, width = 10, height = 12, dpi = 300)

# Optionally, print to view immediately
print(v_immune)


########################################################################
#### checking RNA expression levels for the same markers ###############

library(Seurat)
library(ggplot2)
library(patchwork)

# Set identities if not already done
Idents(se) <- "sample_id"

# Define RNA markers for each immune cell type
rna_tcell_genes <- c("CD3E", "CD4", "CD8A")
rna_bcell_genes <- c("CD19", "MS4A1")
rna_macrophage_genes <- c("CD68", "CD163")

# Function to display violin plots
view_violin_group <- function(se, gene_list, group_name) {
  VlnPlot(se, features = gene_list, group.by = "sample_id", pt.size = 0.1) +
    plot_layout(ncol = 1) +
    plot_annotation(title = paste("RNA Expression -", group_name))
}

# 🧪 Display RNA expression for each immune group
view_violin_group(se, rna_tcell_genes, "T Cells")
view_violin_group(se, rna_bcell_genes, "B Cells")
view_violin_group(se, rna_macrophage_genes, "Macrophages")



#################################################
##### checking for other markers now ############


# List of immune markers (RNA) to check
immune_genes_to_check <- c(
  # Dendritic cells
  "ITGAX", "HLA-DRA", 
  # NK cells
  "NCAM1", "KLRD1", 
  # Neutrophils
  "CEACAM8", "FCGR3B", 
  # Monocytes
  "CD14", "LYZ", "FCGR3A"
)

# Check which of these genes are present in the RNA (Spatial) assay
DefaultAssay(se) <- "Spatial"
available_genes <- immune_genes_to_check[immune_genes_to_check %in% rownames(se)]
missing_genes   <- setdiff(immune_genes_to_check, available_genes)

cat("✅ Available RNA genes:\n")
print(available_genes)

cat("\n❌ Missing from dataset:\n")
print(missing_genes)


#Now violin plots for protein levels 


# Define immune marker groups (protein)
protein_markers <- list(
  dendritic = c("ITGAX.1"),
  neutrophil = c("CEACAM8.1", "FCGR3B.1"),
  monocyte = c("CD14.1", "LYZ.1", "FCGR3A.1")
)

# Generate violin plots for each immune group
for (group_name in names(protein_markers)) {
  print(paste("📊 Violin plot for:", group_name))
  markers <- protein_markers[[group_name]]
  
  # Plot and show
  p <- VlnPlot(se, features = markers, group.by = "sample_id", pt.size = 0.1, assay = "AbCapture") +
    ggtitle(paste("Protein Expression -", group_name))
  
  print(p)
}


# Now violin plots for RNA expression levels


# Define the RNA markers to inspect
rna_markers_to_plot <- c("CD14", "ITGAX") #because visual differences were seen here

# Set the default assay to Spatial (RNA)
DefaultAssay(se) <- "Spatial"

# Generate and display violin plots for RNA markers
p_rna_diff <- VlnPlot(
  object = se,
  features = rna_markers_to_plot,
  group.by = "sample_id",
  pt.size = 0.1
) +
  ggtitle("RNA Expression: CD14 and ITGAX (MYO vs CTL)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Print the plot
print(p_rna_diff)



#checking of protein levels are statostcially different

# Set protein markers of interest
markers_to_test <- c("CD14.1", "ITGAX.1")

# Set assay to AbCapture (protein)
DefaultAssay(se) <- "AbCapture"

# Perform Wilcoxon test for each marker
for (gene in markers_to_test) {
  expr <- FetchData(se, vars = gene)
  groups <- se$sample_id
  
  # Perform Wilcoxon rank-sum test
  test_result <- wilcox.test(expr[groups == "myo", 1],
                             expr[groups == "ctl", 1])
  
  cat("\n", gene, "Wilcoxon test p-value:", signif(test_result$p.value, 4), "\n")
}



#######################################################
########## New stratergy ##############################
#######################################################

library(Seurat)
library(semla)
library(patchwork)

# Define sample info
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

# Load data and images
se <- semla::ReadVisiumData(infoTable)
se <- semla::LoadImages(se)

# Split samples and filter
se_ctl_filt <- subset(se, sample_id == "ctl" & nFeature_Spatial > 100 & nCount_Spatial > 500)
se_myo_filt <- subset(se, sample_id == "myo" & nFeature_Spatial > 40 & nCount_Spatial > 400)

# Merge back
se_filtered <- merge(se_ctl_filt, y = se_myo_filt, add.cell.ids = c("ctl", "myo"))

# Reattach image data for spatial plotting
se_filtered@images$ctl <- se@images$ctl
se_filtered@images$myo <- se@images$myo

# Optional: visualize filtering outcome
p_ctl <- SpatialFeaturePlot(subset(se_filtered, sample_id == "ctl"), features = "nFeature_Spatial") + ggtitle("ctl")
p_myo <- SpatialFeaturePlot(subset(se_filtered, sample_id == "myo"), features = "nFeature_Spatial") + ggtitle("myo")
p_ctl + p_myo

# Proceed with Seurat workflow
DefaultAssay(se_filtered) <- "Spatial"
se_filtered <- NormalizeData(se_filtered)
se_filtered <- FindVariableFeatures(se_filtered)
se_filtered <- ScaleData(se_filtered)
se_filtered <- RunPCA(se_filtered)
se_filtered <- RunUMAP(se_filtered, dims = 1:30)
se_filtered <- FindNeighbors(se_filtered, dims = 1:30)
se_filtered <- FindClusters(se_filtered, resolution = 0.8)

###############################################################
#Visualize MYO sample clusters
# Set identities to clusters
Idents(se_filtered) <- "seurat_clusters"

# Split UMAP by sample to inspect MYO clusters
DimPlot(se_filtered, reduction = "umap", label = TRUE, split.by = "sample_id") + ggtitle("Clusters by sample")

#display these clusters on histology image

library(semla)
library(Seurat)

# Step 1: Load data with images
se <- semla::ReadVisiumData(infoTable)
se <- semla::LoadImages(se)

# Step 2: Apply filtering directly to original object (do NOT split or merge)
se <- SubsetSTData(
  se,
  expression = case_when(
    se$sample_id == "ctl" ~ se$nFeature_Spatial > 100 & se$nCount_Spatial > 500,
    se$sample_id == "myo" ~ se$nFeature_Spatial > 40 & se$nCount_Spatial > 400,
    TRUE ~ FALSE  # drop anything else
  )
)



# Step 3: Proceed with processing
se <- NormalizeData(se)
se <- FindVariableFeatures(se)
se <- ScaleData(se)
se <- RunPCA(se)
se <- RunUMAP(se, dims = 1:30)
se <- FindNeighbors(se, dims = 1:30)
se <- FindClusters(se, resolution = 0.8)

# Step 4: Assign numeric labels
Idents(se) <- "seurat_clusters"
se$seurat_clusters_numeric <- as.numeric(as.character(se$seurat_clusters))

# Step 5: Spatial cluster mapping
MapFeatures(
  object = se,
  features = "seurat_clusters_numeric",
  image_use = "raw",
  pt_size = 2,
  override_plot_dims = TRUE,
  colors = viridis::inferno(length(unique(se$seurat_clusters)))
) + ggtitle("Spatial mapping of clusters") +
  theme(legend.position = "right")

##############################################################
########### try masking ######################################
# why because table(se$sample_id)

#ctl  myo 
#1417 1138 

#AND

#table(se_filtered$sample_id)

#ctl  myo 
#1057  155 
#although the filtering runs just fine but the number of spots left on the myo sample are very few


# Re-attach images from original object if they were removed
library(semla)

# Step 1: Read and load images (you've likely done this already)
se <- semla::ReadVisiumData(infoTable)
se <- semla::LoadImages(se)

# Step 2: Apply automated masking (masks background tissue)
se <- MaskImages(
  se,
  thresholding = TRUE,
  iso.blur = 2,          # controls smoothing
  channels.use = 1:3,    # typical for Visium: RGB channels
  compactness = 1,       # adjust superpixel density
  add.contrast = TRUE,   # enhance contrast before SLIC
  verbose = TRUE
)

# Step 3: Visualize masked image
# Just call ImagePlot() directly on the Seurat object
ImagePlot(se)
#####################################################
#Observed that masking is just relevant for visual enhancement, it does not remove the spots
#which are not overlapping with tissues sample, the numbers remain same!

## Now checking something else 
#I will remove the clusters and corresponding spots that are off tissue

# After masking the object:
se <- NormalizeData(se)
se <- FindVariableFeatures(se)
se <- ScaleData(se)
se <- RunPCA(se)
se <- RunUMAP(se, dims = 1:30)
se <- FindNeighbors(se, dims = 1:30)
se <- FindClusters(se, resolution = 0.8)
Idents(se) <- "seurat_clusters"


se$seurat_clusters_numeric <- as.numeric(as.character(se$seurat_clusters))
semla::MapFeatures(
  object = se,
  features = "seurat_clusters_numeric",
  image_use = "raw",  # valid option in your version
  pt_size = 2,
  override_plot_dims = TRUE,
  colors = viridis::inferno(length(unique(se$seurat_clusters)))
)

#Quantify Spot Counts per Cluster per Sample:
table(se$sample_id, se$seurat_clusters)
#       0    1    2    3    4
#ctl 1114  222   37   24   20
#myo   97  504  233  235   69

#Manually Remove Off-Tissue Clusters (i want to remove cluster 2 from myo , 2 and 4 from ctl):
# Ensure correct identity is set
Idents(se) <- "seurat_clusters"

# Create a metadata column for clarity if not already present
if (!"sample_id" %in% colnames(se@meta.data)) {
  stop("Missing 'sample_id' column in metadata.")
}

# Flag clusters to exclude (specific to each sample)
se$exclude <- with(se@meta.data,
                   (sample_id == "myo" & seurat_clusters == "2") |
                     (sample_id == "ctl" & seurat_clusters %in% c("2", "4")))

# Use semla's safer method to subset while keeping spatial image alignment
se_clean <- semla::SubsetSTData(se, expression = !se$exclude)


se_clean$seurat_clusters_numeric <- as.numeric(as.character(se_clean$seurat_clusters))

semla::MapFeatures(
  object = se_clean,
  features = "seurat_clusters_numeric",
  image_use = "raw",  # or "masked"
  pt_size = 2,
  override_plot_dims = TRUE,
  colors = viridis::inferno(length(unique(se_clean$seurat_clusters)))
)

table(se_clean$sample_id, se_clean$seurat_clusters)
#       0    1    3    4
#ctl 1114  222   24    0
#myo   97  504  235   69


#subset myo sample
se_myo <- semla::SubsetSTData(se_clean, expression = sample_id == "myo")

Idents(se_myo) <- "seurat_clusters"

myo_markers <- FindAllMarkers(
  se_myo,
  assay = "Spatial",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  return.thresh = 0.05
)

# Get top gene from clusters


# Top 30 genes per cluster
top_genes_per_cluster <- myo_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

top_genes_per_cluster

#isloate cluster 3 (looks immune related)
immune_genes_cluster3 <- myo_markers %>%
  filter(cluster == "3") %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30)


#manually checking for known cytokines in myo sample:
cytokine_genes <- c("TNF", "IL1B", "IL6", "CCL2", "CXCL8", "CXCL10", "IFNG", "IL10", "TGFB1")

MapFeatures(
  object = se_myo,
  features = cytokine_genes,
  image_use = "raw",
  pt_size = 2,
  override_plot_dims = TRUE
)

#To catch low-level but significant cytokine expression:
markers_cluster3 <- FindMarkers(se_myo, ident.1 = "3", ident.2 = NULL)
markers_cluster3 %>% filter(gene %in% cytokine_genes)

#Cluster 3 is immunologically distinct from other clusters in the myo sample 
#and likely represents an inflamed region with active immune responses, 
#supported by the overexpression of genes like S100A9, CD36, and CFLAR.

###checking across sample (myo vs ctl)

# Subset se_clean by sample
se_myo <- subset(se_clean, subset = sample_id == "myo")
se_ctl <- subset(se_clean, subset = sample_id == "ctl")

# Set identity class for both
Idents(se_myo) <- "seurat_clusters"
Idents(se_ctl) <- "seurat_clusters"

#Extract cluster 3 from both
# Subset only cluster 3 from both samples
c3_myo <- subset(se_myo, idents = "3")
c3_ctl <- subset(se_ctl, idents = "3")

# Compare markers between C3 in myo and ctl
# Option A: Just compare mean expression of top immune genes
top_immune_genes <- c("S100A9", "HSPB6", "CD36", "TMSB4X", "CFLAR", "S100A6")  # adjust based on your prior list

# Average expression in each C3 subset
avg_myo <- AverageExpression(c3_myo, features = top_immune_genes, assays = "Spatial")$Spatial
avg_ctl <- AverageExpression(c3_ctl, features = top_immune_genes, assays = "Spatial")$Spatial

# Combine for visualization
comparison_df <- data.frame(
  Gene = rownames(avg_myo),
  Myo = avg_myo[, 1],
  Ctl = avg_ctl[, 1]
)

# Visualize as barplot
library(ggplot2)
comparison_df_long <- tidyr::pivot_longer(comparison_df, cols = c("Myo", "Ctl"), names_to = "Sample", values_to = "Expression")

ggplot(comparison_df_long, aes(x = Gene, y = Expression, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Cluster 3: Myo vs Ctl", y = "Avg. Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Plot clearly demonstrates that Cluster 3 in the myo sample shows enhanced 
#immune-related gene expression compared to the same cluster in ctl



#Spatially visualize key DE genes (e.g., S100A9, TMSB4X, HSPB6)
se_myo <- semla::SubsetSTData(se_clean, expression = sample_id == "myo")

immune_genes <- c("S100A9", "TMSB4X", "HSPB6")

plot_list <- list()

# Generate each plot and store
for (gene in immune_genes) {
  p <- MapFeatures(
    object = se_myo,
    features = gene,
    image_use = "raw",
    pt_size = 2,
    override_plot_dims = TRUE
  ) + ggtitle(paste("Spatial Expression of", gene))
  
  plot_list[[gene]] <- p
}

# Combine using patchwork: side-by-side
combined_plot <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
  plot_layout(ncol = 3)

# Display
print(combined_plot)
