########################################################### Spatial Transcriptomics Analysis ########################################################################################
#Load libraries
library(Seurat)
library(SeuratObject)

Sample_1A <- readRDS("path to file /Sample_1Aannot.rds")
Sample_1D <- readRDS("path to file /Sample_1Dannot.rds")
Sample_2A <- readRDS("path to file /Sample_2Aannot.rds")
Sample_2D <- readRDS("path to file /Sample_2Dannot.rds")
Sample_3A <- readRDS("path to file /Sample_3Aannot.rds")
Sample_3D <- readRDS("path to file /Sample_3Dannot.rds")
Sample_4A <- readRDS("path to file /Sample_4Aannot.rds")
Sample_4D <- readRDS("path to file /Sample_4Dannot.rds")

############ Normalise all Seurat objects ############
# List of Seurat objects
seurat_objects <- list(Sample_1A, Sample_1D, Sample_2D, Sample_2A, Sample_3A, Sample_3D, Sample_4A, Sample_4D)  # Add all your Seurat objects here

# Loop through each Seurat object
for (seurat_obj in seurat_objects) {
  # Set the DefaultAssay to "Spatial"
  DefaultAssay(seurat_obj) <- "Spatial"
  
  # Normalize the data for the Spatial assay
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE, assay = "Spatial")
}


######## SCALOP analysis to annotate tumour regions ###########
#Load libraries
library(scalop)
library(readr)

# List of Seurat object names
seurat_names <- c("Sample_1A", "Sample_1D", "Sample_2A", "Sample_2D", "Sample_3A", "Sample_3D", "Sample_4A", "Sample_4D")

#load signatures from Puchalski et al. (2018)
Tumour_region_signatures <- read_csv("path to file /Puchalski_sig.csv")

LE_genes <- Tumour_region_signatures[["LE"]]
CTmvp_genes <- Tumour_region_signatures[["CTmvp"]]
CTpan_genes <- Tumour_region_signatures[["CTpan"]]
CT_genes <- Tumour_region_signatures[["CT"]]

# List of signatures
sigs <- list(CT_genes, CTmvp_genes, CTpan_genes, LE_genes)
names(sigs) <- c("CT_genes", "CTmvp_genes", "CTpan_genes", "LE_genes")

# Loop through each sample in seurat_names
for (sample_name in seurat_names) {
  # Load the Seurat object if not already in the environment
  seurat_obj <- get(sample_name)
  
  # Basic Scoring of Matrix by Gene sigs
  m <- as.matrix(GetAssayData(object = seurat_obj[["Spatial"]], slot = "data"))
  
  m_scores <- baseScores(m, sigs, conserved.genes = 0.7)
  
  # Score a Matrix by Gene sigs (Signatures)
  m_sigscores <- sigScores(m, sigs, groups = NULL, center.rows = TRUE,
                           center = TRUE,
                           expr.center = TRUE,
                           expr.bin.m = NULL,
                           expr.bins = NULL,
                           expr.sigs = NULL,
                           expr.nbin = 30,
                           expr.binsize = 100,
                           conserved.genes = 0.7,
                           replace = FALSE
  )
  
  vector_max_sig <- maxcol_strict(m_sigscores, min = NULL, diff = NULL, splitByCol = FALSE)
  
  highest_Sig_Column <- as.data.frame(colnames(m_sigscores)[max.col(m_sigscores, ties.method = "first")])
  row.names(highest_Sig_Column) <- row.names(m_sigscores)
  highest_Sig_Column$Spot <- rownames(highest_Sig_Column)
  
  # Match the order of Seurat barcodes with your data
  myBarcode <- rownames(seurat_obj@meta.data)
  highest_Sig_Column_2 <- highest_Sig_Column[match(myBarcode, highest_Sig_Column$Spot), ]
  
  # Add the calculated values to the Seurat object's metadata
  seurat_obj$Tumour_region_SpotIDs <- highest_Sig_Column_2[, 1] # Adjust the column index if needed
  
  # Optionally, save or assign the modified Seurat object back to the environment
  assign(sample_name, seurat_obj)
  
  # Progress update (optional)
  print(paste("Processed:", sample_name))
}


############# Add SCALOP cell type sig scores as metadata columns #####################################

#Load libraries
library(scalop)

# List of Seurat object names
seurat_names <- c("Sample_1A", "Sample_1D", "Sample_2A", "Sample_2D", "Sample_3A", "Sample_3D", "Sample_4A", "Sample_4D")

# Gene lists for Niche features of interest
NPC.features <- c("STMN2", "CD24", "RND3", "HMP19", "TUBB3", "MIAT", "DCX", "NSG1", "ELAVL4", "MLLT11", "DLX6-AS1", "SOX11", "NREP", "FNBP1L", "TAGLN3", "STMN4", "DLX5", "SOX4", "MAP1B", "RBFOX2", "IGFBPL1", "STMN1", "HN1", "TMEM161B-AS1", "DPYSL3", "SEPT3", "PKIA", "ATP1B1", "DYNC1I1", "CD200", "SNAP25", "PAK3", "NDRG4", "KIF5A", "UCHL1", "ENO2", "KIF5C", "DDAH2", "TUBB2A", "LBH", "LOC150568", "TCF4", "GNG3", "NFIB", "DPYSL5", "CRABP1", "DBN1", "NFIX", "CEP170", "BLCAP","DLL3", "DLL1", "SOX4", "TUBB3", "HES6", "TAGLN3", "NEU4", "MARCKSL1", "CD24", "STMN1", "TCF12", "BEX1", "OLIG1", "MAP2", "FXYD6", "PTPRS", "MLLT11", "NPPA", "BCAN", "MEST", "ASCL1", "BTG2", "DCX", "NXPH1", "HN1", "PFN2", "SCG3", "MYT1", "CHD7", "GPR56", "TUBA1A", "PCBP4", "ETV1", "SHD", "TNR", "AMOTL2", "DBN1", "HIP1", "ABAT", "ELAVL4", "LMF1", "GRIK2", "SERINC5", "TSPAN13", "ELMO1", "GLCCI1", "SEZ6L", "LRRN1", "SEZ6", "SOX11")
AC.features <- c("CST3", "S100B", "SLC1A3", "HEPN1", "HOPX", "MT3", "SPARCL1", "MLC1", "GFAP", "FABP7", "BCAN", "PON2", "METTL7B", "SPARC", "GATM", "RAMP1", "PMP2", "AQP4", "DBI", "EDNRB", "PTPRZ1", "CLU", "PMP22", "ATP1A2", "S100A16", "HEY1", "PCDHGC3", "TTYH1", "NDRG2", "PRCP", "ATP1B2", "AGT", "PLTP", "GPM6B", "F3", "RAB31", "PPAP2B", "ANXA5", "TSPAN7")
OPC.features <- c("BCAN", "PLP1", "GPR17", "FIBIN", "LHFPL3", "OLIG1", "PSAT1", "SCRG1", "OMG", "APOD", "SIRT2", "TNR", "THY1", "PHYHIPL", "SOX2-OT", "NKAIN4", "LPPR1", "PTPRZ1", "VCAN", "DBI", "PMP2", "CNP", "TNS3", "LIMA1", "CA10", "PCDHGC3", "CNTN1", "SCD5", "P2RX7", "CADM2", "TTYH1", "FGF12", "TMEM206", "NEU4", "FXYD6", "RNF13", "RTKN", "GPM6B", "LMF1", "ALCAM", "PGRMC1", "HRASLS", "BCAS1", "RAB31", "PLLP", "FABP5", "NLGN3", "SERINC5", "EPB41L2", "GPR37L1")
MES.features <- c("CHI3L1", "ANXA2", "ANXA1", "CD44", "VIM", "MT2A", "C1S", "NAMPT", "EFEMP1", "C1R", "SOD2", "IFITM3", "TIMP1", "SPP1", "A2M", "S100A11", "MT1X", "S100A10", "FN1", "LGALS1", "S100A16", "CLIC1", "MGST1", "RCAN1", "TAGLN2", "NPC2", "SERPING1", "C8orf4", "EMP1", "APOE", "CTSB", "C3", "LGALS3", "MT1E", "EMP3", "SERPINA3", "ACTN1", "PRDX6", "IGFBP7", "SERPINE1", "PLP2", "MGP", "CLIC4", "GFPT2", "GSN", "NNMT", "TUBA1C", "GJA1", "TNFRSF1A", "WWTR1", "HILPDA", "ADM", "DDIT3", "NDRG1", "HERPUD1", "DNAJB9", "TRIB3", "ENO2", "AKAP12", "SQSTM1", "MT1X", "ATF3", "NAMPT", "NRN1", "SLC2A1", "BNIP3", "LGALS3", "INSIG2", "IGFBP3", "PPP1R15A", "VIM", "PLOD2", "GBE1", "SLC2A3", "FTL", "WARS", "ERO1L", "XPOT", "HSPA5", "GDF15", "ANXA2", "EPAS1", "LDHA", "P4HA1", "SERTAD1", "PFKP", "PGK1", "EGLN3", "SLC6A6", "CA9", "BNIP3L", "RPL21", "TRAM1", "UFM1", "ASNS", "GOLT1B", "ANGPTL4", "SLC39A14", "CDKN1A", "HSPA9")
Endothelial <- c("CLDN5", "CCL14", "C7", "FLG", "ITM2A", "ISLR", "COL11A1", "CPZ", "LUM", 
                 "IGFBP7", "COL4A1", "VWF", "BMX", "CHRNA1", "DCN", "ANGPT2", "COL1A1", 
                 "CCN2", "IGKV1-27", "ABCG2", "SEMA3G", "COL3A1", "IGHV4-39", "CAVIN2", 
                 "EPAS1", "COL6A3", "IGHV3-30", "SPAAR", "FLT1", "SPARC", "FCN3", "CCL13", 
                 "APOLD1", "EDN1", "SRARP", "MEOX1", "EDN3", "CLEC3B", "IFI27", "FOXC2", 
                 "IGKV3-15", "CA2", "HSPG2", "ABCA9", "IGHV3-23", "LINC01391", "CLEC14A", 
                 "ADGRL4", "CFHR3", "SOX18")
Pericyte <- c("DCN", "COL3A1", "COL1A1", "MGP", "COL1A2", "APOD", "TAGLN", "IGFBP7", 
              "TIMP1", "PTGDS", "CALD1", "TPM2", "COL4A1", "FN1", "IFITM3", "OGN", 
              "CYP1B1", "BGN", "SFRP2", "ACTA2", "RGS5", "COL6A3", "SPARC", "COL6A2", 
              "FBLN1", "MYL9", "TPM1", "CCN2", "LUM", "SLC26A2", "IGFBP5", "COL6A1", 
              "COL4A2", "NDUFA4L2", "IGFBP4", "PCOLCE", "C1S", "FSTL1", "COL5A2", 
              "LEPR", "SPARCL1", "CTHRC1", "TIMP3", "SELENOM", "UACA", "NNMT", "C7", 
              "NID1", "PDGFRB", "HIGD1B")
Myeloid <- c("TGFBI", "S100A9", "RNASE1", "FTL", "CCL3", "TK1", "CH25H", "HLA-DRB5", 
             "FCN1", "S100A12", "DHRS9", "BIRC5", "SELENOP", "MT1G", "FN1", "APOC1", 
             "CXCL3", "CXCL2", "IBSP", "CCL18", "S100A8", "VCAN", "SPP1", "PPBP", 
             "H4C3", "IGKV1-39", "SERPINB10", "ADAMTS5", "LYZ", "C1QB", "UBE2C", "CCL4L2", 
             "CCNB1", "SFRP4", "SRGN", "ADM", "C1QA", "CCL4", "MKI67", "STMN1", "THBS1", 
             "C15orf48", "CDC20", "MT1H", "AREG", "HLA-DRB1", "CTSD", "C1QC", "EREG", "IL1B")
Lymphocyte <- c("CCL5", "ISG20", "GZMH", "IL7R", "GZMK", "GNLY", "CD40LG", "KLRB1", 
                "IL32", "LTB", "KLRC1", "FGFBP2", "SH2D1B", "FKBP11", "IFITM1", "CD79A", 
                "IGKC", "GZMA", "CD8A", "NKG7", "LAG3", "KRT81", "IGLC3", "CD8B", "TNFRSF4", 
                "GZMB", "CCR7", "MAL", "TRDC", "IGHG1", "SPON2", "FLT3LG", "CD2", "XCL2", 
                "XCL1", "IL2RA", "CD3E", "FOXP3", "TRDV2", "IGHA2", "KLF2", "TRBV27", "ANXA1", 
                "CD6", "PRF1", "AQP3", "IGLC2", "CTLA4", "KLRC2", "GPR171")

sigs <- list(NPC.features, AC.features, OPC.features, MES.features, Endothelial, Pericyte, Myeloid, Lymphocyte)
names(sigs) <- c("NPC.features", "AC.features", "OPC.features", "MES.features", "Endothelial", "Pericyte", "Myeloid", "Lymphocyte")

# Loop through each Seurat object in seurat_names
for (seurat_name in seurat_names) {
  cat("Processing", seurat_name, "\n")
  
  # Load the Seurat object
  seurat_obj <- get(seurat_name)
  
  # Extract spatial assay data
  m <- as.matrix(GetAssayData(object = seurat_obj[["Spatial"]], slot = "data"))
  
  # Compute signature scores
  m_sigscores <- sigScores(m, sigs, groups = NULL, center.rows = TRUE,
                           center = TRUE,
                           expr.center = TRUE,
                           expr.bin.m = NULL,
                           expr.bins = NULL,
                           expr.sigs = NULL,
                           expr.nbin = 30,
                           expr.binsize = 100,
                           conserved.genes = 0.7,
                           replace = FALSE)
  
  # Convert m_sigscores to a dataframe and rename columns for Seurat metadata
  m_sigscores_df <- as.data.frame(m_sigscores)
  colnames(m_sigscores_df) <- paste0(colnames(m_sigscores_df), "_SCALOP")
  
  # Ensure row names (barcodes) match Seurat metadata row names
  m_sigscores_df$Spot <- rownames(m_sigscores_df)
  seurat_metadata <- seurat_obj@meta.data
  seurat_metadata$Spot <- rownames(seurat_metadata)
  
  # Merge metadata and m_sigscores
  seurat_metadata <- merge(seurat_metadata, m_sigscores_df, by = "Spot", all.x = TRUE)
  
  # Restore row names after merging
  rownames(seurat_metadata) <- seurat_metadata$Spot
  seurat_metadata$Spot <- NULL  # Remove temporary Spot column
  
  # Assign back updated metadata
  seurat_obj@meta.data <- seurat_metadata
  
  # Assign updated Seurat object
  assign(seurat_name, seurat_obj)
}

# Check the metadata for Sample_1A
head(Sample_1A@meta.data)

############################ Harmony integration ##################################################
# Load necessary libraries
library(Seurat)
library(SeuratDisk)
library(harmony)

# List of spatial Seurat objects
seurat_list <- list(
  "Sample_1A" = Sample_1A, "Sample_1D" = Sample_1D,
  "Sample_2A" = Sample_2A, "Sample_2D" = Sample_2D,
  "Sample_3A" = Sample_3A, "Sample_3D" = Sample_3D,
  "Sample_4A" = Sample_4A, "Sample_4D" = Sample_4D
)

# Ensure SCT is the active assay and add "Sample" metadata to each object
seurat_list <- lapply(names(seurat_list), function(sample_name) {
  obj <- seurat_list[[sample_name]]
  DefaultAssay(obj) <- "SCT"
  obj$Sample <- sample_name  # Add sample name as metadata
  return(obj)
})

# Check if SCT normalization exists and re-run if needed
seurat_list <- lapply(seurat_list, function(x) {
  if (ncol(x@assays$SCT@scale.data) == 0) {
    x <- SCTransform(x, verbose = FALSE)
  }
  return(x)
})

# **Extract variable features before merging**
variable_features_list <- lapply(seurat_list, VariableFeatures)
common_variable_features <- Reduce(intersect, variable_features_list)

# Merge all Seurat objects into one, ensuring unique cell names
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1])

# Ensure SCT is the default assay after merging
DefaultAssay(merged_seurat) <- "SCT"

# Set "Sample" metadata correctly after merging
merged_seurat$Sample <- factor(merged_seurat$Sample)

# **Manually Set Variable Features from Individual Objects**
VariableFeatures(merged_seurat) <- common_variable_features

# Print number of selected variable features
print(paste("Number of variable features after merging:", length(VariableFeatures(merged_seurat))))

# Run PCA
merged_seurat <- RunPCA(merged_seurat, npcs = 30, verbose = FALSE)

# Run Harmony batch correction using the "Sample" metadata
merged_seurat <- RunHarmony(merged_seurat, group.by.vars = "Sample", dims.use = 1:30)

# Run UMAP
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)

# Find neighbors and clusters
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.51)

# Specify the file path and name for the RDS file
file_path <- "C:/Users/shardcl/OneDrive - University of South Australia/Peter Janes/Bioinformatic/New_analysis_2025/Rds_files/Harmony_Integrated_Spatial_Peter.rds"

# Save the Seurat object as an RDS file
saveRDS(merged_seurat, file = file_path)

# Spatial visualization
library(ggplot2)

# Default UMAP Plot (Clusters)
p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP - Clustered")

# UMAP Plot Colored by Sample (Batch Effect Check)
p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "Sample", label = FALSE) +
  ggtitle("UMAP - Colored by Sample")

# UMAP Plot Colored by Tumour region
p3 <- DimPlot(merged_seurat, reduction = "umap", group.by = "Tumour_region_SpotIDs", label = FALSE) +
  ggtitle("UMAP - Colored by Sample")

# UMAP Plot Colored by EPHA3 expression
p4 <- FeaturePlot(merged_seurat, features = c("EPHA3"))

# **Display Plots**
print(p1)
print(p2)
print(p3)
p1 + p2
p1 + p4

############### Individual sample spatial plots #######################################################

# Define the sample name you want to plot
selected_sample <- "Sample_1A"  # Change this to your sample of interest

# Subset the merged Seurat object to keep only cells from this sample
subset_seurat <- subset(merged_seurat, subset = Sample == selected_sample)

# Spatial DimPlot for just this sample
SpatialDimPlot(subset_seurat, group.by = "Tumour_region_SpotIDs", pt.size.factor = 103)
SpatialDimPlot(subset_seurat, group.by = "Pericyte_SCALOP_thresholded", pt.size.factor = 103)
SpatialFeaturePlot(subset_seurat, features = "EPHA3", min.cutoff = 0, max.cutoff = 1.5, pt.size.factor = 103)
SpatialFeaturePlot(subset_seurat, features = "MES.features_SCALOP", pt.size.factor = 103, min.cutoff = -0.4, max.cutoff = 1)


################ EPHA3+ spots per cluster per sample ##############################################

library(Seurat)
library(tidyr)
library(pheatmap)
library(tibble)

# Create an empty list
epha3_distribution_list <- list()

# Loop through each sample
for (sample_name in unique(merged_seurat$Sample)) {
  # Subset by sample
  obj_sub <- subset(merged_seurat, subset = Sample == sample_name)
  
  # Get metadata and EPHA3 expression
  meta <- obj_sub@meta.data
  meta$EPHA3_expr <- FetchData(obj_sub, vars = "EPHA3")[, 1]
  meta$EPHA3_positive <- ifelse(meta$EPHA3_expr > 0, 1, 0)
  
  # Total number of EPHA3+ spots in the sample
  total_EPHA3_pos <- sum(meta$EPHA3_positive)
  
  # Group by cluster to calculate each cluster's contribution
  df_summary <- meta %>%
    group_by(seurat_clusters) %>%
    summarise(
      EPHA3_pos_spots = sum(EPHA3_positive),
      region_contrib_to_total_EPHA3_pos = ifelse(
        total_EPHA3_pos > 0,
        100 * EPHA3_pos_spots / total_EPHA3_pos,
        NA
      ),
      .groups = "drop"
    )
  
  df_summary$Sample <- sample_name
  epha3_distribution_list[[sample_name]] <- df_summary
}

# Combine into a single data frame
EPHA3_cluster_distribution <- do.call(rbind, epha3_distribution_list)

# OPTIONAL: Check if each sample sums to ~100
check_sums <- EPHA3_cluster_distribution %>%
  group_by(Sample) %>%
  summarise(sum_percent = sum(region_contrib_to_total_EPHA3_pos, na.rm = TRUE))

print(check_sums)

# Prepare matrix for heatmap
heatmap_matrix <- EPHA3_cluster_distribution %>%
  select(Sample, seurat_clusters, region_contrib_to_total_EPHA3_pos) %>%
  pivot_wider(names_from = seurat_clusters, values_from = region_contrib_to_total_EPHA3_pos) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# Plot
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize_row = 10,
         fontsize_col = 10,
         color = colorRampPalette(c("white", "red"))(100),
         main = "% Contribution of Each Cluster to Total EPHA3+ Spots")


# Calculate average across tissues
average_percent_per_cluster <- colMeans(heatmap_matrix, na.rm = TRUE)

average_percent_per_cluster


############### Cell type signature quantification per cluster (all samples merged) ###########################

library(dplyr)

# Define the metadata columns of interest (signature scores)
signature_columns <- c("NPC.features_SCALOP", "AC.features_SCALOP", "OPC.features_SCALOP", 
                       "MES.features_SCALOP", "Endothelial_SCALOP", "Pericyte_SCALOP", 
                       "Myeloid_SCALOP", "Lymphocyte_SCALOP")

# Extract metadata and cluster info
signature_scores_df <- merged_seurat@meta.data %>%
  select(seurat_clusters, all_of(signature_columns))

# Compute average signature score per cluster
avg_signature_scores_per_cluster <- signature_scores_df %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(signature_columns), ~mean(.x, na.rm = TRUE)))

# View results
print(avg_signature_scores_per_cluster)



#################################################################### Single cell Analysis ########################################################################################

#Load necessary libraries
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

#load datasets
#Ebert
obj_Ebert <- readRDS("path to file/Ebert.rds")

#LeBlanc
obj_LeBlanc <- readRDS("path to file/LeBlanc.rds")
#remove recurrent samples
Idents(obj_LeBlanc) <- "sample"
obj_LeBlanc <- subset(obj_LeBlanc, idents = c("JK124_reg1_tis_1", "JK124_reg1_tis_2", "JK124_reg2_tis_1", "JK124_reg2_tis_2", "JK125_reg1_tis_1.1", "JK125_reg1_tis_1.2", "JK125_reg2_tis_1", "JK125_reg2_tis_2_r1", "JK126_reg1_tis_1.1", "JK126_reg1_tis_1.2", "JK126_reg2_tis_1", "JK134_reg1_tis_1", "JK134_reg2_tis_1", "JK136_reg1_tis_1", "JK136_reg2_tis_1", "JK136_reg2_tis_2_br", "JK142_reg1_tis_1", "JK142_reg2_tis_1", "JK142_reg2_tis_2.1_br", "JK142_reg2_tis_2.2_br", "JK152_reg1_tis_1", "JK152_reg2_tis_1", "JK153_reg1_tis_1", "JK153_reg2_tis_1", "JK156_reg1_tis_1", "JK156_reg2_tis_1", "JK156_reg2_tis_2_br", "JK163_reg1_tis_1", "JK163_reg2_tis_1"))

#Abdelfattah
obj_Abdelfattah <- readRDS("path to file/Abdelfattah.rds")
#remove recurrent samples
Idents(obj_Abdelfattah) <- "Type"
obj_Abdelfattah <- subset(obj_Abdelfattah, idents = c("GBM"))

# Define a mapping for each dataset's metadata column to standardized cell types
map_to_pooled_celltype <- function(original_celltype, dataset_name) {
  mapping <- list(
    "Ebert" = c(
      "Endothelial cells_2" = "Endothelial",
      "Macrophages_5" = "Myeloid",
      "Pericyte cells_4" = "Pericyte",
      "T cells_1" = "Lymphocyte",
      "Tumor cells_0" = "Malignant",
      "Tumor/prolif. cells_3" = "Myeloid"
    ),
    "LeBlanc" = c(
      "endothelial" = "Endothelial",
      "fibroblast" = "Pericyte",
      "immune" = "Myeloid",
      "malignant" = "Malignant",
      "neuron" = "Other",
      "oligodendrocyte" = "Oligodendrocyte"
    ),
    "Abdelfattah" = c(
      "BCells" = "Lymphocyte",
      "Endo" = "Endothelial",
      "Glioma" = "Malignant",
      "Myeloid" = "Myeloid",
      "Oligo" = "Oligodendrocyte",
      "Other" = "Other",
      "Pericytes" = "Pericyte",
      "TCells" = "Lymphocyte"
    )
    )
  
  # Handle NA values as "Other"
  if (is.na(original_celltype)) {
    return("Other")
  }
  
  # Perform the mapping
  standardized <- mapping[[dataset_name]][original_celltype]
  
  # Return mapped value or "Other" for unmapped categories
  return(ifelse(is.na(standardized), "Other", standardized))
}

# Add Pooled_celltype to each dataset
add_pooled_celltype <- function(seurat_obj, original_metadata, dataset_name) {
  seurat_obj@meta.data$Pooled_celltype <- sapply(
    seurat_obj@meta.data[[original_metadata]],
    map_to_pooled_celltype,
    dataset_name
  )
  return(seurat_obj)
}

# Apply to all datasets
obj_Ebert <- add_pooled_celltype(obj_Ebert, "combined_clusters_annot", "Ebert")
obj_LeBlanc <- add_pooled_celltype(obj_LeBlanc, "cell_type", "LeBlanc")
obj_Abdelfattah <- add_pooled_celltype(obj_Abdelfattah, "Assignment", "Abdelfattah")

# View results
table(obj_Ebert$Pooled_celltype)
table(obj_LeBlanc$Pooled_celltype)
table(obj_Abdelfattah$Pooled_celltype)

FeaturePlot(obj_Ebert, features = "EPHA3", pt.size = 1)
DimPlot(obj_Ebert, group.by = "Pooled_celltype")

########### Subset datasets to only include EPHA3+ cells (counts>0) #####################

# Subset to EPHA3+ cells
obj_Ebert_EPHA3pos <- subset(obj_Ebert, subset = EPHA3 > 0)
obj_LeBlanc_EPHA3pos <- subset(obj_LeBlanc, subset = EPHA3 > 0)
obj_Abdelfattah_EPHA3pos <- subset(obj_Abdelfattah, subset = EPHA3 > 0)

########### calculate percentages per cell type (EPHA3+) #####################

# Function to calculate percentages per cell type
get_celltype_percent <- function(seurat_obj) {
  celltype_counts <- table(seurat_obj$Pooled_celltype)
  celltype_percent <- prop.table(celltype_counts) * 100
  return(as.data.frame(celltype_percent))
}

# Calculate for each dataset
pct_Ebert <- get_celltype_percent(obj_Ebert_EPHA3pos)
pct_LeBlanc <- get_celltype_percent(obj_LeBlanc_EPHA3pos)
pct_Abdelfattah <- get_celltype_percent(obj_Abdelfattah_EPHA3pos)

pct_Abdelfattah

#plot pie chart
library(ggplot2)

plot_pie_chart <- function(data, title) {
  ggplot(data, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y") +
    theme_void() +
    labs(title = title, fill = "Cell Type") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# Add meaningful titles
pie_Ebert <- plot_pie_chart(pct_Ebert, "EPHA3⁺ Cell Types - Ebert")
pie_LeBlanc <- plot_pie_chart(pct_LeBlanc, "EPHA3⁺ Cell Types - LeBlanc")
pie_Abdelfattah <- plot_pie_chart(pct_Abdelfattah, "EPHA3⁺ Cell Types - Abdelfattah")

# Print the plots
print(pie_Ebert)
print(pie_LeBlanc)
print(pie_Abdelfattah)

total_Ebert <- ncol(obj_Ebert_EPHA3pos)
total_LeBlanc <- ncol(obj_LeBlanc_EPHA3pos)
total_Abdelfattah <- ncol(obj_Abdelfattah_EPHA3pos)
