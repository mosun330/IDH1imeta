
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(MAST)
library(tidyverse)
library(reshape2)
library(ggalluvial)
library(scales)
library(ggpubr)
library(rstatix)
library(CellChat)
library(ggsci)
library(VennDiagram)
library(ggsignif)

#1.Integration

seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(colon)
  colon <- ScaleData(colon, features = VariableFeatures(colon))
  colon <- RunPCA(colon, features = VariableFeatures(object = colon))
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:20)
  return(colon)
}

read_10x_data <- function(data_dir,geneuse=2) {
  tryCatch({
    Read10X(data.dir = data_dir,gene.column = geneuse)
  }, error = function(e) {
    cat("Standard Read10X failed, trying manual reading...\n")
    files <- list.files(data_dir)
    matrix_file <- grep("matrix.mtx", files, value = TRUE)
    features_file <- grep("features.tsv|genes.tsv", files, value = TRUE)
    barcodes_file <- grep("barcodes.tsv", files, value = TRUE)
    
    if(length(matrix_file) == 0 || length(features_file) == 0 || length(barcodes_file) == 0) {
      stop(paste("Unable to find necessary files. Directory contents:", paste(files, collapse = ", ")))
    }
    mat <- Matrix::readMM(file.path(data_dir, matrix_file))
    features <- read.delim(file.path(data_dir, features_file), header = FALSE, stringsAsFactors = FALSE)
    barcodes <- read.delim(file.path(data_dir, barcodes_file), header = FALSE, stringsAsFactors = FALSE)[[1]]
    if(ncol(features) >= 2) {
      rownames(mat) <- features$V2
    } else {
      rownames(mat) <- features$V1
    }
    colnames(mat) <- barcodes
    
    return(mat)
  })
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name,geneuse = 2){
  # function for basic seurat based qc and doubletfinder based doublet removal
  setwd("/media/desk/")
  colon.data <- read_10x_data(data_directory,geneuse)
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
  currentSample$orig.ident <- factor(project_name)  # Explicitly set as factor
  currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^mt-")
  erythrocyte_genes <- c("Hba-a1", "Hbb-bs", "Hbb-bt", "Hbb-b1", "Hbb-b2", "Alas2")
  available_genes <- erythrocyte_genes[erythrocyte_genes %in% rownames(currentSample)]
  if(length(available_genes) > 0) {
    currentSample[["percent.ery"]] <- PercentageFeatureSet(currentSample, features = available_genes)
  } else {
    currentSample[["percent.ery"]] <- 0
    warning("No erythrocyte genes found in dataset: ", paste(erythrocyte_genes, collapse = ", "))
  }
  # qc plot-pre filtering
  setwd("/media/desk/")
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ery"),
                ncol = 4, pt.size = 0.05))
  dev.off()
  pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ery"),
                ncol = 4, pt.size = 0))
  dev.off()
  # filter everything to 400 unique genes/cell
  currentSample <- subset(currentSample,
                          subset = nFeature_RNA > 300 &
                            nFeature_RNA < 7000 &
                            nCount_RNA > 500 &
                            percent.mt < 5 &
                            percent.ery < 0.5)
  # ============================= #
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  # Run doublet finder
  nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  seu_colon <- doubletFinder(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  print(head(seu_colon@meta.data))
  # rename columns
  seu_colon$doublet.class <- seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
  seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
  pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
  seu_colon$pANN <- seu_colon[[pann]]
  seu_colon[[pann]] <- NULL
  # plot pre and post doublet finder results
  pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
  dev.off()
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
  dev.off()
  # Remove extra stuff and return filtered Seurat object
  seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
  cat("Processed sample:", project_name, "Cells:", ncol(seu_colon), "\n")
  return(seu_colon)
}

seurat_qc_plots <- function(colon, sample_name){
  # Make some basic qc plots
  pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
  dev.off()
}

set.seed(1)
setwd("/media/desk/")
data_directory=c("S1/","S2/","S3/","S4/","S5/","S6/")
project_name=c("s1","s2","s3","s4","s5","s6")
samples <- project_name
sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1],geneuse = 1)
seu_list <- sample1
for (i in seq_along(samples)){
  sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i],geneuse = 1)
  seu_list=merge(seu_list,sc.i)
}

#2.Annotation

scRNA_harmony <- readRDS("/media/desk/scRNA_integrated.rds")
DimPlot(scRNA_harmony , reduction = "umap",label = T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident')
DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
table(scRNA_harmony$orig.ident)  
Idents(scRNA_harmony)="seurat_clusters"
table(scRNA_harmony@meta.data$seurat_clusters)
table(scRNA_harmony@active.ident)
scs=subset(scRNA_harmony,downsample=500)
table(scs@active.ident)
scs <- JoinLayers(scs)
markers <- FindAllMarkers(scs, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(markers %>%
       group_by(cluster) %>%
       slice_max(n = 20, order_by = avg_log2FC))



#3.Cell Proportion

pB2_df <- table(
  scRNA_harmony@meta.data$celltype.1, 
  scRNA_harmony@meta.data$group.1
) %>% 
  melt()
colnames(pB2_df) <- c("Cluster", "Sample", "Number")

pB2_df <- table(
  scRNA_harmony@meta.data$celltype.4, 
  scRNA_harmony@meta.data$group.1
) %>% 
  melt()
colnames(pB2_df) <- c("celltype", "sample", "freq")

sample_totals <- aggregate(freq ~ sample, data = pB2_df, sum)
colnames(sample_totals) <- c("sample", "total")
pB2_df <- merge(pB2_df, sample_totals, by = "sample")
pB2_df$percent <- pB2_df$freq / pB2_df$total

colour <- c(
  "T_cells_NK" = "#D45A5A",        
  "B_cells" = "#E8A8A8",            
  "Macrophages/Monocytes" = "#F2C899", 
  "Neutrophils" = "#F9E0C7",        
  "Proliferating/cycling_Neutrophils" = "#F0C060", 
  "Hepatocytes_tumor" = "#D2B48C" , 
  "Endothelial_tumor" = "#7BB9E0", 
  "Endothelial_normal" = "#A3D0F0" 
  )

sample_order <- c("1", "2", "3")
pB2_df$sample <- factor(pB2_df$sample, levels = sample_order)
pB2_df$celltype <- factor(pB2_df$celltype, levels = names(colour))

pp <- ggplot(pB2_df, 
             aes(x = sample, y = percent, 
                 stratum = celltype, alluvium = celltype, 
                 fill = celltype)) +
  scale_fill_manual(values = colour) +
  scale_y_continuous(labels = percent_format(), 
                     breaks = seq(0, 1, by = 0.25),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_classic()

pB2 <- pp + 
  geom_col(width = 0.6, color = NA, size = 0.5) +
  geom_flow(width = 0.6, alpha = 0.22, knot.pos = 0.35,
            color = 'white', size = 0.5) + 
  geom_alluvium(width = 0.6, alpha = 1, knot.pos = 0.35,
                fill = NA, color = 'white', size = 0.5) +
  
labs(
  title = "Cell Type Proportion Distribution",
  x = "Experimental Group", 
  y = "Cell Proportion (%)",
  fill = "Cell Type"
) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_blank(),  # Remove x-axis line
    axis.line.y = element_line(color = "black", linewidth = 0.8),  # Thicken y-axis line
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.ticks.y = element_line(color = "black", linewidth = 0.8),  # Thicken y-axis ticks
    axis.ticks.length.y = unit(0.2, "cm"),  # Lengthen y-axis tick marks
    axis.text.x = element_text(
      size = 11, 
      color = "black",
      family = "Arial",
      angle = 45,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      size = 11, 
      color = "black",
      family = "Arial",
      margin = margin(r = 5)
    ),
    axis.title.x = element_text(
      size = 12, 
      color = "black",
      family = "Arial",
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 12, 
      color = "black", 
      family = "Arial",
      face = "bold",
      margin = margin(r = 15)
    ),
    legend.position = "right",
    legend.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 11, face = "bold", family = "Arial"),
    legend.key.size = unit(1, "lines"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    plot.title = element_text(
      size = 14, 
      face = "bold", 
      hjust = 0.5,
      family = "Arial",
      margin = margin(b = 20)
    ),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(clip = "off")

print(pB2)

#4.CGAS Expression Analysis

scRNA_harmony <- readRDS("/media/desk/scRNA_annotated.rds")
DimPlot(scRNA_harmony,label = T,reduction = "umap")
Idents(scRNA_harmony)="orig.ident"
table(scRNA_harmony$orig.ident)
scRNA_harmony=RenameIdents(scRNA_harmony,"s1"="WT",  "s2"="WT",
                           "s3"="R132H", "s4"="R132H", 
                           "s5"="R132H_Flc",  "s6"="R132H_Flc")
table(scRNA_harmony@active.ident)
scRNA_harmony@meta.data$group.1=scRNA_harmony@active.ident
DimPlot(scRNA_harmony,split.by = "orig.ident",group.by = "celltype.1", label = F,)
Idents(scRNA_harmony) <- "group.1"

exp_data <- FetchData(scRNA_harmony, vars = c("Cgas", "group.1"))
colnames(exp_data) <- c("expression", "group")
cgas_positive <- subset(scRNA_harmony, subset = Cgas > 0)
vln_positive <- VlnPlot(cgas_positive, 
                        features = "Cgas",
                        group.by = "group.1",
                        pt.size = 0.3,
                        cols = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  ggtitle("Cgas Positive Cells Expression Distribution") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(vln_positive)

cgas_positive <- subset(scRNA_harmony, subset = Cgas > 0)
positive_exp_data <- FetchData(cgas_positive, vars = c("Cgas", "group.1"))
colnames(positive_exp_data) <- c("expression", "group")

vln_with_stats <- ggplot(positive_exp_data, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.3, color = "black") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  labs(title = "Cgas Positive Cells Expression Distribution", 
       x = "Experimental Group", 
       y = "Expression Level") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        axis.text = element_text(size = 12))
print(vln_with_stats)

vln_with_stats <- vln_with_stats +
  stat_compare_means(
    comparisons = list(
      c("WT", "R132H"),
      c("R132H", "R132H_Flc")
    ),
    method = "wilcox.test",  
    label = "p.format",     
    tip.length = 0.01,
    bracket.size = 0.5,
    label.y = c(3.5, 3.8)  
  )
print(vln_with_stats)

#5. CELLCHAT

data<-readRDS("/media/desk/scRNA_annotated.rds")
meta <- data@meta.data
Control_meta <- subset(meta,group.1=="1")
Model_meta <- subset(meta,group.1=="2")
Flc_meta <- subset(meta,group.1=="3")

data_joined <- JoinLayers(data, assay = "RNA") 

Control_data <- as.matrix(GetAssayData(object = data_joined, assay = "RNA", layer = "data")[, rownames(Control_meta)])
Model_data <- as.matrix(GetAssayData(object = data_joined, assay = "RNA", layer = "data")[, rownames(Model_meta)])
Flc_data <- as.matrix(GetAssayData(object = data_joined, assay = "RNA", layer = "data")[, rownames(Flc_meta)])

cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "celltype.4")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 40) # Use multisession mode for parallel processing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) # Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) # Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ScRNA_Control_cellchat.RDS")

cellchat <- createCellChat(object = Model_data, meta = Model_meta, group.by = "celltype.4")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 40) # Use multisession mode for parallel processing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ScRNA_Model_cellchat.RDS")

cellchat <- createCellChat(object = Flc_data, meta = Flc_meta, group.by = "celltype.4")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 40) # Use multisession mode for parallel processing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ScRNA_Flc_cellchat.RDS")

control<-readRDS( "Control_cellchat.RDS")
model<-readRDS( "Model_cellchat.RDS") 
Flc<-readRDS( "Flc_cellchat.RDS")

object.list <- list(Control=control,Model=model,Flc=Flc)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

dev.new(width = 10, height = 8)
p <- netVisual_heatmap(cellchat, color.use = color)
plot(p)

#6.Cytotoxicity and exhaustion Score

scRNA_harmony <- readRDS("/media/desk15/T_NK_annotation.rds")

exhaustion_genes <- c("Pdcd1", "Layn",  "Lag3", "Cd244", 
                      "Ctla4", "Lilrb1", "Tigit", "Tox",  
                      "Entpd1", "Cd160")
cytotoxicity_genes <- c("Gzma",  "Ifng", "Ifng","Prf1","Ccl4")

existing_exhaustion <- exhaustion_genes[exhaustion_genes %in% rownames(scRNA_harmony)]
missing_exhaustion <- exhaustion_genes[!exhaustion_genes %in% rownames(scRNA_harmony)]
existing_cytotoxicity <- cytotoxicity_genes[cytotoxicity_genes %in% rownames(scRNA_harmony)]
missing_cytotoxicity <- cytotoxicity_genes[!cytotoxicity_genes %in% rownames(scRNA_harmony)]

scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])
scRNA_harmony <- AddModuleScore(
  scRNA_harmony,
  features = list(exhaustion = existing_exhaustion),
  name = "exhaustion"
)
scRNA_harmony <- AddModuleScore(
  scRNA_harmony,
  features = list(cytotoxicity = existing_cytotoxicity),
  name = "cytotoxicity"
)

colnames(scRNA_harmony@meta.data)[grep("exhaustion1", colnames(scRNA_harmony@meta.data))] <- "exhaustion_score"
colnames(scRNA_harmony@meta.data)[grep("cytotoxicity1", colnames(scRNA_harmony@meta.data))] <- "cytotoxicity_score"

metadata <- scRNA_harmony@meta.data
metadata$cell_barcode <- rownames(metadata)
score_data <- data.frame(
  cell_barcode = metadata$cell_barcode,
  celltype = metadata$celltype6,
  group = metadata$group.1,
  exhaustion_score = metadata$exhaustion_score,
  cytotoxicity_score = metadata$cytotoxicity_score
)
score_data <- score_data[!is.na(score_data$exhaustion_score) & !is.na(score_data$cytotoxicity_score), ]

p1 <- ggplot(score_data, aes(x = group, y = exhaustion_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Exhaustion Score Across Groups",
       x = "Group", 
       y = "Exhaustion Score") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  stat_compare_means(
    method = "kruskal.test",  # Non-parametric test, suitable for non-normal distribution
    label = "p.format",
    label.y = max(score_data$exhaustion_score) * 1.1
  )

p2 <- ggplot(score_data, aes(x = group, y = cytotoxicity_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Cytotoxicity Score Across Groups",
       x = "Group", 
       y = "Cytotoxicity Score") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y = max(score_data$cytotoxicity_score) * 1.1
  )

combined_plot <- p1 + p2 + plot_layout(ncol = 2)
print(combined_plot)

cd8_celltypes <- grep("CD8|cd8", unique(score_data$celltype), value = TRUE, ignore.case = TRUE)
print(cd8_celltypes)

if(length(cd8_celltypes) > 0) {
  cd8_data <- score_data[score_data$celltype %in% cd8_celltypes, ]
  
  p3 <- ggplot(cd8_data, aes(x = group, y = exhaustion_score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "CD8+ T Cell Exhaustion Score",
         x = "Group", 
         y = "Exhaustion Score",
         subtitle = paste("Cell types:", paste(cd8_celltypes, collapse = ", "))) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    stat_compare_means(
      comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")),
      method = "wilcox.test",
      label = "p.signif"
    )
  
  p4 <- ggplot(cd8_data, aes(x = group, y = cytotoxicity_score, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "CD8+ T Cell Cytotoxicity Score",
         x = "Group", 
         y = "Cytotoxicity Score",
         subtitle = paste("Cell types:", paste(cd8_celltypes, collapse = ", "))) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    stat_compare_means(
      comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")),
      method = "wilcox.test",
      label = "p.signif"
    )
  
  cd8_plot <- p3 + p4 + plot_layout(ncol = 2)
  print(cd8_plot)
}

group_scores <- score_data %>%
  group_by(group) %>%
  summarise(
    mean_exhaustion = mean(exhaustion_score, na.rm = TRUE),
    mean_cytotoxicity = mean(cytotoxicity_score, na.rm = TRUE),
    sd_exhaustion = sd(exhaustion_score, na.rm = TRUE),
    sd_cytotoxicity = sd(cytotoxicity_score, na.rm = TRUE)
  )

print(group_scores)

melted_scores <- melt(group_scores[, 1:3], id.vars = "group")

heatmap_plot <- ggplot(melted_scores, aes(x = variable, y = group, fill = value)) +
  geom_tile(color = "white", lwd = 1) +
  geom_text(aes(label = round(value, 3)), color = "white", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(melted_scores$value)) +
  labs(title = "Average Scores by Group",
       x = "Score Type", 
       y = "Group",
       fill = "Score") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  scale_x_discrete(labels = c("Exhaustion", "Cytotoxicity"))

print(heatmap_plot)

kruskal_exhaustion <- kruskal.test(exhaustion_score ~ group, data = score_data)
print(kruskal_exhaustion)
kruskal_cytotoxicity <- kruskal.test(cytotoxicity_score ~ group, data = score_data)
print(kruskal_cytotoxicity)

