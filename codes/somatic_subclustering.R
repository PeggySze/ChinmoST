## Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr) 
library(RColorBrewer)
library(xlsx)

## Creat output directory 
out_dir <- "~/ChinmoST/output/scRNA/seurat/subclustering/somatic_subclustering/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

# Load Seurat object 
testis_scRNA.integrated <- readRDS("~/ChinmoST/output/scRNA/seurat/integrated_analysis/integrated_PC30_testis_scRNA.rds")
DefaultAssay(testis_scRNA.integrated) <- "RNA"
testis_scRNA.integrated <- DietSeurat(
  object=testis_scRNA.integrated,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )

# Subset somatic cells
somatic_cells <- rownames(testis_scRNA.integrated@meta.data)[which(testis_scRNA.integrated@meta.data$lineage=="Somatic cells")]
somatic_cells.scRNA <- subset(testis_scRNA.integrated,cells=somatic_cells)

# split the combined object into a list, with each dataset as an element
somatic_cells.scRNA.ls <- SplitObject(somatic_cells.scRNA,split.by = "ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(somatic_cells.scRNA.ls)){
  somatic_cells.scRNA.ls[[i]] <- NormalizeData(somatic_cells.scRNA.ls[[i]],verbose = FALSE)
  somatic_cells.scRNA.ls[[i]] <- FindVariableFeatures(somatic_cells.scRNA.ls[[i]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# identify anchors using the FindIntegrationAnchors function
somatic_cells.scRNA.anchors <- FindIntegrationAnchors(object.list = somatic_cells.scRNA.ls,anchor.features = 2000,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
somatic_cells.scRNA.integrated <- IntegrateData(anchorset = somatic_cells.scRNA.anchors, dims = 1:30,features.to.integrate = rownames(somatic_cells.scRNA.ls[[1]]))

# switch to integrated assay
DefaultAssay(somatic_cells.scRNA.integrated) <- "integrated"

# scale and center features in the dataset
somatic_cells.scRNA.integrated <- ScaleData(somatic_cells.scRNA.integrated, features =rownames(somatic_cells.scRNA.integrated))

# Perform linear dimensional reduction
somatic_cells.scRNA.integrated <- RunPCA(somatic_cells.scRNA.integrated, npcs = 50, verbose = FALSE)


# Determine the ‘dimensionality’ of the dataset
# JackStraw 
somatic_cells.scRNA.integrated <- JackStraw(somatic_cells.scRNA.integrated, num.replicate = 100, dims =50)
somatic_cells.scRNA.integrated  <- ScoreJackStraw(somatic_cells.scRNA.integrated, dims = 1:50)
pdf(str_c(out_dir,"integrated_pc.pdf"))
JackStrawPlot(somatic_cells.scRNA.integrated, dims = 1:50)
ElbowPlot(somatic_cells.scRNA.integrated,ndims=50)
dev.off()

## tune the number of PCs
FlyCellAtlas_markers <- c(
  "CadN", "Fas3", "Socs36E", "hh", "magu", "ptc", "upd1", # germinal proliferation center hub
  "esg", "stg", "tj", "zfh1", # cyst stem cell
  "N", "piwi", "tj", "zfh1", # early cyst cell 1
  "esg", "eya", "so", "tj", "zfh1", # early cyst cell 2
  "eya", "rdo",  # spermatocyte cyst cell branch a
  "Gs2", "Gyg", "eya", # spermatocyte cyst cell branch b
  "tj", "so", "rdo", # cyst cell branch a,cyst cell branch b
  "so",  # spermatid associated cyst cell
  "CG18628", "Fas3", "MtnA", "abd-A", # head cyst cell
  "Nep2", "Nep4", # terminal epithelial cell of testis,somatic cell at base of testis
  "Nep5", "Nep1" # seminal vesicle & testis epithelia)
)

our_curated_markers <- c(
  "org-1","Fas3",# Hub cells
  "InR","aop","Socs36E", # hub cells & CySCs
  "zfh1","ptc","piwi","tj","puc","stg", # CySCs
  "Wnt4","CG31676", # early cyst cells 
  "sog","Six4","Pdcd4", # intermediate cyst cells, spermatocytes associated CCs
  "CG8665", "CG3376", # late cyst cells , Spermatid associated CCs  
  "Nep2","Nep4","mnd", # Terminal epithelial cells
  "Ubx"# "Somatic cells at base of testis
)

for ( nPCs in seq(20,35,5)){
  DefaultAssay(somatic_cells.scRNA.integrated) <- "integrated"
  somatic_cells.scRNA.integrated <- FindNeighbors(somatic_cells.scRNA.integrated, dims = 1:nPCs)
  somatic_cells.scRNA.integrated <- FindClusters(somatic_cells.scRNA.integrated, resolution = 0.2)
  somatic_cells.scRNA.integrated <- RunUMAP(somatic_cells.scRNA.integrated, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(somatic_cells.scRNA.integrated, reduction = "umap",label=TRUE)
  p2 <- DimPlot(somatic_cells.scRNA.integrated, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(somatic_cells.scRNA.integrated) <- "RNA"
  pdf(str_c(out_dir,"integrated_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(somatic_cells.scRNA.integrated,features=unique(c(FlyCellAtlas_markers,our_curated_markers)),combine=FALSE)
  print(p)
  dev.off()
}

# select 30 PCs
DefaultAssay(somatic_cells.scRNA.integrated) <- "integrated"
somatic_cells.scRNA.integrated <- FindNeighbors(somatic_cells.scRNA.integrated, dims = 1:30)
somatic_cells.scRNA.integrated <- RunUMAP(somatic_cells.scRNA.integrated, dims = 1:30)

# tune the resolution 
for (i in seq(0.2,1.2,0.2)){
  somatic_cells.scRNA.integrated <- FindClusters(somatic_cells.scRNA.integrated, resolution = i)
  somatic_cells.scRNA.integrated <- RunUMAP(somatic_cells.scRNA.integrated, dims = 1:30)
  pdf(str_c(out_dir,"integrated_PC30_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(somatic_cells.scRNA.integrated, reduction = "umap",label=TRUE)
  p2 <- DimPlot(somatic_cells.scRNA.integrated, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}

# select resolution 0.6
somatic_cells.scRNA.integrated <- FindClusters(somatic_cells.scRNA.integrated, resolution = 0.6)

# Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
somatic_cells.scRNA.integrated <- BuildClusterTree(somatic_cells.scRNA.integrated,dims =1:30,reorder=TRUE,reorder.numeric=TRUE)

# Figure S3A
pdf(str_c(out_dir,"integrated_PC30_resolution0.6_reordered_UMAP.pdf"))
newpalette <- colorRampPalette(brewer.pal(10,"Set3"))(length(unique(Idents(somatic_cells.scRNA.integrated))))
DimPlot(somatic_cells.scRNA.integrated, reduction = "umap",label=TRUE,cols=newpalette)+theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"))
dev.off()


# Figure S3B - dotplot of cluster markers
Idents(integrated_somatic_cells.scRNA) <- factor(Idents(integrated_somatic_cells.scRNA),levels=c(4,2,16,15,7,6,5,3,17,14,12,11,10,9,8,13,1))
DefaultAssay(integrated_somatic_cells.scRNA) <- "RNA"
markers <- c(
  "org-1","Fas3",# Hub cells
  "InR","aop","Socs36E", # hub cells & CySCs
  "zfh1","ptc","piwi","tj","puc","stg", # CySCs
  "Wnt4","CG31676", # early cyst cells 
  "sog","Six4","Pdcd4", # intermediate cyst cells
  "CG8665", "CG3376", # late cyst cells , Spermatid associated CCs  
  "Nep2","Nep4","mnd", # Terminal epithelial cells
  "Ubx",# "Somatic cells at base of testis"
  "btl", "ect" # adult tracheocyte
  )
pdf(str_c(out_dir,"somatic_lineage_cluster_markers_dotplot.pdf"),width=10,height=6)
DotPlot(integrated_somatic_cells.scRNA,features=markers,cols =c("lightblue","red"))+RotatedAxis()+theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"))
dev.off()


# 




