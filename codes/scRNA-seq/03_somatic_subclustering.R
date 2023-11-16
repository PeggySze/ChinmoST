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


# add annotation
new.cluster.ids <- c("Hub cells & CySCs","Terminal epithelium","Late cyst cells","Somatic cells at base of testis",rep("Late cyst cells",3),rep("Intermediate cyst cells",5),"Early cyst cells","Intermediate cyst cells",rep("Late cyst cells",2),"Intermediate cyst cells")
names(new.cluster.ids) <- levels(Idents(integrated_somatic_cells.scRNA))
integrated_somatic_cells.scRNA <- RenameIdents(integrated_somatic_cells.scRNA, new.cluster.ids)
rough_cell_type_df <-  data.frame(row.names=names(Idents(integrated_somatic_cells.scRNA )),rough_cell_type=Idents(integrated_somatic_cells.scRNA ))
integrated_somatic_cells.scRNA <- AddMetaData(integrated_somatic_cells.scRNA, rough_cell_type_df)


## further subclustering analysis (exclude late cyst cells,Terminal epithelium,Somatic cells at base of testis)
DefaultAssay(integrated_somatic_cells.scRNA) <- "RNA" 
Idents(integrated_somatic_cells.scRNA) <- integrated_somatic_cells.scRNA$orig.ident
cyst_cells.scRNA <- subset(integrated_somatic_cells.scRNA,subset=rough_cell_type %in% c("Hub cells & CySCs","Early cyst cells","Intermediate cyst cells"))
cyst_cells.scRNA <- DietSeurat(
  object=cyst_cells.scRNA,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )

# split the combined object into a list, with each dataset as an element
cyst_cells.scRNA.ls <- SplitObject(cyst_cells.scRNA,split.by = "ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(cyst_cells.scRNA.ls)){
  cyst_cells.scRNA.ls[[i]] <- NormalizeData(cyst_cells.scRNA.ls[[i]],verbose = FALSE)
  cyst_cells.scRNA.ls[[i]] <- FindVariableFeatures(cyst_cells.scRNA.ls[[i]],selection.method = "vst", nfeatures = 1000, verbose = FALSE)
}

# identify anchors using the FindIntegrationAnchors function
cyst_cells.scRNA.anchors <- FindIntegrationAnchors(object.list = cyst_cells.scRNA.ls,anchor.features = 1000,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
cyst_cells.scRNA.integrated <- IntegrateData(anchorset = cyst_cells.scRNA.anchors, dims = 1:30,features.to.integrate = rownames(cyst_cells.scRNA.ls[[1]]))

# switch to integrated assay
DefaultAssay(cyst_cells.scRNA.integrated) <- "integrated"

# scale and center features in the dataset
cyst_cells.scRNA.integrated <- ScaleData(cyst_cells.scRNA.integrated, features =rownames(cyst_cells.scRNA.integrated))

# Perform linear dimensional reduction
cyst_cells.scRNA.integrated <- RunPCA(cyst_cells.scRNA.integrated, npcs = 50, verbose = FALSE)
pdf(str_c(out_dir,"exclude_late_cyst_cells_subclustering_pc.pdf"))
ElbowPlot(cyst_cells.scRNA.integrated,ndims=50)
dev.off()


# select 20 PCs and 0.8 resolution
DefaultAssay(cyst_cells.scRNA.integrated) <- "integrated"
cyst_cells.scRNA.integrated <- FindNeighbors(cyst_cells.scRNA.integrated, dims = 1:20)
cyst_cells.scRNA.integrated <- FindClusters(cyst_cells.scRNA.integrated, resolution = 0.8)
cyst_cells.scRNA.integrated <- RunUMAP(cyst_cells.scRNA.integrated, dims = 1:20)


# expression of markers in UMAP
markers <- c(
  "org-1","Fas3",# Hub cells
  "InR","aop","Socs36E", # hub cells & CySCs
  "zfh1","ptc","piwi","tj","puc","stg", # CySCs
  "Wnt4","CG31676", # early cyst cells 
  "sog","Six4","Pdcd4" # intermediate cyst cells
  )
DefaultAssay(cyst_cells.scRNA.integrated) <- "RNA"
pdf(str_c(out_dir,"exclude_late_cyst_cells_subclustering_PC20_markers_UMAP.pdf"))
FeaturePlot(cyst_cells.scRNA.integrated,features=markers,cols=c("lightgrey","red"),combine=FALSE,order=TRUE)
dev.off()

# Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
cyst_cells.scRNA.integrated <- BuildClusterTree(cyst_cells.scRNA.integrated,dims =1:20,reorder=TRUE,reorder.numeric=TRUE)


# rename tree id
renamed_tree.ident <- c(1,2,11,12,13,3,4,5,6,7,8,9,10)
names(renamed_tree.ident) <- 1:13
cyst_cells.scRNA.integrated <- RenameIdents(cyst_cells.scRNA.integrated,renamed_tree.ident)
cyst_cells.scRNA.integrated$renamed_tree.ident <- Idents(cyst_cells.scRNA.integrated)


# add cell type annotation
Idents(cyst_cells.scRNA.integrated) <- cyst_cells.scRNA.integrated$tree.ident
cell_types <- c("Hub cells","CySCs","ICCs 2","ICCs 2","ICCs 2","ECCs",rep("ICCs 1",7))
names(cell_types) <- 1:13
cyst_cells.scRNA.integrated <- RenameIdents(cyst_cells.scRNA.integrated,cell_types)
cyst_cells.scRNA.integrated$cell_types <- Idents(cyst_cells.scRNA.integrated)


## Figure S3C - exclude late cyst cells subclustering UMAP with clusters
Idents(cyst_cells.scRNA.integrated) <- cyst_cells.scRNA.integrated$renamed_tree.ident
Idents(cyst_cells.scRNA.integrated) <- factor(Idents(cyst_cells.scRNA.integrated),levels=1:13)
newpalette <- brewer.pal(12,"Set3")[1:12]
newpalette[2] <- brewer.pal(9,"Set1")[2] 
newpalette[9] <- brewer.pal(9,"Set1")[9]
newpalette[13] <- brewer.pal(8,"Set2")[2]
pdf(str_c(out_dir,"exclude_late_cyst_cells_cyst_cells_UMAP_with_clusters.pdf"),width=7.5)
DimPlot(cyst_cells.scRNA.integrated, reduction = "umap",label=TRUE,cols=newpalette,label.size=5)+theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"))+labs(x="UMAP Dimension 1",y="UMAP Dimension 2")
dev.off()


## Figure S3D - dotplot of markers for late cyst cells subclustering 
Idents(cyst_cells.scRNA.integrated) <- factor(Idents(cyst_cells.scRNA.integrated),levels=c(13:9,5:4,8:6,3:1))
DefaultAssay(cyst_cells.scRNA.integrated) <- "RNA"
markers <- c(
  "org-1","Fas3",# Hub cells
  "InR","aop","Socs36E", # hub cells & CySCs
  "zfh1","ptc","piwi","tj","stg", # CySCs
  "Wnt4","CG31676", # early cyst cells 
  "sog","Six4" # intermediate cyst cells  
  )
pdf(str_c(out_dir,"exclude_late_cyst_cells_cluster_markers_dotplot.pdf"),width=10,height=6)
DotPlot(cyst_cells.scRNA.integrated,features=markers,cols =c("lightblue","red"))+RotatedAxis()+theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"))
dev.off()


## Figure 2A - exclude late cyst cells subclustering UMAP split by sample
Idents(cyst_cells.scRNA.integrated) <- str_c("C",cyst_cells.scRNA.integrated$renamed_tree.ident,"-",cyst_cells.scRNA.integrated$cell_types)
Idents(cyst_cells.scRNA.integrated) <- factor(Idents(cyst_cells.scRNA.integrated),levels=str_c("C",1:13,"-",c("Hub cells","CySCs","ECCs",rep("ICCs 1",7),rep("ICCs 2",3))))
pdf(str_c(out_dir,"exclude_late_cyst_cells_subclustering_UMAP_split_by_sample.pdf"),width=12)
DimPlot(cyst_cells.scRNA.integrated, reduction = "umap",split.by="orig.ident",cols=newpalette,ncol=3,label=FALSE)
dev.off()


## Figure 2B - percentage of each sample in each cluster
df <- data.frame(sample=cyst_cells.scRNA.integrated$orig.ident,clsuter=str_c("C",cyst_cells.scRNA.integrated$renamed_tree.ident,"-",cyst_cells.scRNA.integrated$cell_types))
counts_df <- df %>% 
  dplyr::group_by(sample,clsuter) %>%
  dplyr::summarise(counts=n()) %>%
  dplyr::mutate(normalized_counts = counts / sum(counts))
counts_df <- counts_df %>%
  dplyr::group_by(clsuter) %>%
  dplyr::mutate(sample_percentage_each_cluster = normalized_counts / sum(normalized_counts)*100)
counts_df$clsuter <- factor(counts_df$clsuter,levels=str_c("C",1:13,"-",c("Hub cells","CySCs","ECCs",rep("ICCs 1",7),rep("ICCs 2",3))))
newpalette <- c(brewer.pal(9, "Blues")[c(2,4,6)],brewer.pal(9, "YlOrRd")[c(2,4,6)])
pdf(str_c(out_dir,"exclude_late_cyst_cells_subclustering_sample_percentage_each_cluster.pdf"),width=7,height=7)
ggplot(data=counts_df,aes(x=clsuter,y=sample_percentage_each_cluster,fill=sample))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  scale_fill_manual(values=newpalette)+
  guides(col = guide_legend(ncol = 3, byrow = TRUE)) +
  labs(x="Cluster",y="Percentage ( % )",fill="Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12,color="black"),axis.text.y=element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5,size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.title =element_text(size=13))
dev.off()


## Figure S3E - Correlation of each cluster of testis and each cell type of WT ovary 
Rust_2020_ovary_scRNA <- readRDS("~/ChinmoST/input/Drosophila_ovary_atlas_Rust_2020.rds")
new.ids <- c("germ cell","germ cell","MB","MB","pFC","MB","MB","hemocyte","stretch cell","MB","stretch cell","MB","MB","MB","TF","FSC","EC","polar cell","muscle cell","pFC","stalk cell","EC","GSC","border cell","EC","cap cell")
names(new.ids) <- unique(Idents(Rust_2020_ovary_scRNA))
Rust_2020_ovary_scRNA <- RenameIdents(Rust_2020_ovary_scRNA, new.ids)


shared_genes <- intersect(rownames(cyst_cells.scRNA.integrated),rownames(Rust_2020_ovary_scRNA))
Rust_2020_ovary_scRNA.ls <- SplitObject(Rust_2020_ovary_scRNA,split.by="ident")
ovary_pseudobulk.ls <- lapply(1:length(Rust_2020_ovary_scRNA.ls),function(i){
  counts <- as.data.frame(Rust_2020_ovary_scRNA.ls[[i]]@assays$RNA@counts)
  normalized_counts <- apply(counts,2,function(x) x/sum(as.numeric(x))*10^6)
  normalized_pseudocounts <- apply(normalized_counts,1,mean)
  df <- data.frame(gene=names(normalized_pseudocounts),counts=as.vector(normalized_pseudocounts))
  df <- df[which(df$gene %in% shared_genes),]
  df
  })
ovary_pseudobulk_df <- do.call(cbind,ovary_pseudobulk.ls)
ovary_pseudobulk_df <- ovary_pseudobulk_df[,c(1,seq(2, 28,2))]
colnames(ovary_pseudobulk_df)[2:15] <- names(Rust_2020_ovary_scRNA.ls)
ovary_pseudobulk_df <- ovary_pseudobulk_df[,c("gene","TF","cap cell","EC","FSC","pFC","MB","polar cell","stalk cell","stretch cell")]

cyst_cells.scRNA.integrated$sample_type <- ifelse(cyst_cells.scRNA.integrated$orig.ident %in% c("CT35","CT68","CT911"),"WT","Mutant")
Idents(cyst_cells.scRNA.integrated) <- paste(cyst_cells.scRNA.integrated$sample_type,str_c("C",cyst_cells.scRNA.integrated$renamed_tree.ident),sep="_")
cyst_cells_scRNA.ls <- SplitObject(cyst_cells.scRNA.integrated,split.by="ident")
cyst_cells_pseudobulk.ls <- lapply(1:length(cyst_cells_scRNA.ls),function(i){
  counts <- as.data.frame(cyst_cells_scRNA.ls[[i]]@assays$RNA@counts)
  normalized_counts <- apply(counts,2,function(x) x/sum(as.numeric(x))*10^6)
  normalized_pseudocounts <- apply(normalized_counts,1,mean)
  df <- data.frame(gene=names(normalized_pseudocounts),counts=as.vector(normalized_pseudocounts))
  df <- df[which(df$gene %in% shared_genes),]
  df
  })
cyst_cells_pseudobulk_df <- do.call(cbind,cyst_cells_pseudobulk.ls)
cyst_cells_pseudobulk_df <- cyst_cells_pseudobulk_df[,c(1,seq(2, 52,2))]
colnames(cyst_cells_pseudobulk_df)[2:27] <- names(cyst_cells_scRNA.ls)

df <- merge(cyst_cells_pseudobulk_df,ovary_pseudobulk_df,by="gene")
df <- df[,-1]
df <- df[,c("TF","cap cell","EC","FSC","pFC","MB","polar cell","stalk cell","stretch cell",paste(rep(c("WT","Mutant"),13),str_c("C",rep(1:13,each=2)),sep="_"))]
library(corrplot)
correlationMatrix <- cor(df)
pdf(str_c(out_dir,"exclude_late_cyst_cells_testis_ovary_cyst_cells_corrplot.pdf"))
corrplot(correlationMatrix,tl.col="black", tl.srt=45,tl.cex=0.9, cl.lim =c(0,1),cl.length=11,method="pie",type="lower")
dev.off()


## Figure 2C - Correlation of WT,Mutant CySCs and WT ovary FSCs
ovary_FSC_pseudobulk_df <- ovary_pseudobulk_df[,c("gene","FSC")]
CySCs_pseudobulk_df <- cyst_cells_pseudobulk_df[,c("gene","WT_2","Mutant_2")]
df <- merge(CySCs_pseudobulk_df,ovary_FSC_pseudobulk_df,by="gene")
df <- df[,-1]
colnames(df) <- c("WT CySCs","Mutant CySCs","WT FSCs")
correlationMatrix <- cor(df)
colmat <- colorRampPalette(brewer.pal(9,"Blues"))
pdf(str_c(out_dir,"WT_Mutant_CySCs_FSCs_corrplot.pdf"))
corrplot(correlationMatrix,order="hclust",tl.col="black", tl.srt=45,tl.cex=1, col.lim =c(0,1),cl.length=5)
dev.off()


## Figure 2E 
