## Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr) 
library(RColorBrewer)

## Creat output directory 
out_dir <- "~/ChinmoST/output/scRNA/seurat/subclustering/germline_subclustering/"
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


# subset germline cells
germline_cells <- rownames(testis_scRNA.integrated@meta.data)[which(testis_scRNA.integrated@meta.data$lineage=="Germline cells")]
germline_cells.scRNA <- subset(testis_scRNA.integrated,cells=germline_cells)

# split the combined object into a list, with each dataset as an element
germline_cells.scRNA.ls <- SplitObject(germline_cells.scRNA,split.by = "ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(germline_cells.scRNA.ls)){
  germline_cells.scRNA.ls[[i]] <- NormalizeData(germline_cells.scRNA.ls[[i]],verbose = FALSE)
  germline_cells.scRNA.ls[[i]] <- FindVariableFeatures(germline_cells.scRNA.ls[[i]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# identify anchors using the FindIntegrationAnchors function
germline_cells.scRNA.anchors <- FindIntegrationAnchors(object.list = germline_cells.scRNA.ls,anchor.features = 2000,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
germline_cells.scRNA.integrated <- IntegrateData(anchorset = germline_cells.scRNA.anchors, dims = 1:30,features.to.integrate = rownames(germline_cells.scRNA.ls[[1]]))

# switch to integrated assay
DefaultAssay(germline_cells.scRNA.integrated) <- "integrated"

# scale and center features in the dataset
germline_cells.scRNA.integrated <- ScaleData(germline_cells.scRNA.integrated, features =rownames(germline_cells.scRNA.integrated))

# Perform linear dimensional reduction
germline_cells.scRNA.integrated <- RunPCA(germline_cells.scRNA.integrated, npcs = 50, verbose = FALSE)

# Determine the ‘dimensionality’ of the dataset
# JackStraw 
germline_cells.scRNA.integrated <- JackStraw(germline_cells.scRNA.integrated, num.replicate = 100, dims =50)
germline_cells.scRNA.integrated  <- ScoreJackStraw(germline_cells.scRNA.integrated, dims = 1:50)
pdf(str_c(out_dir,"integrated_pc.pdf"))
JackStrawPlot(germline_cells.scRNA.integrated, dims = 1:50)
ElbowPlot(germline_cells.scRNA.integrated,ndims=50)
dev.off()

# tune the number of PCs
FlyCellAtlas_markers <- c(
  "esg", "nos", "stg", # spermatogonium
  "bam", "nos", "stai", # mid-late proliferating spermatogonia
  "Rbp4", "bam", "stg", "vas",# spermatogonium-spermatocyte transition
  "Rbp4", "bam", # spermatocyte 1
  "Rbp4", # spermatocyte 3,4
  "CycB", "kl-3", # spermatocyte 5
  "CycB",# spermatocyte 6
  "CG3927", "CycB", "aly", # maturing primary spermatocyte
  "kl-2", "kl-3", "kl-5", # late primary spermatocyte
  "sunz", "whip", # early-mid elongation-stage spermatid
  "c-cup", "orb", "p-cup", "soti", "whip"# mid-late elongation-stage spermatid)
)

our_curated_markers <- c(
  "aub","zpg","vas","nos","stg","bam", # GSCs/spermatogonia
  "Rbp4","kmg","wuc","nht","Taf12L", # early spermatocytes
  "nht","Taf12L","sa","CycB","fzo","twe","CG3927", # intermediate spermatocytes
  "CycB","fzo","CG3927","bol", # late spermatocytes
  "bol","Dpy-30L2", # early spermatids
  "schuy","hale","boly","whip","p-cup","wa-cup","r-cup","d-cup" # late spermatids
)

for ( nPCs in c(10,15,20)){
  DefaultAssay(germline_cells.scRNA.integrated) <- "integrated"
  germline_cells.scRNA.integrated <- FindNeighbors(germline_cells.scRNA.integrated, dims = 1:nPCs)
  germline_cells.scRNA.integrated <- FindClusters(germline_cells.scRNA.integrated, resolution = 0.2)
  germline_cells.scRNA.integrated <- RunUMAP(germline_cells.scRNA.integrated, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(germline_cells.scRNA.integrated, reduction = "umap",label=TRUE)
  p2 <- DimPlot(germline_cells.scRNA.integrated, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(germline_cells.scRNA.integrated) <- "RNA"
  pdf(str_c(out_dir,"integrated_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(germline_cells.scRNA.integrated,features=unique(c(FlyCellAtlas_markers,our_curated_markers)),combine=FALSE)
  print(p)
  dev.off()
}

# select 20 PCs
DefaultAssay(germline_cells.scRNA.integrated) <- "integrated"
germline_cells.scRNA.integrated <- FindNeighbors(germline_cells.scRNA.integrated, dims = 1:20)
germline_cells.scRNA.integrated <- RunUMAP(germline_cells.scRNA.integrated, dims = 1:20)


# tune the resolution 
for (i in seq(0.4,1.2,0.2)){
  germline_cells.scRNA.integrated <- FindClusters(germline_cells.scRNA.integrated, resolution = i)
  pdf(str_c(out_dir,"integrated_PC20_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(germline_cells.scRNA.integrated, reduction = "umap",label=TRUE)
  p2 <- DimPlot(germline_cells.scRNA.integrated, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}

# select 0.4 resolution
germline_cells.scRNA.integrated <- FindClusters(germline_cells.scRNA.integrated, resolution = 0.4)

# Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
germline_cells.scRNA.integrated <- BuildClusterTree(germline_cells.scRNA.integrated,dims =1:20,reorder=TRUE,reorder.numeric=TRUE)
pdf(str_c(out_dir,"integrated_20PCs_resolution0.4_ClusterTree.pdf"))
PlotClusterTree(germline_cells.scRNA.integrated,use.edge.length=FALSE)
dev.off()

newpalette <- colorRampPalette(brewer.pal(10,"Set3"))(length(unique(germline_cells.scRNA.integrated@meta.data$tree.ident)))
pdf(,out_dir"integrated_PC20_resolution0.4_reorder_UMAP.pdf")
DimPlot(germline_cells.scRNA.integrated, reduction = "umap",label=TRUE,group.by="tree.ident",cols=newpalette)
dev.off()

# annotate cell types
Idents(germline_cells.scRNA.integrated) <- germline_cells.scRNA.integrated$tree.ident
new.cluster.ids <- c("GSCs/spermatogonia","Late spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatocytes","Early spermatids","Early spermatids","Early spermatids","Early spermatids","Late spermatocytes","Late spermatids","Late spermatids","Early spermatocytes","GSCs/spermatogonia")
names(new.cluster.ids) <- levels(Idents(germline_cells.scRNA.integrated))
germline_cells.scRNA.integrated <- RenameIdents(germline_cells.scRNA.integrated, new.cluster.ids)
celltype_df <- data.frame(cell_type=as.vector(Idents(germline_cells.scRNA.integrated)),row.names=names(Idents(germline_cells.scRNA.integrated)))
germline_cells.scRNA.integrated <- AddMetaData(germline_cells.scRNA.integrated,celltype_df)


# Figure S4B - UMAP with cell types
Idents(germline_cells.scRNA.integrated) <- germline_cells.scRNA.integrated$cell_type
df <- FetchData(germline_cells.scRNA.integrated,c("UMAP_1","UMAP_2","cell_type"))
df$cell_type <- factor(df$cell_type,levels=c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids"))
newpalette <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF","#357EBDFF","#9632B8FF")
pdf(str_c(out_dir,"Germline_UMAP_with_cell_types.pdf"))
ggplot(df,aes(UMAP_1,UMAP_2,color=cell_type))+
  geom_point(size=0.2)+
  theme_bw() +
  scale_color_manual(values = newpalette)+
  labs(x="UMAP Dimension 1",y="UMAP Dimension 2",col="Cell type")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = 'none',panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=14))
dev.off()


# Figure 3A - UMAP split by sample 
Idents(germline_cells.scRNA.integrated) <- germline_cells.scRNA.integrated$cell_type
Idents(germline_cells.scRNA.integrated) <- factor(Idents(germline_cells.scRNA.integrated),levels=c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids"))
newpalette <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF","#357EBDFF","#9632B8FF")
pdf(str_c(out_dir,"Germline_UMAP_split_by_sample.pdf"),width=12)
DimPlot(germline_cells.scRNA.integrated, reduction = "umap",group.by="rough_cell_type",split.by="orig.ident",cols=newpalette,ncol=3,label=FALSE,pt.size=0.4)
dev.off()


# Figure 3B - percentage of each sample in each cell type
df <- data.frame(sample=germline_cells.scRNA.integrated$orig.ident,cell_type=germline_cells.scRNA.integrated$cell_type)
counts_df <- df %>% 
  group_by(sample,cell_type) %>%
  summarise(counts=n()) %>%
  mutate(normalized_counts = counts / sum(counts))
counts_df <- counts_df %>%
  group_by(cell_type) %>%
  mutate(sample_percentage_each_cell_type = normalized_counts / sum(normalized_counts)*100)
counts_df$cell_type <- factor(counts_df$cell_type,levels=c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids"))
newpalette <- c(brewer.pal(9, "Blues")[c(2,4,6)],brewer.pal(9, "YlOrRd")[c(2,4,6)])
pdf(out_dir,"germline_lineage_sample_percentage_each_cell_type.pdf",width=6.5,height=6.5)
ggplot(data=counts_df,aes(x=cell_type,y=sample_percentage_each_cell_type,fill=sample))+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  scale_fill_manual(values=newpalette)+
  guides(col = guide_legend(ncol = 3, byrow = TRUE)) +
  labs(x="Cell type",y="Percentage ( % )",fill="Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color="black",size=12),axis.text.y=element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"))
dev.off()



## PCA of FlyCellAtlas WT ovary (10x+ss2), WT testis and Mutant testis germline pseudobulk 
library(SeuratDisk)
FlyCellAtlas_10x_ss2_ovary.seurat <- LoadH5Seurat("~/ChinmoST/input/FlyCellAtlas/10x_ss2_ovary.h5seurat")
FlyCellAtlas_10x_ss2_ovary_data <- log1p(FlyCellAtlas_10x_ss2_ovary.seurat@assays$RNA@counts)
FlyCellAtlas_10x_ss2_ovary.seurat@assays$RNA@data <- FlyCellAtlas_10x_ss2_ovary_data
FlyCellAtlas_10x_ss2_ovary.seurat$lineage <- ifelse(grepl("germ",FlyCellAtlas_10x_ss2_ovary.seurat$annotation),"Germline","Somatic")
Idents(FlyCellAtlas_10x_ss2_ovary.seurat) <- FlyCellAtlas_10x_ss2_ovary.seurat$lineage
FlyCellAtlas_10x_ss2_ovary_df <- AverageExpression(FlyCellAtlas_10x_ss2_ovary.seurat,assays="RNA")$RNA
FlyCellAtlas_10x_ss2_ovary_df <- FlyCellAtlas_10x_ss2_ovary_df %>%
  dplyr::select(Germline) %>%
  dplyr::rename(WT_ovary=Germline)


germline_cells.scRNA.integrated$lineage <- "Germline"
Idents(germline_cells.scRNA.integrated) <- paste(germline_cells.scRNA.integrated$orig.ident,germline_cells.scRNA.integrated$lineage,sep="_")
MaLab_testis_df <- AverageExpression(germline_cells.scRNA.integrated,assays="RNA")$RNA
merged_df <- merge(MaLab_testis_df,FlyCellAtlas_10x_ss2_ovary_df,by="row.names")
merged_df <- merged_df[,-1]
pca_matrix <- t(merged_df)
pca <- prcomp(pca_matrix)

# variance explained by each PC
pc_eigenvalues <- pca$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
pdf(str_c(out_dir,"testis_FlyCellAtlas_ovary_germline_PCA_scree_plot.pdf"))
pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")
dev.off()

# Figure S3C - Visualising samples on 2D PC space
pc_scores <- data.frame(sample=rownames(pca$x),PC1=pca$x[,1],PC2=pca$x[,2])
pc_scores$sample <- factor(pc_scores$sample,levels=c("CT35_Germline","CT68_Germline","CT911_Germline","Mu35_Germline","Mu68_Germline","Mu911_Germline","WT_ovary"),labels=c("WT testis 3-5D","WT testis 6-8D","WT testis 9-11D","Mutant testis 3-5D","Mutant testis 6-8D","Mutant testis 9-11D","WT ovary"))
newpalette <- c(brewer.pal(9, "Blues")[c(2,4,6)],brewer.pal(9, "YlOrRd")[c(2,4,6)],brewer.pal(9, "Greens")[5])
pdf(str_c(out_dir,"testis_FlyCellAtlas_ovary_germline_PCA_plot.pdf"),width=7,height=5)
ggplot(data=pc_scores,aes(x=PC1,y=PC2,color=sample))+
  geom_point(size=5)+
  scale_color_manual(values=newpalette)+
  theme_bw()+
  xlim(-125,350)+
  ylim(-125,350)+
  labs(x="PC1 (79.9%)",y="PC2 (14.0%)",color="Sample")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),axis.text = element_text(size=13),axis.title=element_text(size=14))
dev.off()


## Correlation of pseudubulk transcriptomes for each germline cell type of WT,Mutant testis and FlyCellAtlas WT ovary germ cells
library(SeuratDisk)
FlyCellAtlas_10x_ss2_ovary.seurat <- LoadH5Seurat("~/ChinmoST/input/FlyCellAtlas/10x_ss2_ovary.h5seurat")
FlyCellAtlas_10x_ss2_ovary_data <- log1p(FlyCellAtlas_10x_ss2_ovary.seurat@assays$RNA@counts)
FlyCellAtlas_10x_ss2_ovary.seurat@assays$RNA@data <- FlyCellAtlas_10x_ss2_ovary_data
FlyCellAtlas_germline.seurat <- subset(FlyCellAtlas_10x_ss2_ovary.seurat,subset=annotation %in% c("young germ cell","post-mitotic germ cell early 16-cell cyst","16-cell germline cyst in germarium region 2a and 2b","germ cell stage 4 and later","germline cell, unknown stage"))
Idents(FlyCellAtlas_germline.seurat) <- FlyCellAtlas_germline.seurat$annotation


# retain shared genes
shared_genes <- intersect(rownames(FlyCellAtlas_germline.seurat),rownames(germline_cells.scRNA.integrated))


# pseudobulk transcriptomes for each  germ cell type of WT ovary 
FlyCellAtlas_germline.seurat.ls <- SplitObject(FlyCellAtlas_germline.seurat,split.by="ident")
FlyCellAtlas_germline_pseudobulk.ls <- lapply(1:length(FlyCellAtlas_germline.seurat.ls),function(i){
  normalized_counts <- as.data.frame(FlyCellAtlas_germline.seurat.ls[[i]]@assays$RNA@counts)
  normalized_pseudocounts <- apply(normalized_counts,1,mean)
  df <- data.frame(gene=names(normalized_pseudocounts),counts=as.vector(normalized_pseudocounts))
  df <- df[which(df$gene %in% shared_genes),]
  df
  })
FlyCellAtlas_germline_pseudobulk_df <- do.call(cbind,FlyCellAtlas_germline_pseudobulk.ls)
FlyCellAtlas_germline_pseudobulk_df <- FlyCellAtlas_germline_pseudobulk_df[,c(1,seq(2,10,2))]
colnames(FlyCellAtlas_germline_pseudobulk_df)[2:6] <- paste("ovary",names(FlyCellAtlas_germline.seurat.ls),sep="_")


# pseudobulk transcriptomes for each germline cell type of WT,Mutant testis
germline_cells.scRNA.integrated$sample_type <-  ifelse(germline_cells.scRNA.integrated$orig.ident %in% c("CT35","CT68","CT911"),"WT","Mutant")
Idents(germline_cells.scRNA.integrated) <- paste(germline_cells.scRNA.integrated$sample_type,germline_cells.scRNA.integrated$cell_type,sep="_")
germline_cells_scRNA.ls <- SplitObject(germline_cells.scRNA.integrated,split.by="ident")
germline_cells_pseudobulk.ls <- lapply(1:length(germline_cells_scRNA.ls),function(i){
  counts <- as.data.frame(germline_cells_scRNA.ls[[i]]@assays$RNA@counts)
  normalized_counts <- apply(counts,2,function(x) x/sum(as.numeric(x))*10^6)
  normalized_pseudocounts <- apply(normalized_counts,1,mean)
  df <- data.frame(gene=names(normalized_pseudocounts),counts=as.vector(normalized_pseudocounts))
  df <- df[which(df$gene %in% shared_genes),]
  df
  })
germline_cells_pseudobulk_df <- do.call(cbind,germline_cells_pseudobulk.ls)
germline_cells_pseudobulk_df <- germline_cells_pseudobulk_df[,c(1,seq(2, 24,2))]
colnames(germline_cells_pseudobulk_df)[2:13] <-  names(germline_cells_scRNA.ls)


# calculate Pearson correlation coefficients
df <- merge(germline_cells_pseudobulk_df,FlyCellAtlas_germline_pseudobulk_df,by="gene")
df <- df[,-1]
df <- df[,c(paste(rep(c("WT","Mutant"),6),rep(c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids"),each=2),sep="_"),paste("ovary",c("young germ cell","post-mitotic germ cell early 16-cell cyst","16-cell germline cyst in germarium region 2a and 2b","germ cell stage 4 and later"),sep="_"))]
correlationMatrix <- cor(df)


# Figure 3C, S4C
library(corrplot)
pdf(str_c(out_dir,"testis_ovary_germline_cells_corrplot.pdf"),width=12,height=12)
corrplot(correlationMatrix,tl.col="black",tl.cex=1, cl.lim =c(0,1),cl.length=5,method="pie")
dev.off()


## ChinmoST vs WT DEGs for each cell type and their enriched GO terms
library(clusterProfiler)
library(org.Dm.eg.db)

Idents(germline_cells.scRNA.integrated) <- paste(germline_cells.scRNA.integrated$sample_type,germline_cells.scRNA.integrated$cell_type,sep="_")
cell_types <- c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids")
WT_cell_types <- paste("WT",cell_types,sep="_")
Mutant_cell_types <- paste("Mutant",cell_types,sep="_")
DefaultAssay(germline_cells.scRNA.integrated) <- "RNA"
germline_Mutant_vs_WT_markers.ls <- lapply(1:length(cell_types),function(i){
	FindMarkers(germline_cells.scRNA.integrated, ident.1 = Mutant_cell_types[i], ident.2 = WT_cell_types[i],assay="RNA",logfc.threshold=-Inf,min.pct=-Inf)
	})
germline_mutant_up_genes.ls <- lapply(1:length(cell_types),function(i){
  df <- germline_Mutant_vs_WT_markers.ls[[i]] %>%
  	filter(p_val<0.05,avg_logFC>=log(1.5)) %>%
  	arrange(desc(avg_logFC),desc(pct.1))
  df
})
names(germline_mutant_up_genes.ls) <- cell_types
germline_mutant_down_genes.ls <- lapply(1:length(cell_types),function(i){
  df <- germline_Mutant_vs_WT_markers.ls[[i]] %>%
  	filter(p_val<0.05,avg_logFC<=-log(1.5)) %>%
  	arrange(avg_logFC,desc(pct.1))
  df
})
names(germline_mutant_down_genes.ls) <- cell_types


get_logFC_with_ENTREZID_id <- function(DEGs_df){
  DEGs <- rownames(DEGs_df)
  DEGs <- bitr(DEGs,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Dm.eg.db",drop=FALSE)
  DE_logFC <- DEGs_df$avg_logFC
  if (any(is.na(DEGs$ENTREZID))==TRUE){
    DE_logFC <- DE_logFC[-which(is.na(DEGs$ENTREZID))]
    names(DE_logFC) <- DEGs$ENTREZID[-which(is.na(DEGs$ENTREZID))]
  } else {
    names(DE_logFC) <- DEGs$ENTREZID
  }
  DE_logFC <- sort(DE_logFC,decreasing =TRUE)
  return(DE_logFC)
}

GO_ORA <- function(logFC_with_ENTREZID_id,kept_term_number=10){
  BP_GO_ORA <- enrichGO(gene=names(logFC_with_ENTREZID_id),
    OrgDb=org.Dm.eg.db,
    ont= "BP",
    readable= TRUE)
  BP_GO_ORA <- simplify(BP_GO_ORA)
  MF_GO_ORA <- enrichGO(gene=names(logFC_with_ENTREZID_id),
    OrgDb=org.Dm.eg.db,
    ont= "MF",
    readable= TRUE)
  MF_GO_ORA <- simplify(MF_GO_ORA)
  CC_GO_ORA <- enrichGO(gene=names(logFC_with_ENTREZID_id),
    OrgDb=org.Dm.eg.db,
    ont= "CC",
    readable= TRUE)
  CC_GO_ORA <- simplify(CC_GO_ORA)
  BP_GO_ORA_result <- as.data.frame(BP_GO_ORA)
  MF_GO_ORA_result <- as.data.frame(MF_GO_ORA)
  CC_GO_ORA_result <- as.data.frame(CC_GO_ORA)
  BP_GO_ORA_result$Type <- "Biological process"
  MF_GO_ORA_result$Type <- "Molecular function"
  CC_GO_ORA_result$Type <- "Cellular component"
  BP_GO_ORA_result <- BP_GO_ORA_result[order(BP_GO_ORA_result$pvalue),]
  MF_GO_ORA_result <- MF_GO_ORA_result[order(MF_GO_ORA_result$pvalue),]
  CC_GO_ORA_result <- CC_GO_ORA_result[order(CC_GO_ORA_result$pvalue),]
  results <- list(BP_GO_ORA=BP_GO_ORA_result,
    MF_GO_ORA=MF_GO_ORA_result,
    CC_GO_ORA=CC_GO_ORA_result,
    plot_df=plot_df)
  return(results)
}


# Figure 3F - DEGs for ChinmoST vs WT GSCs/spermatogonia and their enriched terms
Mutant_GSCs_up_DEGs <- get_logFC_with_ENTREZID_id(germline_mutant_up_genes.ls[["GSCs/spermatogonia"]])
Mutant_GSCs_up_DEGs_GO_ORA <- GO_ORA(Mutant_GSCs_up_DEGs)

Mutant_GSCs_down_DEGs <- get_logFC_with_ENTREZID_id(germline_mutant_down_genes.ls[["GSCs/spermatogonia"]])
Mutant_GSCs_down_DEGs_GO_ORA <- GO_ORA(Mutant_GSCs_down_DEGs)


# kegg
library(KEGGREST)
get_keggConv_geneList <- function(logFC_with_ENTREZID_id){
  genes <- paste("ncbi-geneid",names(logFC_with_ENTREZID_id),sep=":")
  kegg_genes <- keggConv("genes",genes)
  gsub("dme:","",kegg_genes)
}
Mutant_GSCs_up_DEGs_keggConv_geneList <- get_keggConv_geneList(Mutant_GSCs_up_DEGs)
Mutant_GSCs_up_DEGs_KEGG_ORA <- enrichKEGG(
	gene=Mutant_GSCs_up_DEGs_keggConv_geneList,
    organism     = 'dme',
    pvalueCutoff = 0.05)
Mutant_GSCs_up_DEGs_KEGG_ORA@result[which(Mutant_GSCs_up_DEGs_KEGG_ORA@result$pvalue < 0.05),]

Mutant_GSCs_down_DEGs_keggConv_geneList <- get_keggConv_geneList(Mutant_GSCs_down_DEGs)
Mutant_GSCs_down_DEGs_KEGG_ORA <- enrichKEGG(
	gene= Mutant_GSCs_down_DEGs_keggConv_geneList,
    organism     = 'dme',
    pvalueCutoff = 0.05)
Mutant_GSCs_down_DEGs_KEGG_ORA@result[which(Mutant_GSCs_down_DEGs_KEGG_ORA@result$pvalue < 0.05),]


## Figure S3G - Volcano plot visualzing ChinmoST vs WT DEGs 
pdf(str_c(out_dir,"mutant_vs_WT_DEGs_volcano_plot.pdf"),width=9)
for ( i in 1:length(cell_types)){
  germline_Mutant_vs_WT_markers.ls[[i]]$mlog10p_val <- -log10(germline_Mutant_vs_WT_markers.ls[[i]]$p_val)
  germline_Mutant_vs_WT_markers.ls[[i]]$type <- "nonsignificant DEGs"
  germline_Mutant_vs_WT_markers.ls[[i]]$type[which(germline_Mutant_vs_WT_markers.ls[[i]]$p_val<0.05 & germline_Mutant_vs_WT_markers.ls[[i]]$avg_logFC>=log(1.5))] <- "up-regulated DEGs in mutant"
  germline_Mutant_vs_WT_markers.ls[[i]]$type[which(germline_Mutant_vs_WT_markers.ls[[i]]$p_val<0.05 & germline_Mutant_vs_WT_markers.ls[[i]]$avg_logFC<=-log(1.5)] <- "down-regulated DEGs in mutant"
  germline_Mutant_vs_WT_markers.ls[[i]]$type <- factor(germline_Mutant_vs_WT_markers.ls[[i]]$type,levels=c("up-regulated DEGs in mutant","down-regulated DEGs in mutant","nonsignificant DEGs"))
  top10 <- c(rownames(germline_mutant_up_genes.ls[[i]])[1:5],rownames(germline_mutant_down_genes.ls[[i]])[1:5])
  germline_Mutant_vs_WT_markers.ls[[i]]$label <- ifelse(rownames(germline_Mutant_vs_WT_markers.ls[[i]]) %in% top10,rownames(germline_Mutant_vs_WT_markers.ls[[i]]),"")
  germline_Mutant_vs_WT_markers.ls[[i]$log2FC <- log2(exp(germline_Mutant_vs_WT_markers.ls[[i]$avg_logFC))
  p1 <- ggplot(data=germline_Mutant_vs_WT_markers.ls[[i]],aes(x=log2FC,y=mlog10p_val))+
    geom_point(aes(color=type))+
    ggrepel::geom_text_repel(data=germline_Mutant_vs_WT_markers.ls[[i]],aes(x=log2FC,y=mlog10p_val,label=label))+
    xlim(-max(abs(germline_Mutant_vs_WT_markers.ls[[i]]$log2FC)),max(abs(germline_Mutant_vs_WT_markers.ls[[i]]$log2FC)))+
    scale_color_manual(values = c("firebrick3","dodgerblue3","lightgrey"))+
    geom_vline(xintercept=c(-log2(1.5),log2(1.5)),linetype="dotted")+
    theme_classic()+
    theme_bw()+
    labs(x="Log2 Fold Change",y="-log10(p_value)",title=str_c("Mutant vs WT ",cell_types[i]," DEGs"))+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=12),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),legend.position ="none")
  print(p1)
}
dev.off()


## Figure 3G - Violin plot showing expression of interested genes in ChinmoST & WT GSCs/spermatogonia, Early spermatocytes
early_germline.scRNA <- subset(germline_cells.scRNA.integrated,subset=cell_type %in% c("GSCs/spermatogonia","Early spermatocytes"))
early_germline.scRNA$day <- gsub("CT|Mu","D",early_germline.scRNA$orig.ident)
early_germline.scRNA$sample_type_day <- paste(early_germline.scRNA$sample_type,early_germline.scRNA$day,sep="_")
genes <- c("wuc","ovo","otu","Hrb27C","eIF4A")
newpalette <- c(brewer.pal(9, "Blues")[c(2,4,6)],brewer.pal(9, "YlOrRd")[c(2,4,6)])
pdf(str_c(out_dir,"germline_development_interested_genes_time_series_expression_VlnPlot.pdf"),height=5,width=6)
for (i in 1:length(genes)){
  df <- FetchData(object=early_germline.scRNA,c(genes[i],"cell_type","day","sample_type","sample_type_day"),slot="data")
  colnames(df)[1] <- "gene" 
  df$sample_type <- factor(df$sample_type,levels=c("WT","Mutant"))
  df$cell_type <- factor(df$cell_type,levels=c("GSCs/spermatogonia","Early spermatocytes"))
  df$sample_type_day <- factor(df$sample_type_day,levels=paste(rep(c("WT","Mutant"),each=3),rep(c("D35","D68","D911"),2),sep="_"))
  p <- ggplot(df,aes(x=sample_type,y=gene,fill=sample_type_day))+
        geom_violin(scale="width",trim=TRUE,position=position_dodge(0.9))+
        geom_boxplot(width=0.1,outlier.size=0.8,position=position_dodge(0.9)) +
        scale_fill_manual(values=newpalette)+
        ylim(0,round(max(df[,1]),1))+
        facet_grid(day ~ .)+
        facet_grid(rough_cell_type ~ .,switch = "y")+
        labs(x="Sample",y="Expression Level",fill="Sample type",title=genes[i]) +
        theme_bw() +
        theme(plot.title = element_text(face = "bold",size = rel(1.5), hjust = 0.5),
        text = element_text(),
        panel.background=element_rect(fill='transparent', color='black',linetype="solid"),
        axis.title = element_text(face = "bold",size = rel(1.2)),
        axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(color="black",size=rel(1.2)), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(face="italic"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="black",fill="#f0f0f0",linetype="solid"),
        strip.text = element_text(face="bold",size=10,margin = margin(1.5,0.3,1.5,0.3, "cm")),
        legend.position = "none")
    print(p)
  
}
dev.off()