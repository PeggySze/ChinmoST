## reference to condiments
# https://github.com/HectorRDB/condiments
# https://hectorrdb.github.io/condimentsPaper/


## Load required packages
library(Seurat)
library(condiments)
library(tidyr)
library(dplyr)
library(stringr) 
library(RColorBrewer)
library(ggplot2)
library(slingshot)
library(tradeSeq)
library(BiocParallel)
library(cowplot)
library(scales)


## Creat output directory 
out_dir <- "~/ChinmoST/output/scRNA/condiments/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load Seurat object 
germline_cells.scRNA.integrated <- readRDS("~/ChinmoST/output/scRNA/seurat/subclustering/germline_subclustering/integrated_germline_cells_scRNA.rds")
Idents(germline_cells.scRNA.integrated) <- germline_cells.scRNA.integrated$cell_type

## Exclude spermatids in the following analysis
exclude_spermatids_scRNA <- subset(germline_cells.scRNA.integrated,subset=cell_type %in% c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes"))
UMAP_df <- as.data.frame(exclude_spermatids_scRNA[["umap"]]@cell.embeddings)
Retained_cells <- rownames(UMAP_df)[which(UMAP_df$UMAP_1<=2.5)]
exclude_spermatids_scRNA <- subset(exclude_spermatids_scRNA,cells=Retained_cells)


## Covert Seurat object to SingleCellExperiment 
exclude_spermatids_scRNA.sce <-  as.SingleCellExperiment(exclude_spermatids_scRNA, assay = "RNA")

## calculate imbalance score
# Regions with a high score indicate that the local cell distribution according to treatment label is unbalanced compared the overall distribution. 
df <- bind_cols(
  as.data.frame(reducedDims(exclude_spermatids_scRNA.sce)$UMAP),
  as.data.frame(colData(exclude_spermatids_scRNA.sce))
  )
df$sample_type <- factor(df$sample_type,levels=c("WT","Mutant"))
scores <- imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = as.vector(df$sample_type))
df$scores <- scores$scaled_scores

# Figure S4E - UMAP with imbalance score
pdf(str_c(out_dir,"exclude_spermatids_germline_UMAP_with_imbalance_score.pdf"))
ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Imbalance scores")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=14))
dev.off()

## Trajectory Inference and Differential Topology
germline_cells_slingshot <- slingshot(exclude_spermatids_scRNA.sce, reducedDim = 'UMAP',clusterLabels = colData(exclude_spermatids_scRNA.sce)$cell_type,start.clus = "GSCs/spermatogonia", end.clus="Late spermatocytes",approx_points = 100)
set.seed(100)
topologyTest(SlingshotDataSet(germline_cells_slingshot), conditions=colData(exclude_spermatids_scRNA.sce)$sample_type, methods = "KS_mean", threshs = .01,rep=100)
#Running KS-mean test
#   method thresh statistic p.value
#1 KS_mean   0.01  0.083118 2.2e-16

## fit an individual trajectory for mutant and wild-type germline cells, respectively 
individual_slingshot <- slingshot_conditions(germline_cells_slingshot, 
  conditions=colData(exclude_spermatids_scRNA.sce)$sample_type, 
  approx_points = FALSE,extend = "n", reweight = FALSE, reassign = FALSE)
individual_slingshot$condition_id <- names(individual_slingshot)
individual_slingshot$mapping <- matrix(c(1,1), nrow = 1, ncol = 2, byrow = TRUE)
sds <- do.call(merge_sds, individual_slingshot)

df <- merge(df,as.data.frame(slingPseudotime(sds)),by="row.names")

# Figure S4E - UMAP with pseudotime and inferred trajectory curves
WT_pseudotime_df <- df[which(df$sample_type=="WT"),]
WT_slingCurves_df <- slingCurves(individual_slingshot[["WT"]])[[1]]$s[slingCurves(individual_slingshot[["WT"]])[[1]]$ord, ] %>% as.data.frame()
Mu_pseudotime_df <- df[which(df$sample_type=="Mutant"),]
Mu_slingCurves_df <- slingCurves(individual_slingshot[["Mutant"]])[[1]]$s[slingCurves(individual_slingshot[["Mutant"]])[[1]]$ord, ] %>% as.data.frame()
pdf(str_c(out_dir,"exclude_spermatids_germline_pseudotime_UMAP.pdf"),width=8)
ggplot() +
  geom_point(data=WT_pseudotime_df,aes(x = UMAP_1, y = UMAP_2, fill = Lineage1),size = 1,pch = 21, color = '#00000000') +
  scale_fill_viridis_c() +
  geom_path(data = WT_slingCurves_df, aes(x = UMAP_1, y = UMAP_2),size = 1,linetype = "dashed")  +
  labs(fill = "Pseudotime",title="WT")+
  ylim(-6,5.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=14))
ggplot() +
  geom_point(data=Mu_pseudotime_df,aes(x = UMAP_1, y = UMAP_2, fill = Lineage1),size = 1,pch = 21, color = '#00000000') +
  scale_fill_viridis_c() +
  geom_path(data = Mu_slingCurves_df, aes(x = UMAP_1, y = UMAP_2),size = 1,linetype = "dashed")  +
  labs(fill = "Pseudotime",title="Mutant")+
  ylim(-6,5.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=14))
dev.off()


# Figure 3E top - proportion of each cell type along pseudotime
WT_df <- df[which(df$sample_type=="WT"),]
Mu_df <- df[which(df$sample_type=="Mutant"),]
newpalette <- rev(c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF"))
WT_df$cell_type <- factor(WT_df$cell_type,levels=rev(c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes")))
Mu_df$cell_type <- factor(Mu_df$cell_type,levels=rev(c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes")))
pdf(str_c(out_dir,"each_cell_type_proportion_along_pseudotime.pdf"),width=15,height=3)
ggplot(data=WT_df,aes(x=Lineage1,fill=cell_type))+
  geom_density(position = "fill",bw=1) +
  scale_fill_manual(values=newpalette) +
  labs(x="Pseudotime",y="Density",fill="Cell type")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=14))
ggplot(data=Mu_df,aes(x=Lineage1,fill=cell_type))+
  geom_density(position = "fill",bw=1) +
  scale_fill_manual(values=newpalette) +
  labs(x="Pseudotime",y="Density",fill="Cell type")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text = element_text(size=13),axis.title=element_text(size=14))
dev.off()

## Differential Differentiation
#  Evaluate the optimal number of knots required for fitGAM
pdf(str_c(out_dir,"exclude_spermatids_germline_evaluateK_plot.pdf",width=9))
set.seed(100)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 4
icMat <- evaluateK(counts = as.matrix(assays(exclude_spermatids_scRNA.sce)$counts),
                   pseudotime = slingPseudotime(sds, na = FALSE),
                   cellWeights = sds@assays@data$weights,
                   conditions = factor(colData(exclude_spermatids_scRNA.sce)$sample_type),
                   nGenes = 500,
                   k = 3:10,
                   parallel = TRUE,
                   BPPARAM = BPPARAM)
dev.off()


# Fit GAM using 8 knots
set.seed(100)
exclude_spermatids_scRNA.sce@int_metadata$slingshot <- sds
exclude_spermatids_scRNA.sce <- fitGAM(counts = exclude_spermatids_scRNA.sce,
               conditions = factor(colData(exclude_spermatids_scRNA.sce)$sample_type),
               parallel = TRUE,
               BPPARAM = BPPARAM,
               nknots = 8)


## Differential expression between conditions
condRes <- conditionTest(exclude_spermatids_scRNA.sce, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
conditionGenes <- rownames(condRes)[condRes$padj < 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]


## Figure 3E below - Heatmaps of dynamic DEGs between conditions
yhatSmooth <- predictSmooth(exclude_spermatids_scRNA.sce, gene = conditionGenes, nPoints = 100, tidy = FALSE) %>% log1p()
yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))
library(ComplexHeatmap)
yhatSmoothScaled <- yhatSmoothScaled[,c(101:200,1:100)]
pdf(str_c(,out_dir"exclude_spermatids_germline_cells_trajectory_DEGs_heatmap.pdf"),width=8)
Heatmap(yhatSmoothScaled,name="Expression",col=rev(brewer.pal(11,"RdBu")[2:10]),show_row_names=FALSE,cluster_columns=FALSE,cluster_rows=TRUE,show_column_names=FALSE,heatmap_legend_param = list(legend_width = unit(4, "cm"), title_position = "topcenter"),column_split=factor(rep(c("WT","Mutant"),each=100),levels=c("WT","Mutant")),column_gap=unit(3, "mm"))
dev.off()
