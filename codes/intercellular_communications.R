## Load required packages
library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(stringr) 
library(RColorBrewer)
library(xlsx)
library(ComplexHeatmap)
library(reshape2)
options(stringsAsFactors = FALSE)

## Creat output directory
out_dir <- "~/ChinmoST/output/scRNA/CellChat/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

# Load Seurat object 
testis_scRNA.integrated <- readRDS("~/ChinmoST/output/scRNA/seurat/integrated_analysis/integrated_PC30_testis_scRNA.rds")

## Remove uninterested cell types
subset_testis_scRNA <- subset(testis_scRNA.integrated,subset=cell_type %in% c("Hub cells","CySCs","Early cyst cells","Intermediate cyst cells","Late cyst cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids"))

## Prepare Cellchat input data
WT.data.input <- subset_testis_scRNA[["RNA"]]@data[,rownames(subset_testis_scRNA@meta.data)[which(subset_testis_scRNA$sample_type=="WT")]]
WT.meta <- subset_testis_scRNA@meta.data[which(subset_testis_scRNA$sample_type=="WT"),"cell_type",drop=FALSE]
Mu.data.input <- subset_testis_scRNA[["RNA"]]@data[,rownames(subset_testis_scRNA@meta.data)[which(subset_testis_scRNA$sample_type=="Mutant")]]
Mu.meta <- subset_testis_scRNA@meta.data[which(subset_testis_scRNA$sample_type=="Mutant"),"cell_type",drop=FALSE]


## Create CellChat objects
WT.cellchat <- createCellChat(object = WT.data.input, meta = WT.meta, group.by = "cell_type")
Mu.cellchat <- createCellChat(object = Mu.data.input, meta = Mu.meta, group.by = "cell_type")

## set "cell_type" as default cell identity
WT.cellchat <- setIdent(WT.cellchat, ident.use = "cell_type") 
Mu.cellchat <- setIdent(Mu.cellchat, ident.use = "cell_type") 
WT.cellchat@idents <- factor(WT.cellchat@idents,levels=c("Late spermatids","Late cyst cells","Intermediate cyst cells","Early cyst cells","CySCs","Hub cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids"))
Mu.cellchat@idents <- factor(Mu.cellchat@idents,levels=c("Late spermatids","Late cyst cells","Intermediate cyst cells","Early cyst cells","CySCs","Hub cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids"))
WT.groupSize <- as.numeric(table(WT.cellchat@idents)) # number of cells in each cell group
Mu.groupSize <- as.numeric(table(Mu.cellchat@idents))


## Set the ligand-receptor interaction database for drosophila
interaction_input <- read.csv(file ='~/ChinmoST/input/CellChat/drosophila_interaction_input_CellChatDB.csv',row.names = 1)
complex_input <- read.csv(file = '~/ChinmoST/input/CellChat/drosophila_complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = '~/ChinmoST/input/CellChat/drosophila_cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = '~/ChinmoST/input/CellChat/drosophila_geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
CellChatDB.use <- CellChatDB
WT.cellchat@DB <- CellChatDB.use
Mu.cellchat@DB <- CellChatDB.use


## Figure S5A - percentage of each interaction type
df <- interaction_input %>%
  group_by(annotation) %>%
  summarise(counts=n()) %>%
  mutate(percentage=counts/sum(counts)*100)
df$annotation <- factor(df$annotation,levels=c("Secreted Signaling","ECM-receptor","Cell-Cell Contact"),labels=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))
newpalette <- brewer.pal(8,"Set2")[1:3]
pdf(str_c(out_dir,"database_interaction_type_percentage_piechart.pdf"))
ggplot(data=df,aes(x='Feature',y=percentage,fill=annotation,color="black")) +
  geom_bar(stat = 'identity', position = 'stack',color="black",size=0.02)  + 
  coord_polar(theta = 'y') + 
  labs(x = '', y = '', title = '') +
  scale_fill_manual(values = newpalette)+
  theme_bw() + 
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid=element_blank())
dev.off()


## Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
WT.cellchat <- subsetData(WT.cellchat)
Mu.cellchat <- subsetData(Mu.cellchat)

# identify over-expressed ligands or receptors in one cell group
WT.cellchat <- identifyOverExpressedGenes(WT.cellchat)
Mu.cellchat <- identifyOverExpressedGenes(Mu.cellchat)

# identify over-expressed ligand-receptor interactions 
WT.cellchat <- identifyOverExpressedInteractions(WT.cellchat)
Mu.cellchat <- identifyOverExpressedInteractions(Mu.cellchat)


## Inference of cell-cell communication network
# CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test.
# By default, CellChat uses a statistically robust mean method called ‘trimean’ to calculate the average gene expression per cell group. we specify type="truncatedMean",trim =0.1 to calculate the average gene expression by discarding 10% from each end of the data
WT.cellchat <- computeCommunProb(WT.cellchat, raw.use = TRUE,type="truncatedMean",trim=0.1) 
Mu.cellchat <- computeCommunProb(Mu.cellchat, raw.use = TRUE,type="truncatedMean",trim=0.1)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
WT.cellchat <- filterCommunication(WT.cellchat, min.cells = 10)
Mu.cellchat <- filterCommunication(Mu.cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame
WT.df.net <- subsetCommunication(WT.cellchat)
Mu.df.net <- subsetCommunication(Mu.cellchat)

# Supplemental Table S3
write.xlsx(WT.df.net,str_c(out_dir,"truncatedMean0.1_WT_cell_communication_network_df.xlsx"),row.names=FALSE)
write.xlsx(Mu.df.net,str_c(out_dir,"truncatedMean0.1_Mutant_cell_communication_network_df.xlsx"),row.names=FALSE)

## Infer the cell-cell communication at a signaling pathway level
# CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway
# Note : The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
WT.cellchat <- computeCommunProbPathway(WT.cellchat)
Mu.cellchat <- computeCommunProbPathway(Mu.cellchat)


## Calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. 
WT.cellchat <- aggregateNet(WT.cellchat)
Mu.cellchat <- aggregateNet(Mu.cellchat)


## Figure S5C - visualize BMP and Hedgehog signaling between WT hub cells & CySCs & GSCs
df <- WT.df.net %>%filter(pathway_name %in% c("BMP","Hedgehog"),source %in% c("Hub cells","CySCs"),target %in% c("CySCs","GSCs"))
df <- df[-which(df$source=="CySCs"&df$target=="CySCs"),]
df <- df[-which(df$interaction_name %in% c("dpp_dally","dpp_tkv")),]
df$source2target <- str_c(df$source," -> ",df$target)
df$source2target <- factor(df$source2target,levels=c("Hub cells -> CySCs","Hub cells -> GSCs","CySCs -> GSCs"))
df$pathway_name <- factor(df$pathway_name,levels=c("BMP","Hedgehog"))
df <- df[,c("interaction_name_2","source2target","prob")]
cast_df <- dcast(df,interaction_name_2~source2target,value.var="prob")
rownames(cast_df) <- cast_df$interaction_name_2
cast_df <- cast_df[,-1]
pdf(str_c(out_dir,"WT_niche_BMP_Hedgehog_signaling.pdf"))
Heatmap(cast_df,name = "Communication probability",na_col = "white",cluster_rows=FALSE,cluster_columns=FALSE,row_names_side = "left",row_names_rot = 0,column_names_rot = 45,rect_gp = gpar(col = "grey", lwd = 0.5),col=colorRamp2(seq(-max(abs(cast_df),na.rm =TRUE), max(abs(cast_df),na.rm =TRUE), length = 3), c("#2488F0","white","#E22929")),row_title="ligand - (receptor)",row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),row_title_gp = gpar(fontsize = 15),column_title_gp = gpar(fontsize = 16),border_gp = gpar(col = "black", lwd = 2))
dev.off()


## Figure S5B - visualize cell-cell communication network
# WT
WT.groupSize <- as.numeric(table(WT.cellchat@idents))
newpalette <- c(colorRampPalette(brewer.pal(8,"Set2"))(16)[1:5],colorRampPalette(brewer.pal(8,"Set2"))(16)[c(9,11,8,13,15,14,16)])[c(12,5:1,7:11)]
pdf(str_c(out_dir,"WT_cell_communication_network_circle_plot.pdf"))
netVisual_circle(WT.cellchat@net$count,vertex.weight = WT.groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions in Control testes",vertex.label.cex=1,color.use=newpalette,arrow.size=0.3)
netVisual_circle(WT.cellchat@net$weight, vertex.weight = WT.groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength in Control testes",vertex.label.cex=1,color.use=newpalette,arrow.size=0.3)
dev.off()


# mutant
Mu.groupSize <- as.numeric(table(Mu.cellchat@idents))
newpalette <- c(colorRampPalette(brewer.pal(8,"Set2"))(16)[1:5],colorRampPalette(brewer.pal(8,"Set2"))(16)[c(9,11,8,13,15,14,16)])[c(12,5:1,7:11)]
pdf(str_c(out_dir,"Mutant_cell_communication_network_circle_plot.pdf"))
netVisual_circle(Mu.cellchat@net$count, vertex.weight = Mu.groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions in Mutant testes",vertex.label.cex=1,color.use=newpalette,arrow.size=0.3)
netVisual_circle(Mu.cellchat@net$weight, vertex.weight = Mu.groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength in Mutant testes",vertex.label.cex=1,color.use=newpalette,arrow.size=0.3)
dev.off()

### Mutant vs WT
## merge different CellChat objects together
object.list <- list(WT = WT.cellchat, Mutant = Mu.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

## Figure S5D - Compare the major sources and targets in 2D space to identify the cell populations with significant changes in sending or receiving signals between different datasets.
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
newpalette <- c(colorRampPalette(brewer.pal(8,"Set2"))(16)[1:5],colorRampPalette(brewer.pal(8,"Set2"))(16)[c(9,11,8,13,15,14,16)])[c(12,5:1,7:11)] 
pdf(str_c(out_dir,"WT_and_Mu_netAnalysis_signalingRole_scatter.pdf"),width=10.5,height=6)
gg <- list()
for (i in 1:length(object.list)) {
  # Compute the network centrality scores
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use=newpalette)
}
my_plot <- list()
WT_signalingRole_df <- gg[[1]]$data
Mu_signalingRole_df <- gg[[2]]$data
WT_signalingRole_df$labels <- factor(WT_signalingRole_df$labels,levels=c("Late spermatids","Late cyst cells","Intermediate cyst cells","Early cyst cells","CySCs","Hub cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids"))
Mu_signalingRole_df$labels <- factor(Mu_signalingRole_df$labels,levels=c("Late spermatids","Late cyst cells","Intermediate cyst cells","Early cyst cells","CySCs","Hub cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids"))
my_plot[[1]] <- ggplot(data = WT_signalingRole_df, aes(x, y)) +
  geom_point(aes(colour = labels, fill = labels),size=3) +
  theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"),plot.title = element_text(size = 12, face = "plain",hjust = 0.5)) + 
  labs(title="Control testes",x = "Outgoing interaction strength", y = "Incoming interaction strength",axis.line.x =element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
  scale_fill_manual(values=newpalette) +
  scale_colour_manual(values=newpalette) +
  xlim(0,2.4)+
  ylim(0,2.9)+
  #scale_x_continuous(breaks = seq(0,5,1)) +
  #scale_y_continuous(breaks = seq(0,3,0.5)) +
  ggrepel::geom_text_repel(mapping = aes(label = labels), size = 3, show.legend = F, 
            segment.size = 0.2, segment.alpha = 0.5) +
  CellChat_theme_opts() +
  theme(legend.position = "none") 
my_plot[[2]] <- ggplot(data = Mu_signalingRole_df, aes(x, y)) +
  geom_point(aes(colour = labels, fill = labels),size=3) +
  theme(text = element_text(size = 10), legend.key.height = grid::unit(0.15, "in"),plot.title = element_text(size = 12, face = "plain",hjust = 0.5)) + 
  labs(title="Mutant testes",x = "Outgoing interaction strength", y = "Incoming interaction strength",axis.line.x =element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) +
  scale_fill_manual(values=newpalette) +
  scale_colour_manual(values=newpalette) +
  xlim(0,2.4)+
  ylim(0,2.9)+
  #scale_x_continuous(breaks = seq(0,5,1)) +
  #scale_y_continuous(breaks = seq(0,3,0.5)) +
  ggrepel::geom_text_repel(mapping = aes(label = labels), size = 3, show.legend = F, 
            segment.size = 0.2, segment.alpha = 0.5) +
  CellChat_theme_opts() +
  theme(legend.position = "bottom") 
patchwork::wrap_plots(plots = my_plot)
dev.off()


## Compare the number of interactions and interaction strength among different cell populations
# Figure 4A - visualize differential number of interactions or interaction strength by heatmap
cell_types <- c("Hub cells","CySCs","Early cyst cells","Intermediate cyst cells","Late cyst cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids")
WT_cell_type_communication_weight <- cellchat@net$WT$weight
WT_cell_type_communication_weight <- WT_cell_type_communication_weight[cell_types,cell_types]
Mu_cell_type_communication_weight <- cellchat@net$Mutant$weight
Mu_cell_type_communication_weight <- Mu_cell_type_communication_weight[cell_types,cell_types]
Mu_vs_WT_communication_weight <- Mu_cell_type_communication_weight-WT_cell_type_communication_weight
Mu_vs_WT_communication_weight[Mu_vs_WT_communication_weight == 0] <- NA
pdf(str_c(out_dir,"Mu_vs_WT_interaction_weight_heatmap.pdf"),height=6,width=7.8)
Heatmap(Mu_vs_WT_communication_weight,name = "Relative strength",na_col = "white",cluster_rows=FALSE,cluster_columns=FALSE,row_names_side = "left",row_names_rot = 0,column_title = "Mutant vs Control differential interaction strength",column_names_rot = 45,rect_gp = gpar(col = "grey", lwd = 0.5),col=colorRamp2(seq(-max(abs(Mu_vs_WT_communication_weight),na.rm =TRUE), max(abs(Mu_vs_WT_communication_weight),na.rm =TRUE), length = 3), c("#2488F0","white","#E22929")),row_title="Sources (Senders)",row_names_gp = gpar(fontsize = 14),column_names_gp = gpar(fontsize = 14),row_title_gp = gpar(fontsize = 15),column_title_gp = gpar(fontsize = 16),border_gp=gpar(col = "black", lwd = 2))
dev.off()


## Figure 4B - visualize differtial interaction strength of niche cells and other cell types  by Circle plot
edited_netVisual_diffInteraction <- function(object, comparison = c(1,2), measure = c("count", "weight", "count.merged", "weight.merged"), color.use = NULL, color.edge = c('#b2182b','#2166ac'), title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                      weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex=1,vertex.label.color= "black",
                                      edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                      edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2,
                                      arrow.width=1,arrow.size = 0.2){
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }

  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  net[abs(net) < stats::quantile(abs(net), probs = 1-top,na.rm = TRUE)] <- 0
  net[is.na(net)] <- 0

  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)

  igraph::E(g)$weight <- abs(igraph::E(g)$weight)

  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }


  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}

newpalette <- c(colorRampPalette(brewer.pal(8,"Set2"))(16)[1:5],colorRampPalette(brewer.pal(8,"Set2"))(16)[c(9,11,8,13,15,14,16)])[c(12,5:1,7:11)]
pdf(str_c(out_dir,"niche_sources_WT_vs_Mu_interactions_circle_plot.pdf"))
edited_netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",vertex.label.cex=0.8,title.name="Mutant vs WT differential interaction strength from hub cells",color.use=newpalette,edge.width.max=12,arrow.size=0.5,margin=0.1,sources.use=c("Hub cells"))
edited_netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",vertex.label.cex=0.8,title.name="Mutant vs WT differential interaction strength from CySCs",color.use=newpalette,edge.width.max=12,arrow.size=0.5,margin=0.1,sources.use=c("CySCs"))
edited_netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",vertex.label.cex=0.8,title.name="Mutant vs WT differential interaction strength from GSCs",color.use=newpalette,edge.width.max=12,arrow.size=0.5,margin=0.1,sources.use=c("GSCs"))
dev.off()

pdf(str_c(out_dir,"niche_targets_WT_vs_Mu_interactions_circle_plot.pdf"))
edited_netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",vertex.label.cex=0.8,title.name="Mutant vs WT differential interaction strength to hub cells",color.use=newpalette,edge.width.max=12,arrow.size=0.5,margin=0.1,targets.use =c("Hub cells"))
edited_netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",vertex.label.cex=0.8,title.name="Mutant vs WT differential interaction strength to CySCs",color.use=newpalette,edge.width.max=12,arrow.size=0.5,margin=0.1,targets.use =c("CySCs"))
edited_netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",vertex.label.cex=0.8,title.name="Mutant vs WT differential interaction strength to GSCs",color.use=newpalette,edge.width.max=12,arrow.size=0.5,margin=0.1,targets.use =c("GSCs"))
dev.off()


## Figure 4C - Compare outgoing (or incoming) signaling associated with each cell population
cell_types <- c("Hub cells","CySCs","Early cyst cells","Intermediate cyst cells","Late cyst cells","GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes","Early spermatids","Late spermatids")
signaling <- c("BMP","Activin","JAK-STAT","EGFR","TNF","Hedgehog","FGFR","Insulin","Wnt-TCF ","Wnt/Ca2+","Notch ","Pvr","Hippo","Toll ")
WT_centr <- slot(object.list[[1]], "netP")$centr
WT_outgoing <- matrix(0, nrow = nlevels(object.list[[1]]@idents), ncol = length(WT_centr))
WT_incoming <- matrix(0, nrow = nlevels(object.list[[1]]@idents), ncol = length(WT_centr))
dimnames(WT_outgoing) <- list(levels(object.list[[1]]@idents), names(WT_centr))
dimnames(WT_incoming) <- dimnames(WT_outgoing)
for (i in 1:length(WT_centr)) {
    WT_outgoing[,i] <- WT_centr[[i]]$outdeg
    WT_incoming[,i] <- WT_centr[[i]]$indeg
  }
WT_outgoing_mat <- t(WT_outgoing)
WT_incoming_mat <- t(WT_incoming)
WT_outgoing_mat <- WT_outgoing_mat[signaling,cell_types]
WT_incoming_mat <- WT_incoming_mat[signaling,cell_types]

Mu_centr <- slot(object.list[[2]], "netP")$centr
Mu_outgoing <- matrix(0, nrow = nlevels(object.list[[2]]@idents), ncol = length(Mu_centr))
Mu_incoming <- matrix(0, nrow = nlevels(object.list[[2]]@idents), ncol = length(Mu_centr))
dimnames(Mu_outgoing) <- list(levels(object.list[[2]]@idents), names(Mu_centr))
dimnames(Mu_incoming) <- dimnames(Mu_outgoing)
for (i in 1:length(Mu_centr)) {
    Mu_outgoing[,i] <- Mu_centr[[i]]$outdeg
    Mu_incoming[,i] <- Mu_centr[[i]]$indeg
  }
Mu_outgoing_mat <- t(Mu_outgoing)
Mu_incoming_mat <- t(Mu_incoming)
Mu_outgoing_mat <- Mu_outgoing_mat[signaling,cell_types]
Mu_incoming_mat <- Mu_incoming_mat[signaling,cell_types]

Mu_vs_WT_outgoing_mat <- Mu_outgoing_mat-WT_outgoing_mat
Mu_vs_WT_incoming_mat <- Mu_incoming_mat-WT_incoming_mat
Mu_vs_WT_outgoing_mat[Mu_vs_WT_outgoing_mat == 0] <- NA
Mu_vs_WT_incoming_mat[Mu_vs_WT_incoming_mat == 0] <- NA


melt_Mu_vs_WT_outgoing_mat <-  melt(Mu_vs_WT_outgoing_mat)
melt_Mu_vs_WT_outgoing_mat$Var2 <- factor(melt_Mu_vs_WT_outgoing_mat$Var2,levels=cell_types)
melt_Mu_vs_WT_outgoing_mat$Var1 <- factor(melt_Mu_vs_WT_outgoing_mat$Var1,levels=rev(signaling))
melt_Mu_vs_WT_incoming_mat <- melt(Mu_vs_WT_incoming_mat)
melt_Mu_vs_WT_incoming_mat$Var2 <- factor(melt_Mu_vs_WT_incoming_mat$Var2,levels=cell_types)
melt_Mu_vs_WT_incoming_mat$Var1 <- factor(melt_Mu_vs_WT_incoming_mat$Var1,levels=rev(signaling))

pdf(str_c(out_dir,"Mu_vs_WT_signalingRole_heatmap.pdf"),height=6)
ggplot(data=melt_Mu_vs_WT_outgoing_mat,aes(x=Var2,y=Var1,fill=value))+
  geom_tile(colour="grey")+
  scale_fill_gradientn(colours = brewer.pal(11, name = "RdBu")[c(10,6,2)],limits=c(-0.05,0.05),oob = scales::squish,na.value = 'white')+
  theme_bw()+
  labs(x="Cell type",y="Signaling",fill="Relative strength",title="Mutant vs WT outgoing signaling")+
   theme(axis.text.x = element_text(angle = 45,color='black',hjust = 1,size=10),axis.text.y=element_text(color='black',size=11),plot.title = element_text(hjust = 0.5,size=13),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggplot(data=melt_Mu_vs_WT_incoming_mat,aes(x=Var2,y=Var1,fill=value))+
  geom_tile(colour="grey")+
  scale_fill_gradientn(colours = brewer.pal(11, name = "RdBu")[c(10,6,2)],limits=c(-0.05,0.05),oob = scales::squish,na.value = 'white')+
  theme_bw()+
  labs(x="Cell type",y="Signaling",fill="Relative strength",title="Mutant vs WT incoming signaling")+
   theme(axis.text.x = element_text(angle = 45,color='black',hjust = 1,size=10),axis.text.y=element_text(color='black',size=11),plot.title = element_text(hjust = 0.5,size=13),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()


## Figure 4D - Visualze relative contribution of each ligand-receptor pair to the overall communication network of insulin signaling pathway target to mutant GSCs/spermatogonia
pairLR <- searchPair(signaling = "Insulin", pairLR.use = Mu.cellchat@LR$LRsig,key = "pathway_name", matching.exact = T, pair.only = T)
net <- Mu.cellchat@net
pairLR.use.name <- dimnames(net$prob)[[3]]
pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
pairLR <- pairLR[pairLR.name, ]
prob <- net$prob
pval <- net$pval
prob[pval > 0.05] <- 0
if (length(pairLR.name) > 1) {
        pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
            3, sum) != 0]
    }
if (length(pairLR.name.use) != 0) {
        pairLR <- pairLR[pairLR.name.use, ]
    }
nRow <- length(pairLR.name.use)
prob <- prob[, , pairLR.name.use]
pval <- pval[, , pairLR.name.use]
Ilp6_InR_df <- as.data.frame(prob[,"GSCs/spermatogonia","Ilp6_InR",drop=FALSE])
colnames(Ilp6_InR_df) <- "GSCs/spermatogonia"
Ilp6_InR_df$ligand_receptor <- "Ilp6-InR"
Ilp8_InR_df <- as.data.frame(prob[,"GSCs/spermatogonia","Ilp8_InR",drop=FALSE])
colnames(Ilp8_InR_df) <- "GSCs/spermatogonia"
Ilp8_InR_df$ligand_receptor <- "Ilp8-InR"
df <- rbind(Ilp6_InR_df,Ilp8_InR_df)
df <- df %>% 
  group_by(ligand_receptor) %>%
  summarise(sum=sum(GSCs)) %>%
  mutate(ratio=sum/sum(sum))
pdf(str_c(out_dir,"insulin_signaling_ligand_receptor_relative_contribution.pdf"),width=4,height=5)
ggplot(data=df,aes(x=ligand_receptor,y=ratio)) +
  geom_bar(stat="identity",width=0.5)+
  labs(x="",y="Relative contribution",title="Contribution of each L-R pair in insulin siganling")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),plot.title = element_text(hjust = 0.5,size=16),axis.text.y = element_text(size=13),axis.title=element_text(size=14),axis.text.x=element_text(size=14))
dev.off()


## Figure 4E - Chord diagram visualizing cell-cell communication mediated by Ilp6-InR
newpalette <- c(colorRampPalette(brewer.pal(8,"Set2"))(16)[1:7],brewer.pal(8,"Accent")[3],colorRampPalette(brewer.pal(8,"Set2"))(16)[c(12,9:11,8,13,15,14,16)])[c(4,3,2,1,12)]
names(newpalette) <- c("Intermediate cyst cells","Early cyst cells","CySCs","Hub cells","GSCs/spermatogonia")
pairLR.insulin <- extractEnrichedLR(WT.cellchat, signaling = "Insulin", geneLR.return = FALSE)
pairLR.insulin <- pairLR.insulin[1,,drop=FALSE]
pdf(str_c(out_dir,"insulin_signaling_chord_diagram.pdf"))
netVisual_chord_gene(WT.cellchat,pairLR.use=pairLR.insulin,sources.use=c(3:7),targets.use=c(5,6,7),legend.pos.x=5,legend.pos.y=8,title.name="Ilp6-InR signaling in control testis",lab.cex=1,color.use=newpalette)
netVisual_chord_gene(Mu.cellchat,pairLR.use=pairLR.insulin,sources.use=c(3:7),targets.use=c(5,6,7),legend.pos.x=5,legend.pos.y=8,title.name="Ilp6-InR signaling in mutant testis",lab.cex=1,color.use=newpalette)
dev.off()


## Figure 4H,S5E - Violin plot showing expression of interested genes involved in insulin signaling pathway
DefaultAssay(testis_scRNA.integrated) <- "RNA"
niche.scRNA <- subset(testis_scRNA.integrated,subset=cell_type %in% c("Hub cells","CySCs","GSCs/spermatogonia"))
niche.scRNA$sample_type <- ifelse(niche.scRNA$orig.ident %in% c("CT35","CT68","CT911"),"WT","Mutant")
genes.ls <- list("genes_01"=c("ImpL2","Ilp6","InR"),"genes_02"=c("Akt1","Pi3K92E","foxo"))
newpalette <- c(brewer.pal(9, "Blues")[4],brewer.pal(9, "YlOrRd")[4])
pdf(str_c(out_dir,"interested_genes_expression_in_hub_CySCs_GSCs_VlnPlot.pdf",width=8,height=7))
for (i in 1:length(genes.ls)){
	df <- FetchData(object=niche.scRNA,c(genes.ls[[i]],"cell_type","sample_type"),slot="data")
	df <- reshape2::melt(df,id=c("cell_type","sample_type"),variable.name="gene",value.name ="expression")
	df$sample_type <- factor(df$sample_type,levels=c("WT","Mutant"))
	df$cell_type <- factor(df$cell_type,levels=c("Hub cells","CySCs","GSCs/spermatogonia"))
	df$gene <- factor(df$gene,levels=genes.ls[[i]])
	p <- ggplot(df,aes(x=cell_type,y=expression,fill=sample_type))+
  		geom_violin(position=position_dodge(0.9),scale="width",trim=TRUE)+
  		geom_boxplot(width=0.1,position=position_dodge(0.9),outlier.size=0.8) +
  		scale_fill_manual(values=newpalette)+
  		ylim(0,3.1)+
  		facet_grid(gene ~ .,switch = "y")+
  		labs(x="Sample",y="Expression Level",fill="Sample type") +
  		theme_bw() +
  		theme(plot.title = element_text(face = "bold",size = rel(1.5), hjust = 0.5),
  			panel.background=element_rect(fill='transparent', color='black',linetype="solid"),
      		axis.title = element_text(face = "bold",size = rel(1.2)),
      		axis.title.y = element_text(angle=90,vjust =2),
      		axis.title.x = element_text(vjust = -0.2),
      		axis.text = element_text(color="black",size=rel(1.2)), 
      		axis.line = element_line(colour="black"),
      		axis.ticks = element_line(),
      		panel.grid.major = element_blank(),
      		panel.grid.minor = element_blank(),
      		legend.title = element_text(),
      		plot.margin=unit(c(10,5,5,5),"mm"),
      		strip.background=element_rect(colour="black",fill="#f0f0f0",linetype="solid"),
      		strip.text = element_text(face="bold",size=14),legend.position = "bottom")
  	print(p)
}
dev.off()
