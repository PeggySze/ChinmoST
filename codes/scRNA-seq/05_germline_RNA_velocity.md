## Generate spliced and unspliced matrices by Velocyto
- script : germline_velocyto.sh

## Analyze RNA velocity of germline cells
```R
## Load required packages
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(stringr)

## Creat output directory 
out_dir <- "~/ChinmoST/output/scRNA/Velocyto/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load RNA Velocity data from a loom file
dirs <- list.files(out_dir,include.dirs=TRUE)
files <- sapply(1:length(dirs),function(i){
  list.files(str_c(out_dir,dirs[i]),full=TRUE)
  })
files <- files[c(4:6,1:3)]
samples <- c("WT35","WT68","WT911","Mu35","Mu68","Mu911")
ldat_obj.ls <- lapply(1:length(files),function(i){
  ldat <- ReadVelocity(file = files[i])
  obj <- as.Seurat(x = ldat)
  obj
  })

## change gene name and cell name 
flybase_id_to_symbol <- read.table("~/DB/dm6/annotation/dmel-all-r6.33-id2symbol.txt",row.names=1,header=FALSE)
for ( i in 1:length(ldat_obj.ls)){
  ldat.id2symbol <- data.frame(flybase_id=rownames(ldat_obj.ls[[i]]),row.names=rownames(ldat_obj.ls[[i]]))
  ldat.id2symbol <- merge(ldat.id2symbol,flybase_id_to_symbol,by="row.names",sort=FALSE)
  # change gene name
  for ( assay in c("spliced","unspliced","ambiguous")){
    convertID.assay <- ldat_obj.ls[[i]][[assay]]@counts
    rownames(convertID.assay) <- ldat.id2symbol$V2
    convertID_assay <- CreateAssayObject(counts=convertID.assay)
    ldat_obj.ls[[i]][[assay]] <- convertID_assay
  }
  # change cell name
  cell_names <-  sapply(1:ncol(ldat_obj.ls[[i]]),function(j){unlist(strsplit(colnames(ldat_obj.ls[[i]])[j],":"))[2]})
  cell_names <- gsub("x","-1",cell_names)  
  ldat_obj.ls[[i]] <- RenameCells(object=ldat_obj.ls[[i]],new.names=cell_names)
}

WT_merged_obj <- merge(x = ldat_obj.ls[[1]], 
  y = c(ldat_obj.ls[[2]],ldat_obj.ls[[3]]), 
  add.cell.ids=samples[1:3])

Mutant_merged_obj <-  merge(x = ldat_obj.ls[[4]], 
  y = c(ldat_obj.ls[[5]],ldat_obj.ls[[6]]), 
  add.cell.ids=samples[4:6])

bm <- merge(x = ldat_obj.ls[[1]], 
  y = c(ldat_obj.ls[[2]],ldat_obj.ls[[3]],ldat_obj.ls[[4]],ldat_obj.ls[[5]],ldat_obj.ls[[6]]), 
  add.cell.ids=samples)

WT_spliced <- CreateAssayObject(GetAssayData(WT_merged_obj, assay = "spliced"))
WT_unspliced <- CreateAssayObject(GetAssayData(WT_merged_obj, assay = "unspliced"))
WT_ambiguous <- CreateAssayObject(GetAssayData(WT_merged_obj, assay = "ambiguous"))
Mutant_spliced <- CreateAssayObject(GetAssayData(Mutant_merged_obj, assay = "spliced"))
Mutant_unspliced <- CreateAssayObject(GetAssayData(Mutant_merged_obj, assay = "unspliced"))
Mutant_ambiguous <- CreateAssayObject(GetAssayData(Mutant_merged_obj, assay = "ambiguous"))

## Load germline Seurat object
germline_cells.scRNA.integrated <- readRDS("~/ChinmoST/output/scRNA/seurat/subclustering/germline_subclustering/integrated_germline_cells_scRNA.rds")


## Exclude spermatids in the following analysis
Idents(germline_cells.scRNA.integrated) <- germline_cells.scRNA.integrated$cell_type
WT_early_germline_cells.scRNA.integrated <- subset(germline_cells.scRNA.integrated,subset=sample_type=="WT" & cell_type %in% c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes") & UMAP_1<=2.5 & tree.ident!=2)
Mutant_early_germline_cells.scRNA.integrated <- subset(germline_cells.scRNA.integrated,subset=sample_type=="Mutant" & cell_type %in% c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes")& UMAP_1<=2.5 & tree.ident!=2)

## filter and order matrix
genes <- rownames(germline_cells.scRNA.integrated)
WT_cells <- colnames(WT_early_germline_cells.scRNA.integrated)
WT_spliced <- WT_spliced[genes,WT_cells]
WT_unspliced <- WT_unspliced[genes,WT_cells]
WT_ambiguous <- WT_ambiguous[genes,WT_cells]
WT_early_germline_cells.scRNA.integrated[["spliced"]] <- CreateAssayObject(counts=WT_spliced)
WT_early_germline_cells.scRNA.integrated[["unspliced"]] <- CreateAssayObject(counts=WT_unspliced)
WT_early_germline_cells.scRNA.integrated[["ambiguous"]] <- CreateAssayObject(counts=WT_ambiguous)
Mutant_cells <- colnames(Mutant_early_germline_cells.scRNA.integrated)
Mutant_spliced <- Mutant_spliced[genes,Mutant_cells]
Mutant_unspliced <- Mutant_unspliced[genes,Mutant_cells]
Mutant_ambiguous <- Mutant_ambiguous[genes,Mutant_cells]
Mutant_early_germline_cells.scRNA.integrated[["spliced"]] <- CreateAssayObject(counts=Mutant_spliced)
Mutant_early_germline_cells.scRNA.integrated[["unspliced"]] <- CreateAssayObject(counts=Mutant_unspliced)
Mutant_early_germline_cells.scRNA.integrated[["ambiguous"]] <- CreateAssayObject(counts=Mutant_ambiguous)


## Run RNA Velocty
WT_early_germline_cells.scRNA.integrated <- RunVelocity(object = WT_early_germline_cells.scRNA.integrated, deltaT = 1, kCells = 10, fit.quantile = 0.02)
Mutant_early_germline_cells.scRNA.integrated <- RunVelocity(object = Mutant_early_germline_cells.scRNA.integrated, deltaT = 1, kCells = 10, fit.quantile = 0.02)


## Figure 3D - Visualize RNA velocities in UMAP
ident.colors <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF")
names(ident.colors) <- c("GSCs/spermatogonia","Early spermatocytes","Intermediate spermatocytes","Late spermatocytes")
WT_cell.colors <- ident.colors[Idents(object = WT_early_germline_cells.scRNA.integrated)]
names(WT_cell.colors) <- colnames(WT_early_germline_cells.scRNA.integrated)
Mutant_cell.colors <- ident.colors[Idents(object = Mutant_early_germline_cells.scRNA.integrated)]
names(Mutant_cell.colors) <- colnames(Mutant_early_germline_cells.scRNA.integrated)
pdf(str_c(out_dir,"WT_Mutant_early_germline_cells_RNA_velocities.pdf"),width=6.5)
show.velocity.on.embedding.cor(emb = Embeddings(object = WT_early_germline_cells.scRNA.integrated, reduction = "umap"),vel=Tool(object = WT_early_germline_cells.scRNA.integrated, slot = "RunVelocity"),cell.colors=ac(x = WT_cell.colors, alpha = 0.5),cex = 0.8, arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
show.velocity.on.embedding.cor(emb = Embeddings(object = Mutant_early_germline_cells.scRNA.integrated, reduction = "umap"),vel=Tool(object = Mutant_early_germline_cells.scRNA.integrated, slot = "RunVelocity"),cell.colors=ac(x = Mutant_cell.colors, alpha = 0.5),cex = 0.8, arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
dev.off()
```