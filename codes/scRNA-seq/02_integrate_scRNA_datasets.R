### Integrate scRNA-seq of WT and ChinmoST testes across three timepoints by Seurat
## Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr) 
library(RColorBrewer)
library(xlsx)

## Creat output directory 
out_dir <- "~/ChinmoST/output/scRNA/seurat/integrated_analysis/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load the dataset
folders <- list.files('~/ChinmoST/output/scRNA/cellranger',full.names = TRUE)
folders <- c(folders[4:6],folders[1:3])
samples <- c("WT35","WT68","WT911","Mu35","Mu68","Mu911")
objList <- lapply(1:length(folders),function(i){
	CreateSeuratObject(counts = Read10X_h5(str_c(folders[i],"/outs/filtered_feature_bc_matrix.h5")),
	project = samples[i],assay = "RNA") 
})

## convert gene names from flybase id to gene symbol
flybase_id_to_symbol <- read.table("~/DB/dm6/annotation/dmel-all-r6.33-id2symbol.txt",row.names=1,header=FALSE)
for ( i in 1:length(objList)){
  scRNA.id2symbol <- data.frame(flybase_id=rownames(objList[[i]]),row.names=rownames(objList[[i]]))
  scRNA.id2symbol <- merge(scRNA.id2symbol,flybase_id_to_symbol,by="row.names",sort=FALSE)
  convertID.counts <- objList[[i]][["RNA"]]@counts
  rownames(convertID.counts) <- scRNA.id2symbol$V2
  renamed_RNA_assay <- CreateAssayObject(counts=convertID.counts)
  objList[[i]][["RNA"]] <- renamed_RNA_assay
}

## QC and filter out low-quality cells
for ( i in 1:length(objList)){
  objList[[i]][["percent.mt"]] = PercentageFeatureSet(objList[[i]],pattern="^mt:")
  pdf(str_c(out_dir,samples[i],"_qc_plot.pdf"),width=9)
  p1=VlnPlot(objList[[i]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2=ggplot(data=objList[[i]][["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
  p3=ggplot(data=objList[[i]][["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
  p4=ggplot(data=objList[[i]][["percent.mt"]],aes(x=percent.mt))+geom_density()
  p5=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p6=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
}

## Retain cells that express more than 200 genes and less than 6000 genes and had mitochondrial content <10%
for (i in 1:length(objList)){
  objList[[i]] <- subset(objList[[i]],subset = nFeature_RNA >=200 & nFeature_RNA <= 6000 &percent.mt < 10 )
}


## Simply merge Seurat objects
merged_obj <- merge(objList[[1]],
                    y=c(objList[[2]],objList[[3]],objList[[4]],objList[[5]],objList[[6]]),
                    add.cell.ids = samples)

## split the combined object into a list, with each dataset as an element
testis_scRNA.list <- SplitObject(merged_obj,split.by = "ident")


## perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(testis_scRNA.list)){
  testis_scRNA.list[[i]] <- NormalizeData(testis_scRNA.list[[i]],verbose = FALSE)
  testis_scRNA.list[[i]] <- FindVariableFeatures(testis_scRNA.list[[i]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

## identify anchors using the FindIntegrationAnchors function
testis_scRNA.anchors <- FindIntegrationAnchors(object.list = testis_scRNA.list,anchor.features = 2000,dims = 1:30)


## pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
testis_scRNA.integrated <- IntegrateData(anchorset = testis_scRNA.anchors, dims = 1:30,features.to.integrate = rownames(testis_scRNA.list[[1]]))


## switch to integrated assay
DefaultAssay(testis_scRNA.integrated) <- "integrated"

## scale and center features in the dataset
testis_scRNA.integrated <- ScaleData(testis_scRNA.integrated, features =rownames(testis_scRNA.integrated))

## Perform linear dimensional reduction
testis_scRNA.integrated <- RunPCA(testis_scRNA.integrated, npcs = 50, verbose = FALSE)

## Determine the ‘dimensionality’ of the dataset
# JackStraw 
testis_scRNA.integrated <- JackStraw(testis_scRNA.integrated, num.replicate = 100, dims =50)
testis_scRNA.integrated  <- ScoreJackStraw(testis_scRNA.integrated, dims = 1:50)
pdf(str_c(out_dir,"integrated_pc.pdf"))
JackStrawPlot(testis_scRNA.integrated, dims = 1:50)
# ‘Elbow plot’
ElbowPlot(testis_scRNA.integrated,ndims=50)
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
  "Nep5", "Nep1", # seminal vesicle & testis epithelia
  "Sox100B", "ems",# testis pigment cell
  "Antp", "Gasp", "Nplp4", "Ubx", "btl", "ect", "sano", "trh", "vvl", "wat", # adult tracheocyte
  "Hml", "NimC1", "Ppn", "Pxn", "pnr", "vkg", "zfh1", # hemocyte
  "AkhR", "CG14990", "CG34166", "CG6415", "FASN1", "Fbp1", "Fbp2", "LRP1", "Lsp1alpha", "Lsp1beta", "Lsp2", "Lst", "apolpp", "edin", "path", "srp", # adult fat body
  "Abd-B", "Act57B", "CG32121", "Mf", "Mhc", "Mlc1", "Mlc2", "Mlp84B", "Neto", "Nlg1", "Prm", "Sh", "Shab", "Strn-Mlck", "TpnC41C", "TpnC73F", "Zasp52", "Zasp66", "bt", "mspo", "nord", "nrm", "sls", "sqa", "tn", "up", "wupA", # muscle cells
  "5-HT2A", "Lim1", "Oamb", "Rdl", "Syt1", "VAChT", "brp", "futsch", "mmd", "nSyb", "nrv3", "para", "retn", "elav", # neuron
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
  "c-cup", "orb", "p-cup", "soti", "whip"# mid-late elongation-stage spermatid
)


our_curated_markers <- c(
  "org-1","Fas3",# Hub cells
  "InR","aop","Socs36E", # hub cells & CySCs
  "zfh1","ptc","piwi","tj","puc","stg", # CySCs
  "Wnt4","CG31676", # early cyst cells 
  "sog","Six4","Pdcd4", # intermediate cyst cells, spermatocytes associated CCs
  "CG8665", "CG3376", # late cyst cells , Spermatid associated CCs  
  "Nep2","Nep4","mnd", # Terminal epithelial cells
  "Ubx",# "Somatic cells at base of testis
  "CG18628","MtnA",# male gonad associated epithelium
  "glob1","CG2233", # Pigment cells
  "Ppn","Hml", # Hemocytes
  "Mp20","Zasp66", # muscle cell
  "aub","zpg","vas","nos","stg","bam", # GSCs/spermatogonia
  "Rbp4","kmg","wuc","nht","Taf12L", # early spermatocytes
  "nht","Taf12L","sa","CycB","fzo","twe","CG3927", # intermediate spermatocytes
  "CycB","fzo","CG3927","bol", # late spermatocytes
  "bol","Dpy-30L2", # early spermatids
  "schuy","hale","boly","whip","p-cup","wa-cup","r-cup","d-cup" # late spermatids
)


for ( nPCs in seq(20,40,5)){
  DefaultAssay(testis_scRNA.integrated) <- "integrated"
  testis_scRNA.integrated <- FindNeighbors(testis_scRNA.integrated, dims = 1:nPCs)
  testis_scRNA.integrated <- FindClusters(testis_scRNA.integrated, resolution = 0.2)
  testis_scRNA.integrated <- RunUMAP(testis_scRNA.integrated, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(testis_scRNA.integrated, reduction = "umap",label=TRUE)
  p2 <- DimPlot(testis_scRNA.integrated, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(testis_scRNA.integrated) <- "RNA"
  pdf(str_c(out_dir,"integrated_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(testis_scRNA.integrated,features=unique(c(FlyCellAtlas_markers,our_curated_markers)),combine=FALSE)
  print(p)
  dev.off()
}

## select 30 PCs
DefaultAssay(testis_scRNA.integrated) <- "integrated"
testis_scRNA.integrated <- FindNeighbors(testis_scRNA.integrated, dims = 1:30)
testis_scRNA.integrated <- RunUMAP(testis_scRNA.integrated, dims = 1:30)

## tune the resolution parameter
for (i in seq(0.4,1.2,0.2)){
  testis_scRNA.integrated <- FindClusters(testis_scRNA.integrated, resolution = i)
  pdf(str_c(out_dir,"integrated_PC30_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(testis_scRNA.integrated, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  p2 <- DimPlot(testis_scRNA.integrated, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}

## select 30 PCs, resolution 0.4
testis_scRNA.integrated <- FindClusters(testis_scRNA.integrated, resolution = 0.4)

## Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
testis_scRNA.integrated <- BuildClusterTree(testis_scRNA.integrated,dims =1:30,reorder=TRUE,reorder.numeric=TRUE)
pdf(str_c(out_dir,"integrated_PC30_resolution0.4_ClusterTree.pdf"))
PlotClusterTree(testis_scRNA.integrated,use.edge.length=FALSE)
dev.off()

pdf(str_c(out_dir,"integrated_PC30_resolution0.4_reordered_UMAP.pdf"))
DimPlot(testis_scRNA.integrated, reduction = "umap",label=TRUE)
dev.off()


## obtain cluster markers
testis_scRNA.integrated$tree.ident <- factor(testis_scRNA.integrated$tree.ident,levels=1:19)
Idents(testis_scRNA.integrated) <- testis_scRNA.integrated$tree.ident
clusters <- levels(Idents(testis_scRNA.integrated))
DefaultAssay(testis_scRNA.integrated) <- "RNA"
Cluster_Markers_ls <- lapply(1:length(unique(Idents(testis_scRNA.integrated))),function(i){
  df <- FindMarkers(testis_scRNA.integrated, ident.1 =clusters[i],assay="RNA",only.pos=TRUE)
  df %>% arrange(desc(avg_logFC),desc(pct.1))
  })


## use BgeeDB to do GO-like enrichment of anatomical terms for cluster markers, thereby facilitating cell type annotation
# reference : https://www.bioconductor.org/packages/release/bioc/vignettes/BgeeDB/inst/doc/BgeeDB_Manual.html
# https://bgee.org/?page=top_anat#/
# Fly Anatomy Ontology (FBbt) : https://github.com/FlyBase/drosophila-anatomy-developmental-ontology
# paper :The Drosophila anatomy ontology
library(BgeeDB)

# retrieve available species in the Bgee database
listBgeeSpecies()

# Create a new Bgee object
setwd("~/DB/dm6/annotation/")
bgee <- Bgee$new(species = "Drosophila_melanogaster")

# Loading calls of expression, loads a mapping from genes to anatomical structures based on calls of expression in anatomical structures
myTopAnatData <- loadTopAnatData(bgee)

# Look at the data
names(myTopAnatData)

# Prepare the gene list vector，a vector with background genes as names, and 0 or 1 values depending if a gene is in the foreground or not
GeneToAnatEntities <- read.table("~/DB/dm6/annotation/Drosophila_melanogaster_Bgee_14_1/topAnat_GeneToAnatEntities_7227_PRESENCESILVER.tsv",header=TRUE,sep="\t")
background_genes <- unique(GeneToAnatEntities$GENE_ID)

top50_geneList.ls <- lapply(1:length(unique(Idents(testis_scRNA.integrated))),function(i){
  geneList <- factor(as.integer(background_genes %in% rownames(Cluster_Markers_ls[[i]])[1:50])
  names(geneList) <- background_genes
  geneList
  })

# Prepare the topGO object
top50_myTopAnatObject.ls <- lapply(1:length(unique(Idents(testis_scRNA.integrated))),function(i){
  topAnat(myTopAnatData, top50_geneList.ls[[i]])
  })

# Launch the enrichment test
top50_results.ls <- lapply(1:length(unique(Idents(testis_scRNA.integrated))),function(i){
  runTest(top50_myTopAnatObject.ls[[i]], algorithm = 'weight', statistic = 'fisher')
  })


# Format the table of results after an enrichment test for anatomical terms
# Display results sigificant at a 1% FDR threshold
tableOver.ls <- lapply(1:length(unique(Idents(testis_scRNA.integrated))),function(i){
  df <- makeTable(myTopAnatData, top50_myTopAnatObject.ls[[i]], top50_results.ls[[i]], cutoff = 0.01)
  df %>% arrange(desc(foldEnrichment),pValue)
  })

## Based on BgeeDB result and known markers, we can annotate peripheral clusters
# cluster 1 : hemocyte (FBbt:00005916-head mesoderm derived embryonic hemocyte)
# cluster 3 : pigment cells
# cluster 13 : Male gonad associated epithelium
# cluster 19 : muscle cells (FBbt:10005369-esophageal visceral muscle)

## rough cell type annotation
lineage <- c("Hemocytes","Somatic cells","Pigment cells","Germline cells","Germline cells","Germline cells","Germline cells","Germline cells","Germline cells","Germline cells","Somatic cells","Somatic cells","Male gonad associated epithelium","Somatic cells","Somatic cells","Somatic cells","Somatic cells","Germline cells","Muscle cells")
names(lineage) <- levels(testis_scRNA.integrated)
testis_scRNA.integrated <- RenameIdents(testis_scRNA.integrated, lineage)
lineage_df <- data.frame(lineage=as.vector(Idents(testis_scRNA.integrated)),row.names=names(Idents(testis_scRNA.integrated)))
testis_scRNA.integrated <- AddMetaData(testis_scRNA.integrated,lineage_df)

## save Seurat object
saveRDS(testis_scRNA.integrated,str_c(out_dir,"integrated_PC30_testis_scRNA.rds"))

## performing somatic cells subclutering and germline cells subclustering analysis 
# somatic_subclustering.R
# germline_subclustering.R

## cell type annotation based on subclustering results







