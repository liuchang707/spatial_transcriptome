setwd("/storage/peiweikeLab/guochenyu/Liuchang/project2/02.STtools/liver/02.all_liver/03.work/05.bin300.normalliver.filtergenematrix")
.libPaths('/storage/peiweikeLab/zhangmingyuan/Rlibrary/')

library('Seurat')
library('Matrix')
library('ggplot2')
library('patchwork')
library('dplyr')
library('RColorBrewer')
########get normal td matrix
datadir<-"/storage/peiweikeLab/guochenyu/Liuchang/project2/02.STtools/liver/02.all_liver/02.binsize300_all_C2/"
Barcodes <- read.csv(paste(datadir,"collapsedBarcodes.csv",sep=""))
Genes <-  read.csv(paste(datadir,"collapsedGenes.csv",sep=""))
data <- readMM(paste(datadir,file = "collapsedMatrix.mtx",sep="")) 

colnames(data) = Barcodes$x  ## X x
rownames(data) = Genes$x    ## X x xiaoxie

wantedgene<- normal10um@assays[["Spatial"]]@counts@Dimnames[[1]]
#write.table(wanted,file="wanted.gene",sep="\t")
temp<-read.table('normal.wanted.collapsedBarcodes.table',header=TRUE,sep='\t')
wantedtile<-temp$name
td<-read.table('td.wanted.collapsedBarcodes.table',header=TRUE,sep='\t')
tdtile<-td$name

normalfiltered<-data[wantedgene,wantedtile]
tdfiltered<-data[wantedgene,tdtile]
source('seurat_to_spatial.R')

normal <- seurat_to_spatial(normalfiltered,temp$x_miseq_expand,temp$y_miseq_expand)

VlnPlot(normal, features = "nCount_Spatial", pt.size = 0.1, y.max = 2000) + NoLegend()
VlnPlot(normal, features = "nFeature_Spatial", pt.size = 0.1, y.max = 1000) + NoLegend()
normal <- subset(normal, subset = nFeature_Spatial > 250)
SpatialFeaturePlot(normal, features = "nCount_Spatial") + theme(legend.position = "right")
saveRDS(normal,"normal.RDS")
normal <- SCTransform(normal, assay = "Spatial", verbose = FALSE)
saveRDS(normal,"normal.RDS")

my<- readRDS('normal.stransform.RDS')
my <- RunPCA(my, assay = "SCT", verbose = FALSE)
DefaultAssay(object = my) <- "SCT"
my <- FindNeighbors(my, reduction = "pca", dims = 1:30)
my<-FindClusters(my, resolution = 0.5)
my <- RunUMAP(my, reduction = "pca", dims = 1:30)
DimPlot(my, reduction = "umap", label = TRUE)
DefaultAssay(object = my) <- "Spatial"
FeaturePlot(my, features = c("Mup20","Cyp2e1","mt-Co1","Mup10","Aqp1","Cd5l","Gstp1","Eef1a1","Dcn","Saa3","Malat1","Hba-a1" ))
want.gene<-c("Mup20","Cyp2e1","mt-Co1","Mup10","Aqp1","Cd5l","Gstp1","Eef1a1","Dcn","Saa3","Malat1","Hba-a1")
want1.gene<-c("Adgre1","Clec4f","Cd5l","F4/80","Timd4")
DotPlot(my, features = want.gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
DotPlot(my, features = want1.gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
pp.gene<-c("Tat","Adh1","Uox","Aldob","Pck1")
DotPlot(my, features = pp.gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
want3.gene<-c("Glul", "Cdh1")
DotPlot(my, features = want3.gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
new.cluster.ids <- c("Hep_PP", "Hep_MT", "Hep_PC", "Hep_Mup10", "Macrophage", "ENDO",
                     "Hep_Eef1a1", "Hep_Gstp", "HSC","Hep_Nuc","Hep_Injured","RBC")
names(new.cluster.ids) <- levels(my)
my <- RenameIdents(my, new.cluster.ids)
DotPlot(my, features = want.gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

