library(rjson)
library(Seurat)
seurat_to_spatial <- function(mat,imgx,imgy){
  spatialObj <- CreateSeuratObject(counts =mat,assay='Spatial')
  #spatialObj$imgx <- (spatialObj$x-min(spatialObj$x))/20
  spatialObj$imgx <- as.numeric(imgx)
  spatialObj$imgy <- as.numeric(imgy)
  #spatialObj$imgx <- gsub('_.*','',colnames(spatialObj))
  #spatialObj$imgy <- gsub('.*_','',colnames(spatialObj))
  length <- max(spatialObj$imgx)- min(spatialObj$imgx)
  width <- max(spatialObj$imgy)-min(spatialObj$imgy)
  tissue_positions_list <- data.frame(row.names = colnames(spatialObj),tissue = 1,row = spatialObj$imgx,
                                      col = spatialObj$imgy,imagerow = spatialObj$imgx, imagecol = spatialObj$imgy)
  img <- new(Class = 'VisiumV1',image = matrix(1,length,width),
             scale.factors = scalefactors(spot = 1,fiducial = 1,hires = 1,lowres = 1),
             coordinates = tissue_positions_list,
             spot.radius = 1 / mean(c(length,width)))
  DefaultAssay(img) <- 'Spatial'
  spatialObj[['slice1']] <- img
  spatialObj
}
