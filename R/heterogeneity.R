#' Title
#'
#' @param object seurat object
#' @param TumorName the variable name of tumor in the spot type information
#' @param celltype_var spot type column names in meta.data
#'
#' @return
#' @export
#'
#' @examples
Heterogeneity <- function(object=object,TumorName=TumorName,celltype_var=celltype_var){
  chrSample <- FetchData(object =object,vars = 'orig.ident')
  Sampleid <- unique(chrSample[,1])

  degree <- sapply(Sampleid,function(n){
    sample.object <- subset(object,subset=orig.ident ==n)
    cat(paste0('>>>',Sys.time(),':    ',"subset ",n,'\n'))

    index <- sapply(sample.object@images, function(x){
      name=NULL
      if (nrow(x@coordinates)==0) {name=TRUE}else{name=FALSE}
      return(name)
    })
    sample.object@images[index]=NULL

    ###transcriptome diversity degree
    suppressMessages(sample.object <-SCTransform(sample.object,assay = 'Spatial'))
    cat(paste0('>>>',Sys.time(),':    ',"calculate transcriptome diversity degree\n"))
    chrCelltype <-FetchData(object =sample.object,vars = celltype_var)
    indexCellname <- chrCelltype[,1] %in% TumorName
    tumorCellname <- rownames(chrCelltype)[indexCellname]

    varFeatures <- sample.object@assays[["SCT"]]@var.features
    matFeatures <- sample.object@assays$SCT@data[varFeatures,] %>% as.data.frame()
    matFeatures <- matFeatures[,tumorCellname] %>% na.omit()


    cor <- cor(matFeatures)
    valueCoefficient <- sapply(1:ncol(cor), function(x){
      indexValue <- (x+1):ncol(cor)
      value <- cor[x,][indexValue] %>% as.numeric()
      return(value)
    }) %>% unlist() %>% na.omit()
    MAD <- median(abs(valueCoefficient-median(valueCoefficient)))*1.4826

    return(c('sample'=n,'transcriptome diversity'=MAD))
  })
  return(degree)
}
