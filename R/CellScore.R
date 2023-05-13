#' Title
#'
#' @param object seurat object
#' @param genelist genelist for calculating score
#' @param parallel.sz number of kernels for calculation
#' @param method Method for calculating score. Include "ssGSEA" , "AddModuleScore"
#' @param abs.ranking logical value. Default TRUE
#'
#' @return
#' @export
#'
#' @examples
CellScore <- function(object=object,genelist=NULL,method=method,parallel.sz,abs.ranking=TRUE){

  if(method=='AddModuleScore'){
    object.Module <- AddModuleScore(object =object,features = immune_list,
                                    name =paste0('AddModuleScore_',names(immune_list)))
    index <- grepl('AddModuleScore_',colnames(object.Module@meta.data))

    AddModuleScore_mat <- object.Module@meta.data[,index]
    object@assays[["ImmuneScore"]]$AddModuleScore<-AddModuleScore_mat
    return(object)
  }

  if(method=='ssGSEA'){
    exp.mat <- object@assays$SCT@counts %>% as.matrix()
    if(is.null(genelist)){
      genelist.file <- system.file("data", "CellReports_genelist.rda", package = "SpotSweeper")
      load(genelist.file)
      genelist <- CellReports_genelist
      #genelist <- apply(genelist ,1,FUN = function(x){
      # gene <- unique(x)
      # gene <- ifelse(gene %in% '',NA,gene) %>% na.omit() %>% as.character()
      #})
      message('no genelist input\nUse default genelist')
    }else{genelist=genelist}
    ssgsea<- GSVA::gsva(expr =exp.mat ,gset.idx.list = genelist, method='ssgsea',
                        kcdf='Poisson',abs.ranking=abs.ranking,parallel.sz=parallel.sz)

    #object@assays[["ImmuneScore"]]=NULL
    object@assays[["ImmuneScore"]]$ssGSEA<-as.data.frame(t(ssgsea))
    return(object)
  }

}
