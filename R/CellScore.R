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
CellScore <- function(object=object,genelist=NULL,method=method,parallel.sz=1,abs.ranking=TRUE){
  if(!method %in% c('ssGSEA','AddModuleScore')){
    stop('The input method is incorrect')
  }


  if(is.null(genelist)){

    data("CellReports_genelist",package = 'SpotSweeper')
    genelist <- CellReports_genelist

    message('no genelist input\nUse default genelist')
  }else{genelist=genelist}

  names(genelist) <- gsub(' ','_',names(genelist))


  if(method=='AddModuleScore'){
    cat('Using AddModuleScore')

    object.Module <- AddModuleScore(object =object,features = genelist,
                                    name =paste0('AddModuleScore_',names(genelist)))
    index <- grepl('AddModuleScore_',colnames(object.Module@meta.data))

    AddModuleScore_mat <- object.Module@meta.data[,index]
    object@assays[["ImmuneScore"]]$AddModuleScore<-AddModuleScore_mat
    return(object)
  }

  if(method=='ssGSEA'){
    cat('Using ssGSEA')
    exp.mat <- object@assays$SCT@counts %>% as.matrix()
    ssgsea<- GSVA::gsva(expr =exp.mat ,gset.idx.list = genelist, method='ssgsea',
                        kcdf='Poisson',abs.ranking=abs.ranking,parallel.sz=parallel.sz)

    object@assays[["ImmuneScore"]]$ssGSEA<-as.data.frame(t(ssgsea))
    return(object)
  }

}
