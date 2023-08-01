#' Title
#'
#' @param object
#' @param cancerType
#' @param method
#' @param parallel.sz
#'
#' @return
#' @export
#'
#' @examples
ImmuneState <- function(object=object ,cancerType=cancerType,method=method,parallel.sz=1){
  state.genelist <- system.file("data", "state.genelist.rda", package = "SpotSweeper")

  strsplit(state.genelist$cancer,'; ')
  index <- sapply(1:nrow(state.genelist),function(i){
    Type=strsplit(state.genelist$cancer[[i]],'; ')[[1]]
    if(length(Type)==1){dex <- (Type== cancerType)}
    else{dex <- !all((Type== cancerType)==FALSE)}

    return(dex)
  })
  geneset <- state.genelist[index,]


  ProImmunity.gene <- geneset$gene[grepl('Promote immunity',geneset$relation)]
  InhImmunity.gene <- geneset$gene[grepl('Inhibit immunity',geneset$relation)]

  genelist <- list('ProImmunity'=ProImmunity.gene,
                   'InhImmunity'=InhImmunity.gene)
  exp.mat <- colon@assays$SCT@counts %>% as.matrix()

  score<- GSVA::gsva(expr =exp.mat ,gset.idx.list = genelist, method='ssgsea',
                      kcdf='Poisson',abs.ranking=abs.ranking,parallel.sz=parallel.sz)
  browser()
  object@assays[["ImmuneScore"]]$ImmuneState@medthod<-as.data.frame(t(score))
  return(object)

}

