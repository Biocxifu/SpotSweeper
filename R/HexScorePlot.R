#' Title
#'
#' @param object seurat object
#' @param color
#' @param align
#' @param ncol
#' @param legend legend position
#' @param size
#' @param alpha
#' @param starshape
#' @param plot.image  logical value.
#' @param common.legend logical value.
#' @param limits color range
#' @param ScoreType The type of score you want to show. Include "ssGSEA" , "AddModuleScore"
#' @param type cell type
#'
#' @return
#' @export
#'
#' @examples
HexScorePlot <- function(object = object,type=type.by,color=NULL,align='v',
                          size=0.8,alpha=1,ncol = 1,legend = "right",starshape=6,
                          ScoreType='ssGSEA',plot.image=TRUE,common.legend=FALSE,
                          limits = NULL){
  if (ScoreType=='ssGSEA') {
    score_meta <-as.matrix(object@assays[["ImmuneScore"]]$ssGSEA)[,type]  %>% as.data.frame()
  }
  if (ScoreType=='AddModuleScore') {
    score_meta <- as.matrix(object@assays[["ImmuneScore"]]$AddModuleScore)[,type]%>% as.data.frame()
  }
  colnames(score_meta)=type
  object <-  AddMetaData(object,score_meta)
  color <- color %||% rev(brewer.pal(11,'Spectral'))
  HexSpatialPlot(object = object,group.by = type,common.legend = F,color = color,size=size,
                 align =align ,ncol = ncol,legend = legend,starshape=starshape,limits = limits)
  }
