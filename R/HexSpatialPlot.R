#' Title
#'
#' @param object
#' @param group.by
#' @param color
#' @param align
#' @param ncol
#' @param legend
#' @param size
#' @param alpha
#' @param starshape
#' @param plot.image
#' @param common.legend
#' @param limits
#'
#' @return
#' @export
#'
#' @examples
HexSpatialPlot <- function(object = object,group.by=group.by,color=NULL,align='v',
                           size=0.8,alpha=1,ncol = 1,legend = "right",starshape=6,
                           plot.image=TRUE,common.legend=FALSE,limits = NULL){


  pp <- lapply(object@images,FUN = function(x){
    #image data
    image <- x@image
    img <- grid::rasterGrob(image, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
    scale.factors <- x@scale.factors[["lowres"]]

    #coordinate data
    coordinates <- x@coordinates
    coordinates$imagerow <- round(coordinates$imagerow*scale.factors,3)
    coordinates$imagecol <- round(coordinates$imagecol*scale.factors,3)
    cell <- rownames(coordinates)

    #meta data
    metadata <- FetchData(object,vars = group.by,cells = cell)
    meta <- cbind(coordinates,metadata)

    #plot
    group.name <- group.by
    group.by <- meta[,group.by]

    xmax <- dim(img$raster)[1] # dim pixel
    ymax <- dim(img$raster)[2]
    meta$imagerow <- xmax-meta$imagerow
    if(plot.image==TRUE){
      ppp <- ggplot(data = meta,aes(x=imagerow,y=imagecol,fill=group.by,col=group.by))+
        annotation_custom(img)+
        geom_star(starshape=starshape,size=size,alpha=alpha)+
        theme_void()+
        theme(legend.position = "top",aspect.ratio = 1,
        )+
        scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0,ymax), expand = c(0, 0))+
        coord_flip()
    }
    if(plot.image==FALSE){
      ppp <- ggplot(data = meta,aes(x=imagerow,y=imagecol,fill=group.by,col=group.by))+
        geom_star(starshape=starshape,size=size,alpha=alpha)+
        theme_void()+
        theme(legend.position = "top",aspect.ratio = 1,
        )+
        scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0,ymax), expand = c(0, 0))+
        coord_flip()
    }



    if (class(group.by)%in%c('integer','numeric')) {

      color <- color %||% viridis(9)
      ppp <- ppp+
        scale_fill_gradientn(colors = color,name=as.character(group.name),
                             limits = limits,
        )+
        scale_color_gradientn(colours =  color,name=as.character(group.name),
                              limits = limits)
    }else{
      if(is.null(color)){
        ppp <- ppp+guides(shape=FALSE,fill=FALSE,colour=guide_legend(title=group.name))
        return(ppp)
      }else{
        ppp <- ppp+
          scale_color_manual(values = color,name=as.character(group.name),limits = limits)+
          scale_fill_manual(values = color,name=as.character(group.name),limits = limits)}
    }
    return(ppp)

  })
  pp <- ggarrange(plotlist=pp, align =  align,ncol = ncol,common.legend = common.legend, legend = legend )
  return(pp)
}
