#' Title
#'
#' @param meta
#' @param color
#' @param limits
#' @param group.by
#' @param img
#' @param size
#' @param alpha
#' @param starshape
#' @param plot.image
#'
#' @return
#' @export
#'
#' @examples
HexSpatPlot <- function(meta =meta,color =NULL,limits = NULL,img=img,size=size,
                        alpha=alpha,group.by=group.by,starshape=starshape,
                        plot.image=plot.image){
  suppressMessages(library(ggstar))
  suppressMessages(library(ggplot2))


  group.name <- group.by
  group.by <- meta[,group.by]

  xmax <- dim(img$raster)[1] # dim pixel
  ymax <- dim(img$raster)[2]
  meta$imagerow <- xmax-meta$imagerow
  browser()
  if(plot.image==TRUE){
    ppp <- ggplot(data = meta,aes(x=imagerow,y=-imagecol,fill=group.by,col=group.by))+
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
    ppp <- ggplot(data = meta,aes(x=imagerow,y=-imagecol,fill=group.by,col=group.by))+
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
      return(ppp)
    }else{
      ppp <- ppp+
        scale_color_manual(values = color,name=as.character(group.name),limits = limits)+
        scale_fill_manual(values = color,name=as.character(group.name),limits = limits)}
  }
  return(ppp)
}
