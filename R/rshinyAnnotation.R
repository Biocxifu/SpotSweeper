#' Title
#'
#' @param object
#' @param group.by
#' @param ident
#' @param metadata
#' @param color
#' @param size
#' @param alpha
#' @param shape
#' @param celltype_var
#' @param plot.image
#' @param common.legend
#' @param limits
#'
#' @return
#' @export
#'
#' @examples
singleplot <- function(object = object,group.by=group.by,ident=ident,metadata= metadata,
                       color=NULL,size=0.8,alpha=1,shape=19,celltype_var=celltype_var,
                       plot.image=TRUE,common.legend=FALSE,limits = NULL){
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(ggiraph)
  library(tidyverse)
  library(viridis)

  object <- subset(object,orig.ident %in% ident)
  plot.image <- as.logical(plot.image)
  #remove
  spatial.data<- lapply(object@images, function(i){
    if(nrow(i@coordinates)==0){}else{i}})
  spatial.data <- Filter(Negate(function(x) is.null(unlist(x))),spatial.data)


  image <- spatial.data[[1]]@image
  img <- grid::rasterGrob(image, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
  scale.factors <- spatial.data[[1]]@scale.factors[["lowres"]]

  #coordinate data
  coordinates <- spatial.data[[1]]@coordinates
  coordinates$imagerow <- round(coordinates$imagerow*scale.factors,3)
  coordinates$imagecol <- round(coordinates$imagecol*scale.factors,3)

  metadata <- data.frame(orig.ident=metadata$orig.ident,
                         celltype=metadata$celltype,
                         cellid=metadata$cellid) %>%
    column_to_rownames('cellid')

  index_ident <- (metadata$orig.ident %in% ident)
  metadata <-  metadata[index_ident,]
  metadata <- metadata[rownames(coordinates),]
  meta <- cbind(coordinates,metadata)
  #tooltip <- paste0(rownames(meta),'_',meta$celltype)


  #Feature
  group.name <- group.by
  group.by <- FetchData(object,vars =group.name )
  if(group.name == celltype_var){
    group.by <- metadata$celltype %>% as.data.frame()
  }

  names(group.by) <- group.name
  xmax <- dim(img$raster)[1] # dim pixel
  ymax <- dim(img$raster)[2]
  meta$imagerow <- xmax-meta$imagerow

  if(celltype_var != group.name){
    if(class(group.by[,1])=='numeric'){info = round(group.by[,1],3)}else
      {info = group.by[,1]}
    tooltip <-paste0('Spot type: ',meta$celltype,'\n',group.name,': ',info)
  }else{tooltip <-meta$celltype}


  if(plot.image==TRUE){
    ppp <- ggplot(data = meta,aes(x=imagerow,y=imagecol,fill=group.by[,1],col=group.by[,1]))+
      annotation_custom(img)+
      geom_point_interactive(size=size,alpha=alpha,data_id = rownames(meta),
                             aes(tooltip =tooltip),shape=shape)+
      theme_void()+
      theme(legend.position = "right",aspect.ratio = 1,
      )+
      scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0,ymax), expand = c(0, 0))+
      coord_flip()
    ppp
  }
  if(plot.image==FALSE){
    ppp <- ggplot(data = meta,aes(x=imagerow,y=imagecol,fill=group.by[,1],col=group.by[,1]))+
      geom_point_interactive(size=size,alpha=alpha,data_id = rownames(meta),
                             aes(tooltip = tooltip),shape=shape)+
      theme_void()+
      theme(legend.position = "right",aspect.ratio = 1
      )+
      scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0,ymax), expand = c(0, 0))+
      coord_flip()
  }



  if (class(group.by[,1])%in%c('integer','numeric')) {

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
}


#####shiny-----
#' Title
#'
#' @param object
#' @param celltype_var
#'
#' @return
#' @export
#'
#' @examples
SpotAnnotation <- function(object=object,celltype_var=celltype_var){
  sample <- unique(object$orig.ident)
  feautures <- c(rownames(object),unique(colnames(object@meta.data)))
  celltype_select <- unique(object@meta.data[,celltype_var])
  ui <-  pageWithSidebar(
    headerPanel('SpotAnnotation'),

    sidebarPanel(width = 3,

                 selectInput(inputId = 'orig.ident', label = 'Select sample(orig.ident)', choices = sample, selected = sample[1]),
                 selectInput(inputId = 'Feauture', label = 'Feautures & metadata', choices = feautures, selected = feautures[1]),
                 shiny::hr(),

                 selectInput(inputId = 'Input.type', label = 'Modify spot type', choices = celltype_select, selected = 'Default'),
                 actionButton(inputId = 'confirm', label='Confirm modification'),
                 shiny::hr(),

                 checkboxInput(inputId = 'plot.image',label = 'Show image',value = TRUE),
                 selectInput(inputId= 'SpotShape', label = 'Spot shape', choices = 1:25, selected = 16),
                 sliderInput(inputId='SpotAlpha', label='Spot alpha', min=0, max=1, value=0.4, step=0.1),
                 sliderInput(inputId='SpotSize', label='Spot size', min=0, max=3, value=1, step=0.1),
                 shiny::hr(),
                 HTML("<br>"),

                 actionButton(inputId = 'Stop', label='Save&Quit')
    ),
    mainPanel(
      ggiraph::ggiraphOutput('plot', width = '100%', height = paste0(1000, 'px'))
    )
  )


  server <- function(input, output,session) {
    metadata <- reactiveValues(orig.ident=FetchData(object,vars = 'orig.ident')[,1],
                               celltype=FetchData(object,vars = celltype_var)[,1],
                               cellid=colnames(object))

    output$plot <- renderGirafe({
      withProgress(message = "Updating...",value = 1,{
        p <- singleplot(object,ident = input$orig.ident,group.by = input$Feauture,shape =as.numeric(input$SpotShape),
                        size = input$SpotSize,metadata = metadata,plot.image =input$plot.image,
                        celltype_var = celltype_var,alpha = input$SpotAlpha) %>%
          girafe(ggobj = .,width_svg = 6,height_svg =4,
                 options = list(opts_tooltip(use_fill = TRUE))) %>%
          ggiraph::girafe_options(x = .,
                                     ggiraph::opts_zoom(max = 5,min = 0.5),
                                     ggiraph::opts_selection(type = "multiple",
                                                             css = "fill:cyan;stroke:cyan;opacity:1;"))
        return(p)

      })

    })
    observeEvent(input$confirm,{
      selected <- as.character(input$plot_selected)
      index <-metadata$cellid %in% selected
      metadata$celltype[index] <- input$Input.type
      cat('\nSelected spots:\n')
      cat(selected)
      cat(paste0('\nAnnotated as: ',input$Input.type,'\n'))
      session$sendCustomMessage(type = 'plot_set', message = character(0))
    })
    observe({
      if(input$Stop > 0){
        final.meta <- data.frame(orig.ident=metadata$orig.ident,
                                 celltype=metadata$celltype,
                                 cellid=metadata$cellid)

        object@meta.data[final.meta$cellid,][,celltype_var] <- metadata$celltype

        cat("\nQiut\n")
        stopApp(returnValue = object)
      }})

  }
  #shinyApp(ui, server)
  runApp(list(ui = ui, server = server), launch.browser = T)
}

