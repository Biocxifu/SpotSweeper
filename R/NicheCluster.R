#' Title Determine the number of subtypes
#'
#' @param object
#' @param K.max
#' @param nstart
#' @param B
#' @param ncol
#' @param ScoreType
#'
#' @return
#' @export
#'
#' @examples
NicheKmeans <- function(object =object,K.max = 10,nstart = 25,B = 50,ncol = 2,
                        ScoreType='ssGSEA'){
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)



  if (ScoreType=='ssGSEA') {
    ImmuneScore <- object@assays[["ImmuneScore"]][['ssGSEA']] %>% t()}
  if (ScoreType=='AddModuleScore') {
    ImmuneScore <- object@assays[["ImmuneScore"]][["AddModuleScore"]] %>% t() }
  normalization<-function(x){return((x-min(x))/(max(x)-min(x)))}
  ImmuneScore <- normalization(ImmuneScore)%>% t()

  for (i in rownames(ImmuneScore)) {
    ImmuneScore[i,] <- (ImmuneScore[i,] -min(ImmuneScore[i,]))/(max(ImmuneScore[i,] )-min(ImmuneScore[i,] ))
   # ImmuneScore[i,] <- (ImmuneScore[i,])/sum(ImmuneScore[i,])
  }


  fviz_plot <- fviz_nbclust(ImmuneScore,kmeans,k.max = 10, method = "wss")
  gap_stat <- clusGap(ImmuneScore,
                      FUN = kmeans,
                      nstart = 25,
                      K.max = 10,
                      B = 50)

  #plot number of clusters vs. gap statistic
  fviz_gap <- fviz_gap_stat(gap_stat)


  pp <- plot_grid(fviz_plot,fviz_gap,ncol = ncol)
  return(pp)
}


#' Title subtype niche cluster
#'
#' @param object
#' @param centers
#' @param nstart
#' @param ScoreType
#' @param color
#'
#' @return
#' @export
#'
#' @examples
NicheCluster <- function(object =object,centers = centers, nstart = 25,
                         ScoreType='ssGSEA',color=NULL){
  if (ScoreType=='ssGSEA') {
    ImmuneScore <- object@assays[["ImmuneScore"]][['ssGSEA']] %>% t()}
  if (ScoreType=='AddModuleScore') {
    ImmuneScore <- object@assays[["ImmuneScore"]][["AddModuleScore"]] %>% t() }
  normalization<-function(x){return((x-min(x))/(max(x)-min(x)))}
  ImmuneScore <- normalization(ImmuneScore)%>% t()

  for (i in rownames(ImmuneScore)) {
    ImmuneScore[i,] <- (ImmuneScore[i,] -min(ImmuneScore[i,]))/(max(ImmuneScore[i,] )-min(ImmuneScore[i,] ))
    #ImmuneScore[i,] <- (ImmuneScore[i,])/sum(ImmuneScore[i,])
  }

  km <- kmeans(ImmuneScore, centers = centers, nstart = nstart)
  clusterP <- fviz_cluster(km, data = ImmuneScore,labelsize = 0,ellipse = F)
  clusterData <- clusterP[["data"]]

  PCA_plot <- ggplot(clusterData, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 1) +
    scale_color_manual(values =  color)+
    labs(title = "Scatter Plot of X and Y Coordinates",
         x = "X Coordinate",
         y = "Y Coordinate",
         color = "Cluster") +
    theme_classic()+
    theme(aspect.ratio = 1,axis.text = element_text(colour = 'black'))

  green_palette <- colorRampPalette(c("#fbfdfc",'#89dab2',"#3494a8", "#3494a8"))

  Heatmap_plot <- Heatmap(t(ImmuneScore),column_split = as.character(km$cluster),use_raster = F,
                          col = colorRamp2(seq(0, max(ImmuneScore), length.out = 100), green_palette(100)),
                          name = "cell proportion",show_column_dend = F,show_row_dend = F,
                          cluster_rows = TRUE,
                          cluster_columns = F,
                          show_row_names = TRUE,
                          show_column_names = F)



  ImmuneNiche <- cbind(ImmuneScore,'Ecotype'=paste0('Niche',as.character(km$cluster))) %>% as.data.frame()
  object$Ecotype <- paste0('Niche',as.character(km$cluster))
  SpatialPlot <- HexSpatialPlot(object = object, group.by = "Ecotype",plot.image = F,color = color)
  return(list(PCA_plot, Heatmap_plot,SpatialPlot,ImmuneNiche))
}
