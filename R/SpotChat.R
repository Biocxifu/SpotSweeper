#' Title
#'
#' @param object seurat object
#' @param subSample sample name. Ensure that the sample name is in the 'orig.ident' column
#' @param cores Number of cores used
#' @param json.path Sample corresponding JSON file path
#' @param min.cells
#'
#' @return
#' @export
#'
#' @examples
SpotChat <- function(object =object,subSample=subSample,cores=cores,json.path=json.path,
                     min.cells = 10){
  sample <- subset(colon,subset=orig.ident%in%subSample)

  index <- sapply(sample@images, function(x){
    name=NULL
    if (nrow(x@coordinates)==0) {name=TRUE}else{name=FALSE}
    return(name)
  })
  sample@images[index]=NULL

  # Prepare input data for CelChat analysis
  data.input = GetAssayData(sample, slot = "data", assay = "SCT") # normalized data matrix

  meta = data.frame(labels = sample$Sweeper_type, row.names = names(sample$Sweeper_type)) # manually create a dataframe consisting of the cell labels
  unique(meta$labels)

  spatial.locs = GetTissueCoordinates(sample, scale = NULL,
                                      cols = c("imagerow", "imagecol"))
  scale.factors = jsonlite::fromJSON(txt = json.path)
  scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres,
                       hires = scale.factors$tissue_hires_scalef,
                       lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )

  #creat cellchat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)

  CellChatDB <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use

  #ligand receptor interaction database
  cellchat <- subsetData(cellchat)
  future::plan("multiprocess", workers = cores)


  #Preprocessing for expression data analysis in intercellular communication
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  #Calculate communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                distance.use = TRUE, interaction.length = 200,
                                scale.distance = 0.01)
  table(sample$Sweeper_type)

  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  #Inferring intercellular communication at the level of signaling pathways
  cellchat <- computeCommunProbPathway(cellchat)

  #Computational Aggregated Cell Cell Cell Communication Network
  cellchat <- aggregateNet(cellchat)

  cat('Drawing a circle plot')
  pathways <- cellchat@netP[["pathways"]]
  pdf(file =paste0(subSample,'_circle_chat.pdf') ,width = 5,height = 5)
  lapply(1:length(pathways), function(i){
    pathways.show=pathways[i]
    p=netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
    return(p)
  })
  dev.off()
  cat('Drawing a circle plot')
  cat(paste0('output file: ',subSample,'_circle_chat.pdf'))

  saveRDS(cellchat, file = paste0("cellchat_",subSample,".rds"))

  cat(paste0('output file: ',"cellchat_",subSample,".rds"))

}
