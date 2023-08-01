#' Title
#'
#' @param object seurat object
#' @param group_var group column names in meta.data
#' @param celltype_var spot type column names in meta.data
#' @param tumor_name the variable name of tumor in the spot type information
#' @param nearby_name the variable name of the spot to be analyzed that is adjacent to the tumor spot
#'
#' @return
#' @export
#'
#' @examples
SpotSweeper <-function(object = object,group_var=group_var,celltype_var=celltype_var,
                       tumor_name=tumor_name,nearby_name=nearby_name) {
  SpotSweeper.res <- lapply(object@images,FUN = function(x){
    coordinates <- x@coordinates
    cell <- rownames(coordinates)
    metadata <- FetchData(object,vars = c(group_var,celltype_var),cells = cell)
    meta <- cbind(coordinates,metadata)
    sweeper_ct <- DetectSpots(meta = meta,tumor_name=tumor_name,nearby_name=nearby_name,
                              group_var=group_var,celltype_var = celltype_var)
    return(sweeper_ct)
  })
  res <- NULL
  for (i in 1:length(SpotSweeper.res)) {
    res <- rbind(SpotSweeper.res[[i]],res)
  }

  #add spot type
  metadata <- FetchData(object,vars = c(group_var,celltype_var))

  index1 <- ifelse(res[,1]>0,TRUE,FALSE)
  index2 <- ifelse(res[,2]>0,TRUE,FALSE)

  index_type <- index1&index2
  index_type[index1] <- paste0('Nearby_',colnames(res)[1])
  index_type[index2] <- paste0('Nearby_',colnames(res)[2])

  index_type <- ifelse(index1&index2,'Nearby_all',index_type)
  index_type[(index1|index2)==FALSE]='Nearby_none'

  Sweeper_type <- data.frame('Sweeper_type'=index_type,row.names = rownames(res))
  object <- AddMetaData(object ,Sweeper_type)

  Tumor_index <- metadata[,celltype_var] %in% tumor_name
  object$Sweeper_type[Tumor_index] <- metadata[,group_var][Tumor_index]
  object$Sweeper_type[is.na(object$Sweeper_type)] <- 'Others'

  return(object)
}
