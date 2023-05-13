#' Title
#'
#' @param meta
#' @param tumor_name
#' @param nearby_name
#' @param group_var
#' @param celltype_var
#'
#' @return
DetectSpots <- function(meta =meta,tumor_name=tumor_name,nearby_name=nearby_name,
                        group_var=group_var,celltype_var){


  suppressMessages(library(tidyverse))
  suppressMessages(library(viridis))

  grouptype <- meta[,group_var]
  meta$coordinate <- paste0(meta$row,'_',meta$col)

  nearby_id <- meta[,celltype_var]%in% nearby_name
  nearby_spot <- meta[nearby_id,]
  nearby_coordinate <- nearby_spot[,c('col','row')]

  cat(paste0('>>>',Sys.time(),':    ',"Sweep spot\n"))
  sweeper_res <-
    sapply(rownames(nearby_coordinate), function(i){
      A <- unique(grouptype)[1]
      B <- unique(grouptype)[2]
      spot <- nearby_coordinate[i,] %>%  mutate(`A`=0,`B`=0)
      colnames(spot)[c(3,4)] <- c(A,B)

      inter <- data.frame(row=c(-1,-1,0,0,1,1),col=c(1,-1,2,-2,1,-1)) %>%
        t() %>% as.data.frame()

      #get inter-spot coordinates
      inter_coordinate <-
        sapply(1:6,function(n){
          spot_row <- inter[1,n]+spot$row
          spot_col <- inter[2,n]+spot$col
          coordinate <- paste0(spot_row,'_',spot_col)
          return(coordinate)
        })
      inter_spot <- meta %>% filter(coordinate %in% inter_coordinate)


      if(all(tumor_name%in%names(table(inter_spot$celltype)) )){
        index=(inter_spot$celltype==tumor_name)
        counts <- table(inter_spot[index,group_var])
        if(length(names(counts))==1){
          spot[,names(counts)]=as.numeric(counts)
        }else{
          spot[,names(counts)[1]]=as.numeric(counts[1])
          spot[,names(counts)[2]]=as.numeric(counts[2])
        }

      }else{NULL}
      group1 <- as.numeric(spot[,3])
      group2 <- as.numeric(spot[,4])
      return(c('group1'=group1,'group2'=group2))
    })

  sweeper_ct <-t(sweeper_res)%>%  as.data.frame()
  group1 <- unique(grouptype)[1]
  group2 <- unique(grouptype)[2]
  colnames(sweeper_ct) <- c(group1,group2)
  return(sweeper_ct)
}
