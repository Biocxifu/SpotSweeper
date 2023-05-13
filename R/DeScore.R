#' Title
#'
#' @param object seurat object
#' @param hide.ns logical value. If TRUE, hide ns symbol when displaying significance levels.
#' @param label   label type. Include "p.signif" , "p.format" (shows the formatted p value).
#' @param ctr.group a list for comparison. list index 2 VS list index 1
#' @param test.method method to be used for comparing means
#' @param ScoreType  The type of score you want to compare. Include "ssGSEA" , "AddModuleScore"
#'
#' @return
#' @export
#'
#' @examples

DeScore <- function(object =object,test.method='wilcox.test',hide.ns = T,ctr.group=ctr.group,
                    label = 'p.signif',ScoreType='ssGSEA'){
  if (ScoreType=='ssGSEA') {cellscore <- object@assays[["ImmuneScore"]][['ssGSEA']]}
  if (ScoreType=='AddModuleScore') {
    cellscore <- object@assays[["ImmuneScore"]][["AddModuleScore"]] %>% t() %>% as.data.frame()}

  normalization<-function(x){return((x-min(x))/(max(x)-min(x)))}
  cellscore <- normalization(cellscore)

  cellscore.nearby <- cellscore[,colnames(object)] %>%
    t() %>% as.data.frame()
  cellscore.nearby$Sweeper_type <- object$Sweeper_type
  cellscore.nearby <- cellscore.nearby %>% filter(Sweeper_type %in% unlist(ctr.group))
  cellscore.nearby$Sweeper_type <- ifelse(cellscore.nearby$Sweeper_type %in% ctr.group[[1]],'group1','group2')

  levels <-c('group1','group2')
  string2 <- sapply(ctr.group[[2]],function(x){
    string=as.character(x)
    return(string)
    })
  string1 <- sapply(ctr.group[[1]],function(x){
    string=as.character(x)
    return(string)
  })

  compare <- paste0(paste(collapse='&',string2),'VS',paste(collapse='&',string1))

  df <- cellscore.nearby %>%
    reshape2::melt(id='Sweeper_type')
  df$Sweeper_type <- factor(df$Sweeper_type,levels = levels)

  p <-
    ggplot(df,aes(x = variable,y=value,color=Sweeper_type,fill=Sweeper_type))+
    geom_boxplot(alpha=0.6,outlier.alpha = 0,cex=0.3)+

    scale_color_manual(values = c('#5aadd0','#af7aa1'))+
    scale_fill_manual(values = c('#5aadd0','#af7aa1'))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30,
                                     vjust = 1,hjust = 1,
                                     size = 6,
                                     color="black")

    )+
    ggpubr::stat_compare_means(method = test.method,hide.ns = hide.ns,label = label)
  data_summary <- function(data, valuename, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(median = median(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))}
    data_sum<-ddply(data, groupnames, .fun=summary_func,valuename)
    data_sum <- rename(data_sum, c("median" = valuename))
    return(data_sum)
  }

  df_summary <- data_summary(df,valuename = 'value',groupnames =c('variable','Sweeper_type'))
  fold <- sapply(unique(df_summary$variable),function(x){
    id <- df_summary$variable %in% as.character(x)
    fold_data <- df_summary[id,]
    idx1 <- fold_data$Sweeper_type==levels[1]
    idx2 <- fold_data$Sweeper_type==levels[2]

    fold <- as.numeric(fold_data[idx2,3])/as.numeric(fold_data[idx1,3])
    return(c('Fold'=fold,'name'=x))
  })
  fold[2,]=as.character(unique(df_summary$variable))
  fold <- t(fold) %>% as.data.frame()


  fold$Fold <- as.numeric(fold$Fold)
  fold$name <- factor(fold$name,levels = levels(df$variable))

  p2 <- ggplot(fold,aes(x = name,y=Fold,col=Fold))+
    geom_segment( aes(x = name, xend =name, y = 1, yend = Fold),linetype=2,
                  cex=0.6,color = "grey30")+
    geom_hline(yintercept = 1,color='grey90')+
    geom_point(aes(size=abs(1-Fold)))+
    scale_color_gradientn(colours =  c('#00BCB4','#FFD804'))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x =  element_blank(),
          axis.ticks.x = element_blank(),axis.title.x = element_blank()
    )+
    ggtitle(compare)

  pp <- aplot::insert_top(p,p2,height = c(1,2))
  return(list('plot'=pp,'Fold'=fold))
}
