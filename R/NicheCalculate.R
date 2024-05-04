#' Title Calculate the mean and log2FC value for every niche
#'
#' @param Ecotype_P
#'
#' @return
#' @export
#'
#' @examples
NicheCalculate <- function(Ecotype_P=Ecotype_P){
  #1 median
  calculate_cell_proportions_median <- function(data, niche_column, cell_columns) {
    result <- data %>%
      group_by({{ niche_column }}) %>%
      summarise(across({{ cell_columns }}, median, na.rm = TRUE))
    return(result)
  }
  #2 logFC
  logFC<- function(data) {
    p_values <- lapply(colnames(data)[-length(colnames(data))],  function(i) {
      x=data[,i]

      Niches <- unique(data$Ecotype)
      logFC <- lapply(Niches, function(y){
        group <- ifelse(data$Ecotype==y,y,'others')
        v1 <- x[group==y] %>% median() %>% as.numeric()
        v2 <- x[group=='others'] %>% median() %>% as.numeric()

        logFC <- -log2(v1/v2)

        return(logFC)
      })

      logFC <- as.numeric(logFC)
      return(c(logFC))

    })
  }

  # plot
  cellscore <- Ecotype_P[[4]]
  cell_score_median <- calculate_cell_proportions_median(cellscore, Ecotype,
                                                         colnames(cellscore)[-length(colnames(cellscore))])
  cell_score_median2 <-cell_score_median %>%  column_to_rownames('Ecotype') %>%
    apply(., 1, as.numeric)
  colnames(cell_score_median2) <- cell_score_median$Ecotype
  rownames(cell_score_median2) <- colnames(cell_score_median)[-1]
  green_palette <- colorRampPalette(c("#fbfdfc",'#89dab2',"#4fabab", "#3494a8"))
  library(viridis)

  median_heat <- Heatmap(cell_score_median2,

                         name = "score",show_column_dend = F,show_row_dend = F,
                         col = mako(100),
                         cluster_rows =T,
                         cluster_columns = T,
                         show_row_names = TRUE,
                         show_column_names =T)



  cell_logFC <- logFC(cellscore)
  cell_logFC_data <- sapply(1:(length(colnames(cell_score_median))-1),
                            FUN = function(i)(return(cell_logFC[[i]])))
  rownames(cell_logFC_data) <-unique(cell_score_median[,'Ecotype'])[[1]]
  colnames(cell_logFC_data) <- colnames(cell_score_median)[-1]
  cell_logFC_data[is.infinite(cell_logFC_data)] <- NA
  cell_logFC_data <- t(cell_logFC_data ) %>% na.omit()

  heat_palette2 <- colorRampPalette( rev(c("#c11658",'white','#aee2f8')))

  logFC_heat <- Heatmap(cell_logFC_data,
          col = colorRamp2(seq(-1, 1, length.out = 100), heat_palette2(100)),
          name = "log2FC",show_column_dend = F,show_row_dend = F,
          cluster_rows = T,
          cluster_columns = T,
          show_row_names = TRUE,
          show_column_names =T)
  return(list(median_heat,logFC_heat))

}
