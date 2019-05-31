#' @title Correlation plot for interested genes
#' @description
#' Using mixed methods to visualize a correlation matrix for each data set and the interested genes.
#'
#' @details
#' We will visulize the correaltion between interested genes by using corrplot.mixed function.
#'
#' @param dataset A data set. The data set used.
#' @param gene_info A character or a number. The gene information. Can be either gene name or gene ID. Gene name can be uppercase or lowercase. For example: The gene name can be "S100A8" or "s100a8".
#' @param gene A character. The gene's name
#' @param number.cex A number. The font size of correlation plot. If not specify, then the default is 1.
#' @param title A character. The title of the correlation plot.
#' @param tl.cex A number. The font size of diagonal text label. If not specify, then the default is 1.
#' @param cl.cex A number. Cex of number-label in colorlabel. If not specify, then the default is 1.
#' @examples
#' \dontrun{
#' cor_plot(chin,c(5294,5290,5133),c("PIK3CG", "PIK3CA", "PDCD1"),receptor = FALSE, number.cex = 2, title="CHIN", tl.cex = 2)
#' }
#'
#' @seealso \link{oneboxplot}
#' @seealso \link{corrplot.mixed}
#'
#' @export
#' @import corrplot
#' @import ggplot2
#' @import gridExtra
#' @import lattice

cor_plot <- function(dataset, gene_info, gene, receptor, number.cex, title, tl.cex, cl.cex){

  entrez_id <- NULL
  for (i in 1:length(gene_info)){
    entrez_id[i] <- gene_translate(gene_info[i])
  }

  number.cex <-ifelse(missing(number.cex),1,number.cex)
  tl.cex <-ifelse(missing(tl.cex),1,tl.cex)
  cl.cex <-ifelse(missing(cl.cex),1,cl.cex)

  a <- NULL
  for(i in 1:length(entrez_id)){
    a[i] <- which(featureNames(dataset) %in% entrez_id[i])
  }

  b <- data.frame(exprs(dataset)[a,])
  c <- t(b)
  d <- data.frame(dataset@phenoData@data$receptor_status)
  row.names(d) <- row.names(dataset@phenoData@data)
  dataused1 <- cbind(c,d)
  colnames(dataused1) <- c(gene,"receptor_status")

  if (sum(receptor)>0) {
    dataused2 <- dataused1[which(dataused1$receptor_status == "TN"),]
  } else {
    dataused2 <- dataused1
  }
  a <- which(colnames(dataused2) == "receptor_status")
  dataused3 <- dataused2[,-a]

  corrplot.mixed(cor(dataused3), lower.col = "black", number.cex = number.cex, title = title, mar=c(0,0,2,0),tl.cex = tl.cex, cl.cex = cl.cex)
}
