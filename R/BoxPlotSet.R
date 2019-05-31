#' @title Combined Box plots for all data sets
#' @description
#' Combine box plots in five different data sets.
#' @details
#' We will use function to put boxplots together for five different data sets.
#'
#' @param gene_info A character or a number. The gene information. Can be either gene name or gene ID. Gene name can be uppercase or lowercase. For example: The gene name can be "S100A8" or "s100a8".
#' @param title A character. The description that should be put on the combined boxplot. The example uses GENE name.
#' @param ana TRUE or FALSE. Whether We should put p-value information on the plot. If ana = T then the p-vlaue info will be placed on the boxplot, if ana=FALSE, the p-value info will not be placed on the boxplots.
#' @param fontsize_plot A number. The font size of plot(title, caption, legend title ... and so on)
#' @param fontsize_pvalue A number. The font size of p-value (comparison of two groups and Omnibus p-value) on the plot.
#' @examples
#' \dontrun{
#' boxplot_combine(29126,"CD274",ana=TRUE,fontsize_plot = 10, fontsize_pvalue = 5, fontsize_legend = 10)
#' boxplot_combine("CD274","CD274",ana=FALSE,fontsize_plot = 10, fontsize_pvalue = 5, fontsize_legend = 10)
#' }
#'
#'
#' @export
#' @import corrplot
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @import lattice

#####Put the boxplots together

boxplot_combine <- function(gene_info,title,ana,fontsize_plot,fontsize_pvalue,fontsize_legend){

  entrez_id <- gene_translate(gene_info)
  datasetname <- c("tcga.brca","chin","ispy1","yau","gse25066")

  yes <- c(sum(which(featureNames(tcga.brca) == entrez_id)),
           sum(which(featureNames(chin) == entrez_id)),
           sum(which(featureNames(ispy1) == entrez_id)),
           sum(which(featureNames(yau) == entrez_id)),
           sum(which(featureNames(gse25066) == entrez_id)))

  which <- which(yes != 0)
  n <- length(which(yes != 0))

  pset <- list()
  for (i in 1:n){
    trans <- ifelse(which[i] == 1 || which[i] == 3, TRUE, FALSE)
    title1 <- toupper(datasetname[which[i]])

    dataset <- if(which[i] == 1) {
      tcga.brca
    } else if (which[i] == 2){
      chin
    } else if (which[i] == 3){
      ispy1
    } else if (which[i] == 4){
      yau
    } else {
      gse25066
    }

    pset[[i]] <- oneboxplot(dataset,entrez_id,title1,trans=trans,ana=ana,fontsize_plot = fontsize_plot,fontsize_pvalue = fontsize_pvalue, fontsize_legend = fontsize_legend)

  }

  pvalue_combine <- matrix(NA,n,4)
  N_combine <- matrix(NA,n,3)

  for(i in 1:n){
    pvalue_combine[i,] <- unlist(pset[[i]]$pvalue[1,])
    N_combine[i,] <- unlist(pset[[i]]$N[1,])
  }

  pvalue_combine <- data.frame(round(pvalue_combine,3))
  N_combine <- data.frame(N_combine)

  colnames(pvalue_combine) <- c("Omnibus", "HER2 VS ERPR", "ERPR VS TN", "HER2 VS TN")
  colnames(N_combine) <- c("HER2", "ERPR", "TN")

  row.names(pvalue_combine) <- datasetname[which]
  row.names(N_combine) <- datasetname[which]

  pset_plot <- list()
  for(i in 1:n){
    pset_plot[[i]] <- pset[[i]]$plot
  }


  #plot_combine <- do.call(grid.arrange,c(pset_plot, nrow=1,top=title,bottom="Receptor Type",left="Median-Centered Log2 Expression"))
  plot_combine <- do.call(arrangeGrob,c(pset_plot, nrow=1))
  title1 = textGrob(title, gp=gpar(fontsize=fontsize_plot))
  title2 = textGrob("Receptor Type", gp=gpar(fontsize=fontsize_plot))
  title3 = textGrob("Median-Centered Log2 Expression/Gene Expression",rot=90, gp=gpar(fontsize=fontsize_plot))
  plot_combine1 <- grid.arrange(plot_combine,top=title1,bottom=title2,left=title3)


  return(list(pvalue = pvalue_combine, N=N_combine))

}
