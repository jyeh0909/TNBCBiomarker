#' @title K-M plot for genes
#' @description
#' K-M plots for a specified gene and a specified data set.
#'
#' @details
#' In this function, it will create a dataset with three variables.
#' The first column of data set is time variable, the second column is censored or not, the third column is group variable.
#' The data set is grouped by the p-value of log-rank test. \cr
#' First: We will order the gene expression.\cr
#' Second: Divide them as two group for any possibility, the minimal size of each group is 10 percent of total sample size.\cr
#' Third: Based on the group, we will have the p-values of log-rank test for each possibility \cr
#' Fourth: Order the p-values from minimal to maximum. \cr
#' Fifth: The smallest p-value will be the cutoff point. The data will be separated based on that cutoff point. If gene expression is bigger than or equal to the gene expression of cutoff point, then groups 1, otherwise groups 0. \cr
#' Then it will create a KM plot based on the group.
#'
#' @param dataset A data set. The data set used.
#' @param gene_info A character or a number. The gene information. Can be either gene name or gene ID. Gene name can be uppercase or lowercase. For example: The gene name can be "S100A8" or "s100a8".
#' @param time A character. The survival time variable name in original dataset.
#' @param event A character. The event variable name in original dataset.
#' @param receptor A logical indicating whether we will use the whole data set or receptor status = TN only. If receptor = TRUE, then use the data with receptor status = TN only; if receptor = FLASE, then use the whole data set.
#' @param criterion A number. A criterion that is used to split two groups. Range from 0 to 1.
#' @param title A character. The title of the KM plot.
#' @param fontsize A number. The font size. The default is 1.
#'
#' @examples
#' \dontrun{
#' kmplot(ispy1, gene_info = 6279, "rfs.t", "rfs.e",receptor=T,criterion=0.2, title="ISPY1",fontsize=20)
#' kmplot(ispy1, gene_info = "S100A8", "rfs.t", "rfs.e",receptor=T, title="ISPY1")
#' kmplot(ispy1, gene_info = "s100a8", "rfs.t", "rfs.e",receptor=T, title="ISPY1")
#'
#' }
#'
#' @export
#' @import survival
#' @import survminer
#' @import ggpubr


kmplot <- function (dataset, gene_info, time, event, receptor, criterion, title, fontsize){

  entrez_id <- gene_translate(gene_info)

  a <- which(colnames(dataset@phenoData@data) == time)
  b <- which(colnames(dataset@phenoData@data) == event)
  c <- which(colnames(dataset@phenoData@data) == "receptor_status")

  datasurv <- data.frame(dataset@phenoData@data[,c(a,b,c)])
  colnames(datasurv) <- c("time","event","receptor_status")

  if (sum(receptor)>0) {
    datasurv1 <- datasurv[which(datasurv$receptor_status == "TN"),]
  } else {
    datasurv1 <- datasurv
  }

  # pick the row (=gene) with the above entrez_id
  if(sum(receptor) > 0){
    gene_expr <- exprs(dataset)[which(featureNames(dataset) == entrez_id),][which(datasurv$receptor_status == "TN")]
  } else {
    gene_expr <- exprs(dataset)[which(featureNames(dataset) == entrez_id),]
  }

  if (missing(criterion)){
    criterion1 <- 0.1
  } else {
    criterion1 <- criterion
  }

  if (missing(fontsize)){
    fontsize <- 10
  } else {
    fontsize <- fontsize
  }

  datasurv2 <- cutoff(datasurv1,gene_expr,criterion1)

  fit <- survfit(Surv(time,event) ~ group, data=datasurv2)
  survplottitle <- ifelse(sum(receptor) > 0, "Data with TN only","Whole data")
  title1 <- paste0(title,", ",survplottitle)
  n1 <- sum(datasurv2$group == "High")
  n2 <- sum(datasurv2$group == "Low")
  conf <- ifelse(min(n1,n2) > 1, TRUE, FALSE)
  ggsurv <- ggsurvplot(
    fit,
    fontsize=fontsize/2.5,
    font.tickslab = fontsize,
    font.x = fontsize,
    font.y = fontsize,
    font.title = fontsize,
    font.caption = fontsize,
    title=title1,
    data=datasurv2,
    risk.table=TRUE,
    pval=TRUE,
    pval.size=fontsize/2.5,
    conf.int=conf,
    palette = c("#FF0000","#009999"),
    ggtheme=theme_minimal(),
    risk.table.y.text.col=TRUE,
    risk.table.y.text=FALSE,
    tables.theme = theme_cleantable(),
    tables.col="strata",
    risk.table.fontsize = fontsize/2.5)

  ggsurv$table <- ggsurv$table +
    theme(plot.title = element_text(size = fontsize))

  ggsurv$plot <- ggsurv$plot +
    theme(legend.text = element_text(size = fontsize))

  ggsurv1 <- ggpar(ggsurv,font.legend = list(size = fontsize,color = "black"),font.xtickslab = fontsize,font.ytickslab = fontsize)

  ggsurv1
}
