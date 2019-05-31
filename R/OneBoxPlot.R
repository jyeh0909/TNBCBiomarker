#' @title Box plots by receptor types
#' @description
#' Comparing the gene expresstions based on receptor types by showing each boxplot together.
#' @details
#' This function will create boxplots based on receptor type (HER2, ERPR, TN) for a gene and a data set. Besides boxplots it also
#' gives us the t-test p-value of possible receptor type combinations and the ANOVA p-value of all receptor types. The p-values,
#' total N, and the means will automatically show on the plot.
#' @param dataset A data set. The data set used in the boxplots.
#' @param gene_info A character or a number. The gene information. Can be either gene name or gene ID. Gene name can be uppercase or lowercase. For example: The gene name can be "S100A8" or "s100a8".
#' @param title A character. The boxplot's title, can be specified as gene's name
#' @param trans TRUE or FALSE. Do log transformation or not. If trans = TRUE, then the gene expression will be tranformed to the base 2 logarithm. If trans = FLASE, then it will use original gene expression values.
#' @param ana TRUE or FALSE. Whether We should put p-value information on the plot. If ana = T then the p-vlaue info will be placed on the boxplot, if ana=FALSE, the p-value info will not be placed on the boxplots.

#' @examples
#' \dontrun{
#' oneboxplot(yau,5293,"YAU",trans = F, ana = T, fontsize_plot = 10, fontsize_pvalue = 5, fontsize_legend = 10)
#' oneboxplot(yau,5293,"YAU",trans = F, ana = F, fontsize_plot = 10, fontsize_pvalue = 5, fontsize_legend = 10)

#' }
#' @seealso \link{cor_plot}
#' @export
#' @import corrplot
#' @import ggplot2
#' @import ggpubr

#####One boxplot
oneboxplot <- function(dataset,gene_info,title,trans,ana,fontsize_plot,fontsize_pvalue, fontsize_legend) {

  entrez_id <- gene_translate(gene_info)

  gene_expression <- data.frame(expression=exprs(dataset)[which(featureNames(dataset) == entrez_id),],
                                receptor_status=as.character(dataset@phenoData@data$receptor_status))

  # Add 1 before log2 transformation to avoid -Inf
  gene_expression$logexpr <- log2(gene_expression$expression+1)

  # Reorder the levels of receptor_status so that HER2, ERPR and TN will be
  # plotted in that order
  gene_expression$receptor_status <- factor(gene_expression$receptor_status,
                                            levels=c("HER2", "ERPR", "TN"))

  # Center log2 transformed expression for each group
  gene_expression$logexpr_centered <- gene_expression$logexpr -
    median(gene_expression$logexpr, na.rm=T)

  # In function "give_n" you need to provide values for "y" and "label", which
  # are the aesthetics used for geom="text"
  # offset the label by 10% of the minimum value
  give_n <- function(x){
    return(c(y = min(x)-abs(min(x))/10, label = length(x)))
  }

  #mark postions
  t1 <- max(aggregate(gene_expression[, ifelse( sum(trans)>0 , "logexpr_centered" , "expression")], list(gene_expression$receptor_status), mean)$x) + 0.2
  t2 <- ifelse(sum(trans)>0,t1 + 0.18,t1 + 0.45)
  t3 <- ifelse(sum(trans)>0,t2 + 0.18,t2 + 0.45)

  my_comparisons <- list( c("HER2", "ERPR"), c("ERPR", "TN"), c("HER2", "TN") )

  fontsize_plot <- ifelse(missing(fontsize_plot), 10, fontsize_plot)
  fontsize_pvalue <- ifelse(missing(fontsize_pvalue), 3, fontsize_pvalue)
  fontsize_legend <- ifelse(missing(fontsize_legend), fontsize_plot, fontsize_legend)

  if(sum(ana) > 0) {
    plot <- ggboxplot(data=subset(gene_expression, !is.na(receptor_status)), x = "receptor_status", y = ifelse(sum(trans) > 0, "logexpr_centered" , "expression"),
                      color = "receptor_status", palette = "jco",lwd=1.5) + geom_jitter(width=0.05, size=1,color="grey") +
      stat_summary(fun.data = give_n, geom = "text") +
      stat_summary(fun.y=mean, color="darkred", geom="point",shape=18, size=3) + ggtitle(title) +
      stat_compare_means(method = "t.test", comparisons = my_comparisons, label.y = c(t1, t2, t3),size = fontsize_pvalue)+
      stat_compare_means(aes(label = paste0("Omnibus p = ", ..p.format..)), method = "anova", label.x.npc = "left", label.y.npc = "top",size = fontsize_pvalue) +
      xlab(" ") + ylab(" ") + labs(color = "Receptor Type")
    plot1 <- plot +
      font("title", size = fontsize_plot)+
      font("subtitle", size = fontsize_plot)+
      font("caption", size = fontsize_plot)+
      font("xlab", size = fontsize_plot)+
      font("ylab", size = fontsize_plot)+
      font("xy.text", size = fontsize_plot)+
      font("xy.title", size = fontsize_plot)+
      font("axis.title", size = fontsize_plot)+
      font("legend.title", size = fontsize_legend)+
      font("xy.text",size = fontsize_plot)+
      font("legend.text", size = fontsize_legend)
  } else {
    plot <- ggboxplot(data=subset(gene_expression, !is.na(receptor_status)), x = "receptor_status", y = ifelse(sum(trans) > 0, "logexpr_centered" , "expression"),
                      color = "receptor_status", palette = "jco",lwd=1.5) + geom_jitter(width=0.05, size=1,color="grey") +
      stat_summary(fun.y=mean, color="darkred", geom="point",shape=18, size=3) + ggtitle(title) +
      xlab(" ") + ylab(" ") + labs(color = "Receptor Type")
    plot1 <- plot +
      font("title", size = fontsize_plot)+
      font("subtitle", size = fontsize_plot)+
      font("caption", size = fontsize_plot)+
      font("xlab", size = fontsize_plot)+
      font("ylab", size = fontsize_plot)+
      font("xy.text", size = fontsize_plot)+
      font("legend.title", size = fontsize_legend)+
      font("legend.text", size = fontsize_legend)
  }


  pvalueHE <- ifelse(sum(trans)>0, (t.test(logexpr_centered ~ receptor_status, data=gene_expression[which(gene_expression$receptor_status == "HER2" | gene_expression$receptor_status == "ERPR"),]))$p.value,(t.test(expression ~ receptor_status, data=gene_expression[which(gene_expression$receptor_status == "HER2" | gene_expression$receptor_status == "ERPR"),]))$p.value)
  pvalueET <- ifelse(sum(trans)>0, (t.test(logexpr_centered ~ receptor_status, data=gene_expression[which(gene_expression$receptor_status == "TN" | gene_expression$receptor_status == "ERPR"),]))$p.value,(t.test(expression ~ receptor_status, data=gene_expression[which(gene_expression$receptor_status == "TN" | gene_expression$receptor_status == "ERPR"),]))$p.value)
  pvalueHT <- ifelse(sum(trans)>0, (t.test(logexpr_centered ~ receptor_status, data=gene_expression[which(gene_expression$receptor_status == "HER2" | gene_expression$receptor_status == "TN"),]))$p.value,(t.test(expression ~ receptor_status, data=gene_expression[which(gene_expression$receptor_status == "HER2" | gene_expression$receptor_status == "TN"),]))$p.value)
  pvalueall <- ifelse(sum(trans)>0,summary(aov(logexpr_centered ~ receptor_status, data=gene_expression))[[1]][1,5],summary(aov(expression ~ receptor_status, data=gene_expression))[[1]][1,5])
  pvalue <- list(pvalueall, pvalueHE, pvalueET, pvalueHT)

  nH <- length(which(gene_expression$receptor_status == "HER2"))
  nE <- length(which(gene_expression$receptor_status == "ERPR"))
  nT <- length(which(gene_expression$receptor_status == "TN"))
  N <- c(nH,nE,nT)

  pvalue <- matrix(pvalue,1,4)
  N <- matrix(N,1,3)

  colnames(pvalue) <- c("Omnibus", "HER2 VS ERPR", "ERPR VS TN", "HER2 VS TN")
  row.names(pvalue) <- "P-value"
  colnames(N) <- c("HER2", "ERPR", "TN")
  row.names(N) <- "n"

  return(list(plot = plot1, pvalue = pvalue, N=N))


}


