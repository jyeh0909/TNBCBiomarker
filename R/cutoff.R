#' @title cutoff
#' @description
#' Group the data based on log rank test.
#'
#' @details
#' The data set is grouped by the p-value of log-rank test. \cr
#' First: We will order the gene expression.\cr
#' Second: Divide them as two group for any possibility, the minimal size of each group is 10 percent of total sample size.\cr
#' Third: Based on the group, we will have the p-values of log-rank test for each possibility \cr
#' Fourth: Order the p-values from minimal to maximum. \cr
#' Fifth: The smallest p-value will be the cutoff point. The data will be separated based on that cutoff point. If gene expression is bigger than or equal to the gene expression of cutoff point, then groups 1, otherwise groups 0. \cr
#' Then it will create a KM plot based on the group.
#'
#' @param datasurv A dataset. A dataset contains time and event information.
#' @param gene_expr A vector. A vector contains gene expression information.
#' @param criterion A number. A criterion that is used to split two groups. Range from 0 to 1.

cutoff <- function(datasurv, gene_expr, criterion){

  logrank_pvalues <- NULL
  datasurv$expr_gt_cutoff <- NULL

  if (missing(criterion)){
    criterion1 <- 0.1
  } else {
    criterion1 <- criterion
  }

  min_no_of_genes <- round(length(gene_expr)*criterion)

  for (i in min_no_of_genes:(length(gene_expr)-(min_no_of_genes-1))){

    datasurv$expr_gt_cutoff <- gene_expr >= gene_expr[order(gene_expr)][i]
    tmp_result <- survdiff(Surv(time, event) ~ expr_gt_cutoff, data=datasurv)
    logrank_pvalues[order(gene_expr)[i]] <- 1 - pchisq(tmp_result$chisq, length(tmp_result$n) - 1)
  }

  logrank_pvalues[is.na(logrank_pvalues)] <- 1

  gene_ind <- which(logrank_pvalues == min(logrank_pvalues))

  datasurv$expr_gt_cutoff <- gene_expr >= gene_expr[gene_ind]

  datasurv$group <- ifelse(datasurv$expr_gt_cutoff == "TRUE", "High","Low")
  return(datasurv)
}
