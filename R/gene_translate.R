#' @title Gene translator
#' @description
#' It can automatically translate gene names to gene ID.
#'
#' @details
#' When you type in no matter gene ID or gene name, this function will automatically give you gene ID eventually.
#'
#' @param gene_info A character or a number. The gene information. Can be either gene name or gene ID. Gene name can be uppercase or lowercase. For example: The gene name can be "S100A8" or "s100a8".
#'
#' @examples
#' \dontrun{
#' gene_translate(6279)
#' gene_translate("s100A8")
#' gene_translate("S100A8")
#'
#' }
#'
#' @export
#' @import AnnotationDbi
#' @import org.Hs.eg.db

gene_translate <- function(gene_info){

  gene_name <- ifelse(is.character(gene_info),TRUE,FALSE)
  gene_info <- ifelse(is.character(gene_info), toupper(gene_info) ,gene_info)
  entrez_id <- ifelse(sum(gene_name) > 0, as.numeric(unlist(mget(x=gene_info, envir=org.Hs.egALIAS2EG))), gene_info)

  return(entrez_id)
}
