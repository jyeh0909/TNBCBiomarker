% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorPlot.R
\name{cor_plot}
\alias{cor_plot}
\title{Correlation plot for interested genes}
\usage{
cor_plot(dataset, gene_info, gene, receptor, number.cex, title, tl.cex,
  cl.cex)
}
\arguments{
\item{dataset}{A data set. The data set used.}

\item{gene_info}{A character or a number. The gene information. Can be either gene name or gene ID. Gene name can be uppercase or lowercase. For example: The gene name can be "S100A8" or "s100a8".}

\item{gene}{A character. The gene's name}

\item{number.cex}{A number. The font size of correlation plot. If not specify, then the default is 1.}

\item{title}{A character. The title of the correlation plot.}

\item{tl.cex}{A number. The font size of diagonal text label. If not specify, then the default is 1.}

\item{cl.cex}{A number. Cex of number-label in colorlabel. If not specify, then the default is 1.}
}
\description{
Using mixed methods to visualize a correlation matrix for each data set and the interested genes.
}
\details{
We will visulize the correaltion between interested genes by using corrplot.mixed function.
}
\examples{
\dontrun{
cor_plot(chin,c(5294,5290,5133),c("PIK3CG", "PIK3CA", "PDCD1"),receptor = FALSE, number.cex = 2, title="CHIN", tl.cex = 2)
}

}
\seealso{
\link{oneboxplot}

\link{corrplot.mixed}
}
