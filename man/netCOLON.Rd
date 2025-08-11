\docType{data}
\name{netCOLON}
\alias{netCOLON}
\title{Colon-specific regulatory network}
\format{A matrix with columns labeling colon-specific TFs
and rows labeling gene-targets, both annotated to Entrez gene IDs}

\usage{
data(netCOLON)
}
\description{
Each column labels an colon-specific TFs, and each column entry is
either a 0, +1 or -1, depending on whether the corresponding gene
(labeled by row) is not regulated, positively or negatively regulated by
the TF, respectively. Thus, the non-zero entries within each column
specify the regulon of the corresponding TF.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{netCOLON.m}: the regulon matrix
}
}
\references{
Teschendorff AE and Wang N. Improved detection of tumor suppressor
events in single-cell RNA-Seq data. NPJ Genomic Medicine 2020 Oct 7;5:43.
doi: 10.1038/s41525-020-00151-y.
}
\keyword{datasets}
