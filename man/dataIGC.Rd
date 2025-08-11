\docType{data}
\name{dataIGC}
\alias{dataIGC}
\title{Data and metadata of intestinal gastric cancer study}
\format{The matrix of estimated TFA values for gastric specific TFs, a
  vector encoding the estimated cellular plasticities and metadata
  encoded as a dataframe. Only epithelial cells are included.
}

\usage{
data(dataIGC)
}
\description{
There is one matrix of TFA-values for 47 stomach specific TFs and 32724
epithelial cells. There is one vector encoding the plasticity (CCAT)
values for the 32724 cells. And a datadrame with metacell information.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{tfaIGC.m}: the TFA-matrix for stomach-specific TFs with
  rows labeling TFs (annotated to NCBI identifiers) and columns labeling
  cells.
  \item \code{ccatIGC.v}: the CCAT plasticity estimates.
  \item \code{phenoIGC.df}: dataframe with metacell information.
}
}
\references{
Zhang P et al. Dissecting the Single-Cell Transcriptome Network
Underlying Gastric Premalignant Lesions and Early Gastric Cancer. Cell
Rep.2019, 17:1934-1947.
}
\keyword{datasets}
