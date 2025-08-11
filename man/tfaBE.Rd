\docType{data}
\name{tfaBE}
\alias{tfaBE}
\title{Estimated TFA values for Barrett's Esophagus dataset}
\format{Matrices containing the TFA-estimates for colon, esophagus and
  gastric specific TFs across 2009 cells, plus the anatomical location
  information of each cell.  
}

\usage{
data(tfaBE)
}
\description{
There are 3 matrices of TFA-values, one for each set of tissue-specific
TFs, across the same 2009 cells. The anatomical location of each cell is
specified as a vector of length 2009.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{tfaCOLON.m}: the TFA-matrix for colon-specific TFs with
  rows labeling TFs (annotated to NCBI identifiers) and columns labeling
  cells.
  \item \code{tfaESOPH.m}: the TFA-matrix for esophageal-specific TFs.
  \item \code{tfaGASTRIC.m}: the TFA-matrix for gastric-specific TFs.
  \item \code{phenoBE.v}: a vector specifying the anatomical location of
  each cell.
}
}
\references{
Owen RP et al. Single cell RNA-seq reveals profound transcriptional similarity between Barrett's oesophagus and oesophageal submucosal glands. Nat.Commun.2018
}
\keyword{datasets}
