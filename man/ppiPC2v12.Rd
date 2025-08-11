\docType{data}
\name{ppiPC2v12}
\alias{ppiPC2v12}
\title{Human PPI network from Pathway Commons PC2 version-12}
\format{Adjacency matrix for the PPI network  
}

\usage{
data(ppiPC2v12)
}
\description{
Adjacency matrix describing a curated human PPI network derived from Pathway
Commons-2 (PC2) version-12 on the 4th August 2020, defined for 12649 proteins.
These have been annotated to NCBI identifiers. Network was further
filtered using main cellular localization information from the Human
Protein Atlas. Diagonal entries are 0. Other entries are zero or 1 if an
interaction has been reported. Network is the maximally connected
component of the one after filtering.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{ppiPC2v12.m}: the adjacency matrix
}
}
\references{
Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of
differentiation potency from a cell’s transcriptome. Nat.Commun.2017
}
\keyword{datasets}
