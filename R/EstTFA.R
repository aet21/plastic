#' @title 
#' Estimation of Transcription Factor Regulatory Activity (TFA)
#' 
#' @aliases EstTFA
#'  
#' @description 
#' This function takes as input a log-normalized scRNA-Seq data matrix and a
#' set of TF-regulons and estimates for each TF the corresponding TF
#' regulatory activity level (TFA) in each cell.
#' 
#' @param data
#' A scRNA-Seq data matrix with rows labeling genes and columns 
#' labeling single cells and with rownames annotated to a gene-identifier.
#' scRNA-Seq data matrix should have undergone prior QC to remove poor quality
#' cells and each cell normalized by library size. If data has not been log-transformed,
#' the function will log2-transform with a pseudocount of +1.
#' 
#' @param reg.m
#' The regulon matrix with rows labeling gene-targets (same gene-identifier as
#' in \code{exp.m}) and columns labeling TFs. Each column has entries +1, -1 or 0,
#' depending on whether the gene is positively, negatively or not regulated by the
#' given TF.
#'
#' @param norm
#' Type of normalization to use if data is provided as a matrix.
#' "z" stands for z-score transformation, "c" just row centres the data matrix.
#'
#' @param ncores
#' Number of parallel cores to use
#' 
#' @return tfa.m
#' The estimated TFA-matrix with rows labeling TFs and columns labeling cells
#'  
#' @references 
#' Teschendorff AE., Wang N.
#' \emph{Improved detection of tumor suppressor events in single cell RNA-Seq data} 
#' NPJ Genom Med. 2020 Oct 7;5:43. doi:\href{https://doi.org/10.1038/s41525-020-00151-y}
#' 
#' 
#' 
#' @examples 
#' 
#' @importFrom parallel mclapply
#' 
#' @export
#'     

EstTFA <- function(data,reg.m,norm=c("z","c"),ncores=4){

   common.v <- intersect(rownames(data), rownames(reg.m))
   map1.idx <- match(common.v, rownames(data))
   map2.idx <- match(common.v, rownames(reg.m))
   ndata <- data[map1.idx, ] - rowMeans(data[map1.idx, ])
   if (norm == "z") {
      sd.v <- apply(data[map1.idx, ], 1, sd)
      nz.idx <- which(sd.v > 0)
      ndata[nz.idx,] <- ndata[nz.idx,]/sd.v[nz.idx]
   }
   idx.l <- as.list(1:ncol(data))
   prl.o <- mclapply(idx.l, InferTFactPRL, ndata, reg.m[map2.idx,], mc.cores = ncores)
   tfa.m <- matrix(unlist(prl.o), nrow = ncol(reg.m), ncol = length(prl.o), byrow = FALSE);
   rownames(tfa.m) <- colnames(reg.m)
   colnames(tfa.m) <- colnames(data)

   return(tfa.m);   
    
} ## EOF


#### Auxilliary functions
InferTFactPRL <- function(idx,tmp.m,reg.m){
  exp.v <- tmp.m[,idx];
  act.v <- apply(reg.m,2,function(tmp.v){lm.o <- lm(exp.v ~ tmp.v); act <- summary(lm.o)$coeff[2,3];return(act);})
  return(act.v);
} ### end of InferTFactPRL function



