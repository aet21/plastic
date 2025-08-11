#' @title 
#' Computation of TFAL
#' 
#' @aliases CompTFAL
#'  
#' @description 
#' This function takes as input the estimated TFA matrix and computes for each
#' cell the transcription factor activation load (TFAL), defined as the number
#' of tissue-specific TFs that are activated in a given cell compared to the
#' normal cells.
#' 
#' @param tfa.m
#' The TFA-matrix defined over tissue-specific TFs (rows) and cells (columns)
#' 
#' @param pheno.v
#' A phenotype integer vector for the cells, with normal cells indexed as 0, and
#' preneoplastic cells indexed with 1. Cells from more advanced stages including
#' cancer are indexed with positive integers larger than 1, with higher the integer
#' the more advanced the disease stage is. Note that this vector must have at least
#' 2 disease stages (0,1) representing normal and preneoplastic cells.
#'
#' @param thDTFA
#' The significance threshold for defining differentially active TFs
#' between normal and preneoplastic cells. By default this is
#' Bonferroni-adjusted P<0.05, but can be numeric, in which case the
#' provided value is used.
#' 
#' @param thHits
#' The significance threshold on P-values for declaring activated TFs
#' for each cell. By default this value is 0.05.
#' 
#' @return A list with the following elements
#'
#' @return tfal
#' The transcription factor activation load for all cells. In practice we are
#' interested in those of the preneoplastic cells only, but for
#' convenience we return the load for all cells.
#'
#' @return tfaf
#' As the TFAL, but expressed as a fraction of the total number
#' of tissue-specific TFs
#'
#' @return iTF
#' A vector specifying the gene-identifiers of the TFs that
#' are activated.
#'
#' @return hitMatrix
#' A binary 0-1 matrix defined over the activated TFs and all cells, with 1
#' indicating that the given TF is significantly activated in the given cell.
#'
#' @return tfa
#' The TFA-matrix input object.
#' 
#' @references 
#' 
#' 
#' 
#' @examples 
#' 
#'  
#' @export
#'     

CompTFAL <- function(tfa.m,pheno.v,thDTFA="BF",thHits=0.05){
   nTF <- nrow(tfa.m);
   n.idx <- which(pheno.v==0);
   preN.idx <- which(pheno.v==1);    
   statTFA.m <- tfa.m[,1:2];
   colnames(statTFA.m) <- c("t","P");
   for(r in 1:nrow(tfa.m)){
    lm.o <- lm(tfa.m[r,c(n.idx,preN.idx)] ~ pheno.v[c(n.idx,preN.idx)]);
    statTFA.m[r,] <- summary(lm.o)$coeff[2,3:4];
   }
   if(thDTFA=="BF"){
     sig.idx <- which(statTFA.m[,2] < 0.05/nTF);
   }
   else if (is.numeric(thDTFA)){
     sig.idx <- which(statTFA.m[,2] < thDTFA);
   }
   sigUP.idx <- intersect(sig.idx,which(statTFA.m[,1]>0));
   print(paste("The number of activated TFs=",length(sigUP.idx),sep=""));
   if(length(sigUP.idx)==0){
     print("No TFs are activated, so no TFAL can be computed!");
     tfal.v <- NULL; tfaf.v <- NULL; itf.v <- NULL;
     stop;
   }
   else {
    itf.v <- rownames(tfa.m)[sigUP.idx];
    hitActTF.m <- matrix(0,nrow=length(sigUP.idx),ncol=ncol(tfa.m));
    rownames(hitActTF.m) <- rownames(statTFA.m)[sigUP.idx];
    for(tf in 1:nrow(hitActTF.m)){
     mu <- mean(tfa.m[sigUP.idx[tf],n.idx]);
     sd <- sd(tfa.m[sigUP.idx[tf],n.idx]);
     prob.v <- pnorm(tfa.m[sigUP.idx[tf],],mu,sd,lower.tail=FALSE);
     hitActTF.m[tf,which(prob.v < thHits)] <- 1;
     tfal.v <- apply(hitActTF.m,2,sum);
     tfaf.v <- tfal.v/nTF;
    }
   }

   return(list(tfal=tfal.v,tfaf=tfaf.v,iTF=itf.v,hitMatrix=hitActTF.m,tfa=tfa.m));   
    
} ## EOF





