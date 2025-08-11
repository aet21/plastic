#' @title 
#' Computation of Cancer Risk score
#' 
#' @aliases CompCancerRisk
#'  
#' @description 
#' This function takes the output from CompTFIL and estimates a cancer-risk score
#' for each preneoplastic cell.
#' 
#' @param tfa.m
#' The TFA-matrix defined over tissue-specific TFs (rows) and cells (columns)
#' 
#' @param pheno.v
#' A phenotype integer vector for the cells, with normal cells indexed as 0,
#' preneoplastic cells indexed with 1, and with all other disease stages
#' indexed with increasing positive integers larger than 1, i.e. with higher
#' integers the more advanced the disease stage is. Invasive cancer is assumed to be
#' the largest integer. Note that in order to compute a Cancer Risk Score, it is
#' necessary to include cells from an invasive cancer state.
#' 
#' @param thDTFA
#' The significance threshold for defining differentially active TFs.
#' By default this is Bonferroni-adjusted P<0.05, but can be numeric
#' in which case the provided value is used.
#'
#' @param option
#' Character specifying if cancer risk score is to be computed using the simple correlation option,
#' or whether the diffusion based version of the cancer risk score should also be computed.
#' 
#' 
#' @return A list with the following elements
#'
#' @return risk
#' The cancer-risk score of all preneoplastic cells, computed by correlation
#' similarity in TFA-space to cancer cells, using only TFs that are significantly
#' inactivated in preneoplastic cells compared to normal cells.
#'
#' @return riskTip
#' The cancer-risk score of all preneoplastic cells, computed by correlation
#' similarity in TFA-space to the cancer and cancer-free tip points, as inferred
#' using Diffusion Maps on the full TFA-matrix.
#'
#' @return dmap
#' The diffusion map object returned by applying Diffusion Maps to the TFA-matrix.
#'
#' @return dc
#' The diffusion component matrix with the number of components estimated
#' via RMT.
#'
#' @return transM
#' The Markov transition matrix.
#'
#' @return root
#' The index of the root cell.
#'
#' @return tipC
#' The index of the cancer-fate tip point, as inferred using Diffusion Maps.
#'
#' @return tipN
#' The index of the cancer-free fate tip point, as inferred using Diffusion Maps.
#'
#'
#' 
#' @references 
#' 
#' 
#' 
#' @examples 
#' 
#' @import destiny
#' @import igraph
#' @import irlba
#' 
#' @export
#'     

CompCancerRisk <- function(tfa.m,pheno.v,thDTFA="BF",option=c("simple","both")){
   sd.v <- apply(tfa.m,1,sd);
   tmp.m <- (tfa.m - rowMeans(tfa.m))/sd.v ;
   preN.idx <- which(pheno.v==1);    
   n.idx <- which(pheno.v==0);
   c.idx <- which(pheno.v==max(pheno.v));
   riskTip.v <- NULL; dc.m <- NULL; tipN.idx <- NULL; tipC.idx <- NULL; root.idx <- NULL; dmap.o <- NULL; transM.m <- NULL; dpt.o <- NULL;
   if(option=="both"){
   rmt.o <- FastEstDimRMT(tmp.m);
   ndim <- rmt.o$dim;
   print("Running Diffusion Maps");
   dmap.o <- DiffusionMap(data=t(tfa.m),k=ndim,verbose=TRUE);
   dc.m <- eigenvectors(dmap.o);
   transM.m <- dmap.o@transitions;
   pN.m <- transM.m[n.idx,n.idx];
   colnames(pN.m) <- paste("Cell", 1:nrow(pN.m), sep = "")
   rownames(pN.m) <- paste("Cell", 1:nrow(pN.m), sep = "")
   gr.o <- graph_from_adjacency_matrix(pN.m, mode = "undirected", weighted = TRUE, diag = TRUE);
   nsteps <- floor(0.25 * nrow(pN.m))
   print(paste("Finding Root State: Running walk-trap community detection with max.step number= ",  nsteps, sep = ""));
   walk.o <- cluster_walktrap(gr.o, steps = nsteps)
   distr.v <- summary(factor(walk.o$member))
   max.idx <- which.max(distr.v)
   maxclID <- max.idx
   selcand.idx <- n.idx[which(walk.o$member == max.idx)]
   print(paste("Walk-trap identified a root-state consisting of ", length(selcand.idx), " cells", sep = ""))
   pos.m <- dc.m[selcand.idx, ]
   med.m <- apply(pos.m, 2, median)
   mad.v <- apply(abs(pos.m - med.m), 1, mean)
   root.idx <- selcand.idx[which.min(mad.v)]
   print("Estimating Diffusion Pseudotime");
   dpt.o <- DPT(dmap.o,tips=root.idx);
   tips.idx <- which(dpt.o@tips[,1]);
   print(paste("Number of tips=",length(tips.idx),sep=""));
   tipsNC.idx <- setdiff(tips.idx,root.idx);
   pccC.m <- cor(tfa.m[,tipsNC.idx],tfa.m[,c.idx]);
   pccC.v <- rowMeans(pccC.m);
   pccN.m <- cor(tfa.m[,tipsNC.idx],tfa.m[,n.idx]);
   pccN.v <- rowMeans(pccN.m);
   tipC.idx <- tipsNC.idx[which(pccC.v>pccN.v)];
   tipN.idx <- tipsNC.idx[which(pccN.v>pccC.v)];    
   if( (length(tipC.idx)==1) && (length(tipN.idx)==1) ){
      riskTip.v <- 0.5*(cor(tfa.m[,preN.idx],tfa.m[,tipC.idx]) - cor(tfa.m[,preN.idx],tfa.m[,tipN.idx]));
   }
   else {
       print("Not identified unambiguous cancer and cancer-free tip points");
       riskTip.v <- NULL;
   }
   }

   nTF <- nrow(tfa.m);
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
   sigDN.idx <- intersect(sig.idx,which(statTFA.m[,1]<0));
   print(paste("The number of inactivated TFs=",length(sigDN.idx),sep=""));
   if(length(sigDN.idx)==0){
     print("No TFs are inactivated in preneoplastic cells, so no cancer-risk score can be computed");
     risk.v <- NULL;
     stop;
   }
   else {
     risk.v <- rowMeans(cor(tfa.m[sigDN.idx,preN.idx],tfa.m[sigDN.idx,c.idx]));
   }   
   return(list(risk=risk.v,riskTip=riskTip.v,dmap=dmap.o,dc=dc.m,transM=transM.m,dpt=dpt.o,root=root.idx,tipC=tipC.idx,tipN=tipN.idx));   
    
} ## EOF

#### Auxilliary functions

FastEstDimRMT <- function(data.m,nEigEst=NULL,plot = TRUE){

    print("Centering and scaling matrix");
    M <- apply(data.m, 2, function(X) {
        (X - mean(X))/sd(X)
    })
    print("Done, now perform SVD");
    sigma2 <- var(as.vector(M))
    Q <- nrow(data.m)/ncol(data.m)
    maxNEig <- min(dim(M));
    if(is.null(nEigEst)){
      nEigEst <- floor(maxNEig/10);
    }
    lambdaMAX <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
    lambdaMIN <- sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))

    if(min(dim(data.m))>=500){    
     print(paste("Using Fast IRLBA to approximate ",nEigEst," top singular values",sep=""));
     workDim <- 2*nEigEst;
     svd.o <- irlba(M, nv = nEigEst, nu = nEigEst, maxit = 1000, work = workDim)
     evals <- (svd.o$d^2)/nrow(M);
     print("Done");
    }
    else {
      print("Performing full SVD since dimensionality of data matrix is not big");
      nEigEst <- maxNEig;
      svd.o <- svd(M);
      evals <- (svd.o$d^2)/nrow(M);
      print("Done");
    }
    estdens.o <- density(evals,from = min(evals),to = max(evals), cut = 0)
    intdim <- length(which(evals > lambdaMAX));

    if (plot) {
        step <- (lambdaMAX-lambdaMIN)/10000;
        lambda.v <- seq(lambdaMIN+10^{-6}, lambdaMAX, by = step)
        dens.v <- (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - lambda.v) * (lambda.v - lambdaMIN))/lambda.v;
        thdens.o <- list(min = lambdaMIN, max = lambdaMAX, step = step, lambda = lambda.v, dens = dens.v);
        minx <- min(min(thdens.o$lambda), min(evals))
        maxx <- max(max(thdens.o$lambda), max(evals))
        miny <- min(min(thdens.o$dens), min(estdens.o$y))
        maxy <- max(max(thdens.o$dens), max(estdens.o$y))
        pdf("RMTplot.pdf", width = 4, height = 4)
        plot(thdens.o$lambda, thdens.o$dens, xlim = c(0.5, maxx), 
            ylim = c(miny, maxy), type = "l", col = "green", 
            xlab = "Folded Eigenvalues", ylab = "density", lwd = 1)
        legend(x=0.5*(lambdaMAX-lambdaMIN),y=0,legend=c("NULL","OBS"),lwd=2,col=c("green","red"),horiz=TRUE,cex=0.5);
        i <- min(which(estdens.o$x > min(evals)))
        f <- max(which(estdens.o$x < max(evals)))
        points(x = estdens.o$x[i:f], y = estdens.o$y[i:f], type = "l", 
               col = "red", lwd=1);
        if(intdim>0){
         for (i in 1:intdim) {
            abline(v = evals[i], col = "red", lwd = 1)
         }
        }
        dev.off()
    }
    return(list(dim = intdim,nEigEst=nEigEst,evals = evals,thmin=lambdaMIN,thmax=lambdaMAX));
}





