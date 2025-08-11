#' @title 
#' Build Waddington Landscape
#' 
#' @aliases Build Waddington
#'  
#' @description 
#' This function takes the estimated diffusion map and cellular plasticity values,
#' and builds a corresponding Waddington landscape. Of note, this function may need
#' to be tweaked to alter the region to be plotted, if particular features of the landscape
#' want to be emphasized.
#'  
#' @param dcALL.m
#' The diffusion map matrix with rows labeling cell and columns labeling DCs.
#' 
#' @param epiCT.idx
#' An index vector specifying the cellular states, which should be of the same length as rows
#' of `dcALL.m`
#' 
#' @param ccat.v
#' The cellular plasticity estimates for the same cells as specified in the previous
#' arguments
#'
#' @param theta
#' Angle for viewing landscape (longitude)
#'
#' @param phi
#' Angle for viewing landscape (polar)
#'
#' @param selCells.idx
#' An index vector specifying the indices of the cells to be shown
#'
#' @return Generates a figure
#'
#' 
#' 
#' 
#' @examples 
#' 
#' @import destiny
#' @import marray
#' @import plot3D
#'
#' 
#' @export
#'     

BuildWaddington <- function(dcALL.m,epiCT.idx,ccat.v,theta,phi,selCells.idx,x.v=seq(-0.01,0.015,0.001),y.v=seq(-0.01,0.015,001)){

epiCT.idx <- epiCT.idx[selCells.idx];
ccat.v <- ccat.v[selCells.idx];
dcALL.m <- dcALL.m[selCells.idx,];
    
xMIN.v <- vector();
zMIN.v <- vector();
yMIN.v <- vector();
nCT <- length(unique(epiCT.idx));
for(ct in 1:nCT){
  ct.idx <- which(epiCT.idx==ct);
  tmp.v <- colMeans(dcALL.m[ct.idx,1:2]);
  xMIN.v[ct] <- tmp.v[1];
  yMIN.v[ct] <- tmp.v[2];
  zMIN.v[ct] <- mean(ccat.v[ct.idx]);
}

zMIN2.v <- (max(zMIN.v)-zMIN.v)/(max(zMIN.v)-min(zMIN.v));
zMIN2.v <- -zMIN2.v - 0.15;

sdMIN.v <- c(0.001,0.001,0.001,0.0005,0.0005,0.0005);

zMIN2.v[6] <- -0.3
z.m <- matrix(nrow=length(x.v),ncol=length(y.v));
for(x in 1:length(x.v)){
  for(y in 1:length(y.v)){
      z.m[x,y] <- mapF(x.v[x],y.v[y],xMIN.v,yMIN.v,zMIN2.v,sdMIN.v,zinf=0);
  }
}

par(mar=c(2,2,2,2));
par(mfrow=c(1,2));
color.v <- maPalette(low="brown",mid="red",high="white",k=10);
breaks.v <- c(-10^6,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,-0.05,-0.025,0.5);
persp3D(x=x.v,y=y.v,z=z.m,xlab="DC1",ylab="DC2",colkey=FALSE,col=color.v,lighting=FALSE,contour=TRUE,zlim=c(-1.12,0),theta=theta,phi=phi,border="black",resfac=1,breaks=breaks.v);

} ## EOF

#### Auxilliary functions

mapF1 <- function(x,y,m,xMIN.v,yMIN.v,zMIN.v,sdMIN.v){
    xMIN <- xMIN.v[m];
    yMIN <- yMIN.v[m];
    zMIN <- zMIN.v[m];
    sdMIN <- sdMIN.v[m];               
    out1 <- zMIN*exp(-((x-xMIN)^2 + (y-yMIN)^2)/(2*sdMIN^2));
    out2 <- prod(1 - exp(-((x-xMIN.v[-m])^2 + (y-yMIN.v[-m])^2)/(2*sdMIN.v[-m]^2)) );
    out <- out1*out2;
    return(out);
}

mapF <- function(x,y,xMIN.v,yMIN.v,zMIN.v,sdMIN.v,zinf){
    out <- zinf;
    nMIN <- length(xMIN.v);
    for(m in 1:nMIN){
       out <- out + mapF1(x,y,m,xMIN.v,yMIN.v,zMIN.v,sdMIN.v);
    }
    return(out);
}

