#' @title 
#' Diffusion Map Plot
#' 
#' @aliases PlotDMAP
#'  
#' @description 
#' This function takes the output of CompCancerRisk function and
#' generates a Diffusion Map plot, indicating root-state and the
#' cancer-free and cancer-fate tip points.
#' 
#' @param crs.o
#' The output of CompCancerRisk.
#'
#' @param type
#' The type of diffusion map desired. If DPT, then
#' cells are colored by diffusion pseudotime. If type
#' is Stage, then user needs to specify the disease
#' stage vector with stages labeled by integers. In this
#' option the user must also specify a color vector to
#' label each disease stage.
#'
#' @param stage
#' An integer vector specifying the stage in cancer progression.
#'
#' @param col
#' A color vector specifying the color for each unique disease stage.
#' 
#' @return out
#' A logical.
#'
#' @references 
#' 
#' 
#' 
#' @examples 
#' 
#' @import destiny
#' @import marray
#' @import plotrix
#' 
#' @export
#'     

PlotDMAP <- function(crs.o,type=c("DPT","Stage"),stage=NULL,col=NULL){
   root <- crs.o$root;
   tipC <- crs.o$tipC;
   tipN <- crs.o$tipN;    
   dc.m <- crs.o$dc;

   if(type=="DPT"){
   dpt.o <- crs.o$dpt;
   dpt.v <- dpt.o[[paste("DPT",root,sep="")]];
   q.v <- quantile(dpt.v,probs=seq(0,1,0.1));
   colorDPT.v <- rep(NA,length(dpt.v));
   colorDPT.idx <- rep(NA,length(dpt.v));
   for(qt in 1:10){
     colorDPT.idx[which(dpt.v >= q.v[qt])] <- qt;
   }    
   colorDPT.v <- maPalette(low="grey",high="darkgreen",mid="green",k=10);
   plot(dc.m[,1],dc.m[,2],col=colorDPT.v[colorDPT.idx],pch=23,bg=colorDPT.v[colorDPT.idx],axes=FALSE,cex=0.5,xlab="",ylab="");
   mtext("DC1",side=1,line=0.1);
   mtext("DC2",side=2,line=0.1);
   text(x=dc.m[root,1],y=dc.m[root,2],"Root",col="red",pos=4,font=2);
   points(dc.m[root,1],dc.m[root,2],col="red",cex=1.5,bg="red",pch=23);   
   text(x=dc.m[tipC,1],y=dc.m[tipC,2],"Cancer-fate",col="red",pos=3,font=2);
   points(dc.m[tipC,1],dc.m[tipC,2],col="red",cex=1.5,bg="red",pch=23);   
   text(x=dc.m[tipN,1],y=dc.m[tipN,2],"CancerFree-fate",col="red",pos=2,font=2);
   points(dc.m[tipN,1],dc.m[tipN,2],col="red",cex=1.5,bg="red",pch=23);   

   rangeX.v <- range(dc.m[,1]);
   rangeY.v <- range(dc.m[,2]);
   diffY <- rangeY.v[2] -rangeY.v[1];
   diffX <- rangeX.v[2] -rangeX.v[1];    
   par(xpd=TRUE);
     plotrix::color.legend(xl=rangeX.v[2]-diffX/10,yb=rangeY.v[2]-diffY/5,xr=rangeX.v[2],yt=rangeY.v[2],legend=c(0,round(median(dpt.v),2),ceiling(max(dpt.v))),rect.col=colorDPT.v,cex=0.75,align="lt",gradient="y")
     text(x=rangeX.v[2]-diffX/20,y=rangeY.v[2]-diffY/5,"DPT",font=2,cex=1,pos=1);
   par(xpd=FALSE);
   }
   else if (type=="Stage"){
   color.v <- col;
   if(min(stage)==0){ stage <- stage+1;}
   plot(dc.m[,1],dc.m[,2],col=color.v[stage],pch=23,bg=color.v[stage],axes=FALSE,cex=0.5,xlab="",ylab="");
   mtext("DC1",side=1,line=0.1);
   mtext("DC2",side=2,line=0.1);
   text(x=dc.m[root,1],y=dc.m[root,2],"Root",col="red",pos=4,font=2);
   points(dc.m[root,1],dc.m[root,2],col="red",cex=1.5,bg="red",pch=23);   
   text(x=dc.m[tipC,1],y=dc.m[tipC,2],"Cancer-fate",col="red",pos=3,font=2);
   points(dc.m[tipC,1],dc.m[tipC,2],col="red",cex=1.5,bg="red",pch=23);   
   text(x=dc.m[tipN,1],y=dc.m[tipN,2],"CancerFree-fate",col="red",pos=2,font=2);
   points(dc.m[tipN,1],dc.m[tipN,2],col="red",cex=1.5,bg="red",pch=23);   

   rangeX.v <- range(dc.m[,1]);
   rangeY.v <- range(dc.m[,2]);
   diffY <- rangeY.v[2] -rangeY.v[1];
   diffX <- rangeX.v[2] -rangeX.v[1];    
   par(xpd=TRUE);
     if(min(stage)==1){ stage <- stage-1;}
     plotrix::color.legend(xl=rangeX.v[2]-diffX/10,yb=rangeY.v[2]-diffY/5,xr=rangeX.v[2],yt=rangeY.v[2],legend=unique(stage),rect.col=color.v,cex=0.75,align="lt",gradient="y")
     text(x=rangeX.v[2]-diffX/20,y=rangeY.v[2]-diffY/5,"Stage",font=2,cex=1,pos=1);
   par(xpd=FALSE);
   }
       

   out <- TRUE;
   return(out);
    
} ## EOF








