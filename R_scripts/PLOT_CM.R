###Circular manhatten###
CMplot(HOFR_original, plot.type="c", chr.labels=paste(c(1:29), sep=""), r=1.5, cir.legend=TRUE,
       outward=FALSE, ,amplify=TRUE, threshold = c(0.05, 0.1), threshold.lty=c(1,2),threshold.col=c("red",
       "blue"), signal.line=1, signal.col=c("red","blue"), signal.cex = 2, cir.legend.col="black", memo = "HOFR_original", 
	   cir.chr.h=1 ,chr.den.col="grey", file="jpg", dpi=300, ylim = 4)
	   
CMplot(HOFR_new, plot.type="c", chr.labels=paste(c(1:29), sep=""), r=1.5, cir.legend=TRUE,
       outward=FALSE, ,amplify=TRUE, threshold = c(0.05, 0.1), threshold.lty=c(1,2),threshold.col=c("red",
       "blue"), signal.line=1, signal.col=c("red","blue"), signal.cex = 2, cir.legend.col="black", memo = "HOFR_new", 
	   cir.chr.h=1 ,chr.den.col="grey", file="jpg", dpi=300, ylim = 4)
	   
CMplot(CH_original, plot.type="c", chr.labels=paste(c(1:29), sep=""), r=1.5, cir.legend=TRUE,
       outward=FALSE, ,amplify=TRUE, threshold = c(0.05, 0.1), threshold.lty=c(1,2),threshold.col=c("red",
       "blue"), signal.line=1, signal.col=c("red","blue"), signal.cex = 2, cir.legend.col="black", memo = "CH_original", 
	   cir.chr.h=1 ,chr.den.col="grey", file="jpg", dpi=300, ylim = 4)
	   
CMplot(CH_new, plot.type="c", chr.labels=paste(c(1:29), sep=""), r=1.5, cir.legend=TRUE,
       outward=FALSE, ,amplify=TRUE, threshold = c(0.05, 0.1), threshold.lty=c(1,2),threshold.col=c("red",
       "blue"), signal.line=1, signal.col=c("red","blue"), signal.cex = 2, cir.legend.col="black", memo = "CH_new", 
	   cir.chr.h=1 ,chr.den.col="grey", file="jpg", dpi=300, ylim = 4)
	   
CMplot(LM_original, plot.type="c", chr.labels=paste(c(1:29), sep=""), r=1.5, cir.legend=TRUE,
       outward=FALSE, ,amplify=TRUE, threshold = c(0.05, 0.1), threshold.lty=c(1,2),threshold.col=c("red",
       "blue"), signal.line=1, signal.col=c("red","blue"), signal.cex = 2, cir.legend.col="black", memo = "LM_original", 
	   cir.chr.h=1 ,chr.den.col="grey", file="jpg", dpi=300, ylim = 4)
	   
CMplot(LM_new, plot.type="c", chr.labels=paste(c(1:29), sep=""), r=1.5, cir.legend=TRUE,
       outward=FALSE, ,amplify=TRUE, threshold = c(0.05, 0.1), threshold.lty=c(1,2),threshold.col=c("red",
       "blue"), signal.line=1, signal.col=c("red","blue"), signal.cex = 2, cir.legend.col="black", memo = "LM_new", 
	   cir.chr.h=1 ,chr.den.col="grey", file="jpg", dpi=300, ylim = 4)
	   
	   
#qqplot

CMplot(CM_HOFR_CH_LM, plot.type="q",col=c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),threshold=1e6,
        signal.pch=19,signal.cex=1.5,signal.col="red",conf.int.col="grey",box=FALSE,multracks=
        TRUE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
