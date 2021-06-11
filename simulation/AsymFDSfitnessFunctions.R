############
# additive #
############

#Figure 2
svglite::svglite("Fig2_AsymFDSadd.svg", width=6, height=8)

par(mfrow=c(3,2))
par(mai=c(0.55, 0.55, 0.20, 0.20))

b2 = -0.2; b12 = 0
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
abline(a=1,b=0,col=grey(0.0,0.66),lwd=1)
curve((b12-b2)*(2*x-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1)
curve((x^2)*((b12+b2)*(2*x-1))+((1-x)^2)*(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,1,pch=16,cex=1.5)

b2 = 0.2; b12 = 0
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
abline(a=1,b=0,col=grey(0.0,0.66),lwd=1)
curve((b12-b2)*(2*x-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1)
curve((x^2)*((b12+b2)*(2*x-1))+((1-x)^2)*(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,1,pch=1,cex=1.5)

b2 = -0.2; b12 = -0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
abline(a=1,b=0,col=grey(0.0,0.66),lwd=1)
curve((b12-b2)*(2*x-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1)
curve((x^2)*((b12+b2)*(2*x-1))+((1-x)^2)*(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=16,cex=1.5)

b2 = 0.2; b12 = -0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
abline(a=1,b=0,col=grey(0.0,0.66),lwd=1)
curve((b12-b2)*(2*x-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1)
curve((x^2)*((b12+b2)*(2*x-1))+((1-x)^2)*(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=1,cex=1.5)

b2 = -0.2; b12 = 0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
abline(a=1,b=0,col=grey(0.0,0.66),lwd=1)
curve((b12-b2)*(2*x-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1)
curve((x^2)*((b12+b2)*(2*x-1))+((1-x)^2)*(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=16,cex=1.5)

b2 = 0.2; b12 = 0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
abline(a=1,b=0,col=grey(0.0,0.66),lwd=1)
curve((b12-b2)*(2*x-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1)
curve((x^2)*((b12+b2)*(2*x-1))+((1-x)^2)*(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=1,cex=1.5)

dev.off()


##########
#dominant#
##########

#Figure S2
svglite::svglite("FigS2_AsymFDSdomi.svg", width=6, height=8)

par(mfrow=c(3,2))
par(mai=c(0.55, 0.55, 0.20, 0.20))
b2 = -0.2; b12 = 0.0
curve((b12+b2)*(4*x-2*x^2-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12+b2)*(4*x-2*x^2-1)+1,add=T,col=grey(0.0,0.66),lwd=1.5)
curve((b12-b2)*(4*x-2*x^2-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1.5)
curve((x^2)*((b12+b2)*(4*x-2*x^2-1))+2*x*(1-x)*((b12+b2)*(4*x-2*x^2-1))+((1-x)^2)*((b12-b2)*(4*x-2*x^2-1))+1,add=T,lty=2)
points(0.3,1,pch=16,cex=1.5)

b2 = 0.2; b12 = 0.0
curve((b12+b2)*(4*x-2*x^2-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12+b2)*(4*x-2*x^2-1)+1,add=T,col=grey(0.0,0.66),lwd=1.5)
curve((b12-b2)*(4*x-2*x^2-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1.5)
curve((x^2)*((b12+b2)*(4*x-2*x^2-1))+2*x*(1-x)*((b12+b2)*(4*x-2*x^2-1))+((1-x)^2)*((b12-b2)*(4*x-2*x^2-1))+1,add=T,lty=2)
points(0.3,1,pch=1,cex=1.5)

b2 = -0.2; b12 = -0.3
curve((b12+b2)*(4*x-2*x^2-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12+b2)*(4*x-2*x^2-1)+1,add=T,col=grey(0.0,0.66),lwd=1.5)
curve((b12-b2)*(4*x-2*x^2-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1.5)
curve((x^2)*((b12+b2)*(4*x-2*x^2-1))+2*x*(1-x)*((b12+b2)*(4*x-2*x^2-1))+((1-x)^2)*((b12-b2)*(4*x-2*x^2-1))+1,add=T,lty=2)
points(0.3,1,pch=16,cex=1.5)

b2 = 0.2; b12 = -0.3
curve((b12+b2)*(4*x-2*x^2-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12+b2)*(4*x-2*x^2-1)+1,add=T,col=grey(0.0,0.66),lwd=1.5)
curve((b12-b2)*(4*x-2*x^2-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1.5)
curve((x^2)*((b12+b2)*(4*x-2*x^2-1))+2*x*(1-x)*((b12+b2)*(4*x-2*x^2-1))+((1-x)^2)*((b12-b2)*(4*x-2*x^2-1))+1,add=T,lty=2)
points(0.3,1,pch=1,cex=1.5)

b2 = -0.2; b12 = 0.3
curve((b12+b2)*(4*x-2*x^2-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12+b2)*(4*x-2*x^2-1)+1,add=T,col=grey(0.0,0.66),lwd=1.5)
curve((b12-b2)*(4*x-2*x^2-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1.5)
curve((x^2)*((b12+b2)*(4*x-2*x^2-1))+2*x*(1-x)*((b12+b2)*(4*x-2*x^2-1))+((1-x)^2)*((b12-b2)*(4*x-2*x^2-1))+1,add=T,lty=2)
points(0.3,1,pch=16,cex=1.5)

b2 = 0.2; b12 = 0.3
curve((b12+b2)*(4*x-2*x^2-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12+b2)*(4*x-2*x^2-1)+1,add=T,col=grey(0.0,0.66),lwd=1.5)
curve((b12-b2)*(4*x-2*x^2-1)+1,add=T,lty=1,col=grey(0.0,0.33),lwd=1.5)
curve((x^2)*((b12+b2)*(4*x-2*x^2-1))+2*x*(1-x)*((b12+b2)*(4*x-2*x^2-1))+((1-x)^2)*((b12-b2)*(4*x-2*x^2-1))+1,add=T,lty=2)
points(0.3,1,pch=1,cex=1.5)

dev.off()


###########
# inbred# #
###########

#Figure S3
svglite::svglite("FigS3_AsymFDSinbred.svg", width=6, height=8)

par(mfrow=c(3,2))
par(mai=c(0.55, 0.55, 0.20, 0.20))

b2 = -0.2; b12 = 0
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,1,pch=16,cex=1.5)

b2 = 0.2; b12 = 0
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,1,pch=1,cex=1.5)

b2 = -0.2; b12 = -0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=16,cex=1.5)

b2 = 0.2; b12 = -0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=1,cex=1.5)

b2 = -0.2; b12 = 0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=16,cex=1.5)

b2 = 0.2; b12 = 0.3
curve((b12+b2)*(2*x-1)+1,ylim=c(0.5,1.5),
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1,add=T,lty=2)
points(0.5,(b12+b2)*(2*0.5-1)+1,pch=1,cex=1.5)

dev.off()
