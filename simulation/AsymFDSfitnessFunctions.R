# Figure presentation for fitness functions

##########
#Main Figure: dominant case

svg("AsymFDSdomi.svg", width=6, height=8)

par(mfrow=c(3,2))
par(mai=c(0.55, 0.55, 0.20, 0.20))

W_AA = function(x) (b12+b2)*(4*x-2*x^2-1)+1
W_Aa = function(x) (b12+b2)*(4*x-2*x^2-1)+1
W_aa = function(x) (b12-b2)*(4*x-2*x^2-1)+1

W_A = function(x) x*W_AA(x)+(1-x)*W_Aa(x)
W_a = function(x) x*W_Aa(x)+(1-x)*W_aa(x)

b2 = -0.2; b12 = 0.0
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.293,0.293*W_A(0.293)+(1-0.293)*W_a(0.293),pch=16,cex=1.5)

b2 = 0.2; b12 = 0.0
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.293,0.293*W_A(0.293)+(1-0.293)*W_a(0.293),pch=1,cex=1.5)

b2 = -0.2; b12 = -0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.293,0.293*W_A(0.293)+(1-0.293)*W_a(0.293),pch=16,cex=1.5)

b2 = 0.2; b12 = -0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.293,0.293*W_A(0.293)+(1-0.293)*W_a(0.293),pch=1,cex=1.5)

b2 = -0.2; b12 = 0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.293,0.293*W_A(0.293)+(1-0.293)*W_a(0.293),pch=16,cex=1.5)

b2 = 0.2; b12 = 0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.293,0.293*W_A(0.293)+(1-0.293)*W_a(0.293),pch=1,cex=1.5)

dev.off()


############
# Supp Figure: additive case

svg("AsymFDSadd.svg", width=6, height=8)

par(mfrow=c(3,2))
par(mai=c(0.55, 0.55, 0.20, 0.20))

W_AA = function(x) (b12+b2)*(2*x-1)+1
W_Aa = function(x) 1
W_aa = function(x) (b12-b2)*(2*x-1)+1

W_A = function(x) x*W_AA(x)+(1-x)*W_Aa(x)
W_a = function(x) x*W_Aa(x)+(1-x)*W_aa(x)

b2 = -0.2; b12 = 0
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.5,0.5*W_A(0.5)+(1-0.5)*W_a(0.5),pch=16,cex=1.5)

b2 = 0.2; b12 = 0
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.5,0.5*W_A(0.5)+(1-0.5)*W_a(0.5),pch=1,cex=1.5)

b2 = -0.2; b12 = -0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.5,0.5*W_A(0.5)+(1-0.5)*W_a(0.5),pch=16,cex=1.5)
p2 = 0.5*(b12-b2)/b12
points(p2,p2*W_A(p2)+(1-p2)*W_a(p2),pch=1,cex=1.5)

b2 = 0.2; b12 = -0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.5,0.5*W_A(0.5)+(1-0.5)*W_a(0.5),pch=1,cex=1.5)
p2 = 0.5*(b12-b2)/b12
points(p2,p2*W_A(p2)+(1-p2)*W_a(p2),pch=16,cex=1.5)

b2 = -0.2; b12 = 0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.5,0.5*W_A(0.5)+(1-0.5)*W_a(0.5),pch=16,cex=1.5)
p2 = 0.5*(b12-b2)/b12
points(p2,p2*W_A(p2)+(1-p2)*W_a(p2),pch=1,cex=1.5)

b2 = 0.2; b12 = 0.3
curve(W_A,ylim=c(0.5,1.5),las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve(W_a, add=T,lty=1,col=grey(0.0,0.5),lwd=1.5)
curve(x*W_A(x)+(1-x)*W_a(x),add=T,lty=2)
points(0.5,0.5*W_A(0.5)+(1-0.5)*W_a(0.5),pch=1,cex=1.5)
p2 = 0.5*(b12-b2)/b12
points(p2,p2*W_A(p2)+(1-p2)*W_a(p2),pch=16,cex=1.5)

dev.off()

###########
# Supp Figure: inbred case

svg("AsymFDSinbred.svg", width=6, height=8)

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

