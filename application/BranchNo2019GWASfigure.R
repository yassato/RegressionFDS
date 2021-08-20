#LMM
png("Fig4_ManhattanLMM.png", width=7, height=4, res=600, units="in")
layout(matrix(c(1,1,4,2,2,4,3,3,4),3,3,byrow=TRUE))
par(mai=c(0.4, 0.3, 0.25, 0.20))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_self)
colnames(res) = c("chr","pos","p")
gaston::manhattan(res,ylim=c(0,8),las=1)
abline(h=-log10(p_adj),lty=2,col=grey(0.5,0.5))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_nei)
colnames(res) = c("chr","pos","p")
gaston::manhattan(res,ylim=c(0,8),las=1)
abline(h=-log10(p_adj),lty=2,col=grey(0.5,0.5))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMMasym.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_sxn)
colnames(res) = c("chr","pos","p")
gaston::manhattan(res,ylim=c(0,8.5),las=1)
abline(h=-log10(p_adj),lty=2,col=grey(0.5,0.5))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))
hist(subset(gwas_out,P_nei<0.0001)$beta_nei,las=1,main="",xlab="",ylab="")

dev.off()



#QQ-plot LMM
png("FigS_QQ_LMM.png", width=8, height=3, res=600, units="in")
par(mfcol=c(1,3))
par(mai=c(0.6, 0.6, 0.25, 0.20))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_self)
colnames(res) = c("chr","pos","p")
gaston::qqplot.pvalues(res,las=1,main="")

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_nei)
colnames(res) = c("chr","pos","p")
gaston::qqplot.pvalues(res,las=1,main="")

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMMasym.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_sxn)
colnames(res) = c("chr","pos","p")
gaston::qqplot.pvalues(res,las=1,main="")

dev.off()



#LM
png("Fig5_ManhattanLM.png", width=7, height=4, res=600, units="in")
layout(matrix(c(1,1,4,2,2,4,3,3,4),3,3,byrow=TRUE))
par(mai=c(0.4, 0.3, 0.25, 0.20))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_self)
colnames(res) = c("chr","pos","p")
gaston::manhattan(res,las=1)
abline(h=-log10(p_adj),lty=2,col=grey(0.5,0.5))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_nei)
colnames(res) = c("chr","pos","p")
gaston::manhattan(res,ylim=c(0,8),las=1)
abline(h=-log10(p_adj),lty=2,col=grey(0.5,0.5))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMasym.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_sxn)
colnames(res) = c("chr","pos","p")
gaston::manhattan(res,ylim=c(0,8.5),las=1)
abline(h=-log10(p_adj),lty=2,col=grey(0.5,0.5))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))
hist(subset(gwas_out,P_nei<0.0001)$beta_nei,las=1,main="",xlab="",ylab="")

dev.off()



#QQ-plot LM
png("FigS_QQ_LM.png", width=8, height=3, res=600, units="in")
par(mfcol=c(1,3))
par(mai=c(0.6, 0.6, 0.25, 0.20))

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_self)
colnames(res) = c("chr","pos","p")
gaston::qqplot.pvalues(res,las=1,main="")

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_nei)
colnames(res) = c("chr","pos","p")
gaston::qqplot.pvalues(res,las=1,main="")

gwas_out = read.csv("210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMasym.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))

res = data.frame(gwas_out$Chr,gwas_out$Position,gwas_out$P_sxn)
colnames(res) = c("chr","pos","p")
gaston::qqplot.pvalues(res,las=1,main="")

dev.off()

