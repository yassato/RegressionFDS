#Figures from output nei_lmm()

#base stat
library(vcfR)
library(pROC)
popStat = function(type=c("NFDS","PFDS","OD","STVS"), iter=30) {
  vec = c()
  for(i in 1:iter) {
    f_name = paste0(i,"_",type)
    f_path = paste0("output/",f_name,".vcf.gz")
    
    d0 = read.vcf(f_path)
    d = select.snps(d0, d0@snps$maf>0.01)
    
    ann = read.vcfR(f_path)
    ann = ann[which(d0@snps$maf>0.01),]
    
    pop = factor(rep(1:10,each=200))
    
    # #all
    # res = genetic_diff(vcf=ann, pops=pop, method="nei")
    # n_snp = ncol(d)
    # meanMAF = mean(d@snps$maf)
    # meanHz = mean(d@snps$hz)
    # meanGst = mean(res$Gst)
    # #replacing "# causal" by "# all" give a summary of the genome-wide structure 
    
    #causal
    res = genetic_diff(vcf=ann[grep("MT=3",ann@fix[,"INFO"])], pops=pop, method="nei")
    n_snp = length(grep("MT=3",ann@fix[,"INFO"]))
    meanMAF = mean(d@snps$maf[grep("MT=3",ann@fix[,"INFO"])])
    meanHz = mean(d@snps$hz[grep("MT=3",ann@fix[,"INFO"])])
    meanGst = mean(res$Gst)
    
    vec = rbind(vec,
            c(n_snp,meanMAF,meanHz,meanGst)
            )
  }
  vec = as.data.frame(vec)
  colnames(vec) = c("snp","MAF","Hz","Gst")
  return(vec)
}

resNFDS = popStat("NFDS",30)
resPFDS = popStat("PFDS",30)
resOD = popStat("OD",30)
resSTVS_all = popStat("STVS",30)
saveRDS(resPFDS, file="output/popStatSTVS_self.rds")

#######################
#Supp Figure: Structure of simualted genomes

svg("GenomeSummary.svg",width=6,height=4)
par(mfrow=c(3,4))
par(mai=c(0.15, 0.4, 0.15, 0.1))

resNFDS = readRDS("output/popStatNFDS_self.rds")
resPFDS = readRDS("output/popStatPFDS_self.rds")
resOD = readRDS("output/popStatOD_self.rds")
resSTVS = readRDS("output/popStatSTVS_self.rds")

boxplot(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp,las=1,col="grey",lwd=1.5,ylim=c(0,400))
boxplot(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst,las=1,col="grey",lwd=1.5,ylim=c(0,1))

resNFDS = readRDS("output/popStatNFDS.rds")
resPFDS = readRDS("output/popStatPFDS.rds")
resOD = readRDS("output/popStatOD.rds")
resSTVS = readRDS("output/popStatSTVS.rds")

boxplot(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp,las=1,col="grey",lwd=1.5,ylim=c(0,50))
boxplot(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst,las=1,col="grey",lwd=1.5,ylim=c(0,1))

resNFDS = readRDS("output/popStatNFDS_all.rds")
resPFDS = readRDS("output/popStatPFDS_all.rds")
resOD = readRDS("output/popStatOD_all.rds")
resSTVS = readRDS("output/popStatSTVS_all.rds")

boxplot(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp,las=1,col="grey",lwd=1.5,ylim=c(0,4500))
boxplot(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst,las=1,col="grey",lwd=1.5,ylim=c(0,1))

dev.off()

###########################
# Main Figure: LMM of beta_2
svg("beta2LMMdomi.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### continuous
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)


#### split
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)

dev.off()

########################
# Supp Figure: LM of beta2

svg("beta2LMdomi.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### continuous
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)


#### split
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)

dev.off()


######################
# Supp Figure: LMM of beta_1

svg("beta1LMMdomi.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### continuous

## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5,ylim=c(-0.05,0.15))
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)


#### split
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5,ylim=c(-0.05,0.15))
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)

dev.off()


#############
# Supp Figure: LM of beta_1

svg("beta1LMdomi.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### continuous
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5,ylim=c(-0.05,0.15))
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)


#### split
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col=grey(0.25,0.75))

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  beta1 = c(beta1,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  beta2 = c(beta2,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  beta3 = c(beta3,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  beta4 = c(beta4,median(results$gwas_out$beta_self[results$ans[,1]==1]))
}

boxplot(beta1,beta2,beta3,beta4,las=1,col=c("red","skyblue",grey(0.25,0.75),"grey"),lwd=1.5,ylim=c(-0.05,0.15))
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),pch=1,cex=2)
abline(h=0.0,lty=2)

dev.off()


