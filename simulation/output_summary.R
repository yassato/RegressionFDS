#output nei_lmm()

#base stat
library(vcfR)
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
#Figure S1: Structure of simualted genomes

svg("FigS1_GenomeSummary.svg",width=6,height=4)
par(mfrow=c(3,4))
par(mai=c(0.15, 0.4, 0.15, 0.1))

resNFDS = readRDS("output/popStatNFDSv2_self.rds")
resPFDS = readRDS("output/popStatPFDS_self.rds")
resOD = readRDS("output/popStatOD_self.rds")
resSTVS = readRDS("output/popStatSTVS_self.rds")

boxplot(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp,las=1,col="grey",lwd=1.5,ylim=c(0,400))
boxplot(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst,las=1,col="grey",lwd=1.5,ylim=c(0,1))

resNFDS = readRDS("output/popStatNFDSv2.rds")
resPFDS = readRDS("output/popStatPFDS.rds")
resOD = readRDS("output/popStatOD.rds")
resSTVS = readRDS("output/popStatSTVS.rds")

boxplot(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp,las=1,col="grey",lwd=1.5,ylim=c(0,50))
boxplot(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst,las=1,col="grey",lwd=1.5,ylim=c(0,1))

resNFDS = readRDS("output/popStatNFDSv2_all.rds")
resPFDS = readRDS("output/popStatPFDS_all.rds")
resOD = readRDS("output/popStatOD_all.rds")
resSTVS = readRDS("output/popStatSTVS_all.rds")

boxplot(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp,las=1,col="grey",lwd=1.5,ylim=c(0,4500))
boxplot(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz,las=1,col="grey",lwd=1.5,ylim=c(0,0.5))
boxplot(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst,las=1,col="grey",lwd=1.5,ylim=c(0,1))

dev.off()

###########################
# Figure 3: LMM of beta_2
svg("Fig3_beta2LMM.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### plant
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantv2.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantV2.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)


## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantv2.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

#### animal
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)

## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

dev.off()

########################
# Figure 4: LM of beta2

svg("Fig4_beta2LM.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### plant
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantLM.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plantLM.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01plantLM.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01plantLM.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantLM.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plantLM.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01plantLM.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01plantLM.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)


## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

#### animal
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalLM.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animalLM.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01animalLM.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01animalLM.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalLM.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animalLM.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01animalLM.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01animalLM.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)

## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

dev.off()


######################
# Figure S2: LMM of beta_1

svg("FigS2_beta1LMM.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### plant
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantv2.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantV2.rds"))
  beta1 = c(beta1,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  beta2 = c(beta2,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  beta3 = c(beta3,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  beta4 = c(beta4,results$gwas_out$beta_self[results$ans[,1]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)


## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantv2.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

#### animal
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  beta1 = c(beta1,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  beta2 = c(beta2,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  beta3 = c(beta3,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  beta4 = c(beta4,results$gwas_out$beta_self[results$ans[,1]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)

## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

dev.off()


#############
# Figure S3: LM of beta_1

svg("FigS3_beta1LM.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### plant
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantLM.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plantLM.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01plantLM.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01plantLM.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantLM.rds"))
  beta1 = c(beta1,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plantLM.rds"))
  beta2 = c(beta2,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01plantLM.rds"))
  beta3 = c(beta3,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01plantLM.rds"))
  beta4 = c(beta4,results$gwas_out$beta_self[results$ans[,1]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)


## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01plantLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

#### animal
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalLM.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animalLM.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01animalLM.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01animalLM.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalLM.rds"))
  beta1 = c(beta1,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animalLM.rds"))
  beta2 = c(beta2,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01animalLM.rds"))
  beta3 = c(beta3,results$gwas_out$beta_self[results$ans[,1]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01animalLM.rds"))
  beta4 = c(beta4,results$gwas_out$beta_self[results$ans[,1]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(0.1,0.1,0.1,0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)

## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_s2_01animalLM.rds"))
  p = c(p, -log10(results$gwas_out$p_self))
  ans = c(ans, results$ans[,1])
}
plot(roc(ans~p),add=T,col="grey")

dev.off()


###########################
# Figure S4: LMM of beta_2 with dominance encoding
svg("FigS4_beta2LMMdomi.svg",width=6,height=4)
par(mfrow=c(2,3))
par(mai=c(0.15, 0.4, 0.15, 0.1))

#### plant
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDSdomi_plantv2.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDSdomi_plant.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_ODdomi_plant.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVSdomi_plant.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDSdomi_plantv2.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDSdomi_plant.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_ODdomi_plant.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVSdomi_plant.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)


## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDSdomi_plantv2.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDSdomi_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_ODdomi_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVSdomi_plant.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

#### animal
## auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDSdomi_animalv2.rds"))
  auc1 = c(auc1,results$AUC[2])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDSdomi_animal.rds"))
  auc2 = c(auc2,results$AUC[2])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_ODdomi_animal.rds"))
  auc3 = c(auc3,results$AUC[2])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVSdomi_animal.rds"))
  auc4 = c(auc4,results$AUC[2])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


## beta
beta1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDSdomi_animalv2.rds"))
  beta1 = c(beta1,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDSdomi_animal.rds"))
  beta2 = c(beta2,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_ODdomi_animal.rds"))
  beta3 = c(beta3,results$gwas_out$beta_nei[results$ans[,2]==1])
}

beta4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVSdomi_animal.rds"))
  beta4 = c(beta4,results$gwas_out$beta_nei[results$ans[,2]==1])
}

boxplot(beta1,beta2,beta3,beta4,las=1,col="grey",lwd=1.5)
points(c(1,2,3,4,4),c(-0.1,0.1,0.1,-0.1,0.1),col="red",pch=20,cex=2)
abline(h=0.0,lty=2)

## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDSdomi_animalv2.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),ylim=c(0,1),xlim=c(1.0,0.0),col="blue",las=1,cex.axis=1,ylab="",xlab="",mar=c(1,2,1,1)+.1)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDSdomi_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="red")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_ODdomi_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="black")

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVSdomi_animal.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
plot(roc(ans~p),add=T,col="grey")

dev.off()





####old#####
#### plant
#auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01plantv2.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_plant.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)

#### animal
#auc
auc1=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
  auc1 = c(auc1,results$AUC[1])
}

auc2=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
  auc2 = c(auc2,results$AUC[1])
}

auc3=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_animal.rds"))
  auc3 = c(auc3,results$AUC[1])
}

auc4=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
  auc4 = c(auc4,results$AUC[1])
}

boxplot(auc1,auc2,auc3,auc4,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)


d1 = read.table("clipboard",header=TRUE)
boxplot(d1$NFDS,d1$PFDS,d1$OD,d1$STVS,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)

d2 = read.table("clipboard",header=TRUE)
boxplot(d2$NFDS,d2$PFDS,d2$OD,d2$STVS,ylim=c(0.4,1.0),las=1,col="grey",lwd=1.5)
abline(h=0.5,lty=2)

#ROC figure
par(mfrow=c(2,2))

library(pROC)

i=1
results = readRDS(paste0("output/",i,"_NFDS_s2_01plantv2.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),ylim=c(0,1),xlim=c(1,0),col="blue",las=1,cex.axis=1.25)

results = readRDS(paste0("output/",i,"_PFDS_s2_01plant.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),add=T,col="red")

results = readRDS(paste0("output/",i,"_OD_plant.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),add=T,col="black")

results = readRDS(paste0("output/",i,"_STVS_plant.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),add=T,col="grey")

i=1
results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),ylim=c(0,1),xlim=c(1,0),col="blue",las=1,cex.axis=1.25)

results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),add=T,col="red")

results = readRDS(paste0("output/",i,"_OD_animal.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),add=T,col="black")

results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
p = -log10(results$gwas_out$p_self)
ans = results$ans[,1]
plot(roc(ans~p),add=T,col="grey")


i=1
results = readRDS(paste0("output/",i,"_NFDS_s2_01animalv2.rds"))
p = -log10(results$gwas_out$p_nei)
ans = results$ans[,2]
plot(roc(ans~p),ylim=c(0,1),xlim=c(1,0),col="blue",las=1,cex.axis=1.25)

results = readRDS(paste0("output/",i,"_PFDS_s2_01animal.rds"))
p = -log10(results$gwas_out$p_nei)
ans = results$ans[,2]
plot(roc(ans~p),add=T,col="red")

results = readRDS(paste0("output/",i,"_OD_animal.rds"))
p = -log10(results$gwas_out$p_nei)
ans = results$ans[,2]
plot(roc(ans~p),add=T,col="black")

results = readRDS(paste0("output/",i,"_STVS_animal.rds"))
p = -log10(results$gwas_out$p_nei)
ans = results$ans[,2]
plot(roc(ans~p),add=T,col="grey")

