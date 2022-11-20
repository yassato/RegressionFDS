#Figures from output nei_lmm()

#base stat
library(vcfR); library(pROC)
library(ggplot2); library(patchwork)

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
    
    #causal, MT=2 for beta2; MT=3 for beta1
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
#Supp Figure: Structure of simulated genomes

gsum = function(inputs, title) {
  resNFDS = readRDS(inputs[1])
  resPFDS = readRDS(inputs[2])
  resOD = readRDS(inputs[3])
  resSTVS = readRDS(inputs[4])
  
  lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
  lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
  b1 = ggplot(NULL,aes(x=lab,y=c(resNFDS$snp,resPFDS$snp,resOD$snp,resSTVS$snp))) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2,pch=16) + theme_classic() + ylab("No. of SNPs") + xlab("") + ylim(0,NA) + ggtitle(title)
  b2 = ggplot(NULL,aes(x=lab,y=c(resNFDS$MAF,resPFDS$MAF,resOD$MAF,resSTVS$MAF))) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2,pch=16) + theme_classic() + ylab("Mean MAF") + xlab("") + ylim(0,NA)
  b3 = ggplot(NULL,aes(x=lab,y=c(resNFDS$Hz,resPFDS$Hz,resOD$Hz,resSTVS$Hz))) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2,pch=16) + theme_classic() + ylab("Mean Ht") + xlab("") + ylim(0,NA)
  b4 = ggplot(NULL,aes(x=lab,y=c(resNFDS$Gst,resPFDS$Gst,resOD$Gst,resSTVS$Gst))) + geom_boxplot(outlier.shape=NA) +
    geom_jitter(alpha=0.2,pch=16) + theme_classic() + ylab("Mean Gst") + xlab("") + ylim(0,NA)
  return(b1|b2|b3|b4)
}

inputs = c("output/popStatNFDS_self.rds",
           "output/popStatPFDS_self.rds",
           "output/popStatOD_self.rds",
           "output/popStatSTVS_self.rds")
bp1 = gsum(inputs=inputs,title=expression("(a) Non-zero "*beta[1]))

inputs = c("output/popStatNFDS.rds",
           "output/popStatPFDS.rds",
           "output/popStatOD.rds",
           "output/popStatSTVS.rds")
bp2 = gsum(inputs=inputs,title=expression("(b) Non-zero "*beta[2]))


inputs = c("output/popStatNFDS_all.rds",
           "output/popStatPFDS_all.rds",
           "output/popStatOD_all.rds",
           "output/popStatSTVS_all.rds")
bp3 = gsum(inputs=inputs,title="(c) All SNPs")

bp = bp1 / bp2 / bp3
ggsave(bp,filename="SimGenomeSummary.pdf",width=12,height=8)


###########################
# Main Figure: LMM of beta_2

#### continuous
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res1 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res2 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res3 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res4 = roc(ans~p)

roc_p1 = ggplot(NULL,aes(x=1-roc_res1$specificities,y=roc_res1$sensitivities)) + 
  geom_line(colour="blue") + labs(title="Continuous setting") +
  geom_line(aes(x=1-roc_res2$specificities,y=roc_res2$sensitivities),colour="red") +
  geom_line(aes(x=1-roc_res3$specificities,y=roc_res3$sensitivities),colour=grey(0.25,0.75)) +
  geom_line(aes(x=1-roc_res4$specificities,y=roc_res4$sensitivities),colour="grey") +
  theme_classic() + ylab("True positive rate") + xlab("False positive rate") + geom_abline(slope=1, intercept=0,lty=2)


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

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
auc_p1 = ggplot(NULL,aes(x=lab,y=c(auc1,auc2,auc3,auc4))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + ylim(0.4,1) + geom_hline(yintercept=0.5,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab("Model performance") + xlab("")


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

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
beta_p1 = ggplot(NULL,aes(x=lab,y=c(beta1,beta2,beta3,beta4))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + geom_hline(yintercept=0.0,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab(expression("Estimated "*beta[2])) + xlab("") + 
  geom_text(aes(x=c(1:4,4),y=c(-0.1,0.1,0.1,-0.1,0.1)),label="X",size=4)

p1 = roc_p1 + auc_p1 + beta_p1


#### split
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res5 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res6 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res7 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res8 = roc(ans~p)

roc_p2 = ggplot(NULL,aes(x=1-roc_res5$specificities,y=roc_res5$sensitivities)) + 
  geom_line(colour="blue") + labs(title="Split setting") +
  geom_line(aes(x=1-roc_res6$specificities,y=roc_res6$sensitivities),colour="red") +
  geom_line(aes(x=1-roc_res7$specificities,y=roc_res7$sensitivities),colour=grey(0.25,0.75)) +
  geom_line(aes(x=1-roc_res8$specificities,y=roc_res8$sensitivities),colour="grey") +
  theme_classic() + ylab("True positive rate") + xlab("False positive rate") + geom_abline(slope=1, intercept=0,lty=2)


## auc
auc5=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  auc5 = c(auc5,results$AUC[2])
}

auc6=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  auc6 = c(auc6,results$AUC[2])
}

auc7=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  auc7 = c(auc7,results$AUC[2])
}

auc8=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  auc8 = c(auc8,results$AUC[2])
}

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
auc_p2 = ggplot(NULL,aes(x=lab,y=c(auc5,auc6,auc7,auc8))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + ylim(0.4,1) + geom_hline(yintercept=0.5,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab("Model performance") + xlab("")

## beta
beta5=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMMsplit.rds"))
  beta5 = c(beta5,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta6=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMMsplit.rds"))
  beta6 = c(beta6,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta7=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMMsplit.rds"))
  beta7 = c(beta7,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta8=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMMsplit.rds"))
  beta8 = c(beta8,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
beta_p2 = ggplot(NULL,aes(x=lab,y=c(beta5,beta6,beta7,beta8))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + geom_hline(yintercept=0.0,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab(expression("Estimated "*beta[2])) + xlab("") + 
  geom_text(aes(x=c(1:4,4),y=c(-0.1,0.1,0.1,-0.1,0.1)),label="X",size=4)

p2 = roc_p2 + auc_p2 + beta_p2
p = (p2 / p1) + plot_annotation(tag_levels="a")
ggsave(p,filename="beta2LMMdomi.pdf",width=10,height=6)


########################
# Supp Figure: LM of beta2

#### continuous
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res1 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res2 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res3 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMcontinuous.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res4 = roc(ans~p)

roc_p1 = ggplot(NULL,aes(x=1-roc_res1$specificities,y=roc_res1$sensitivities)) + 
  geom_line(colour="blue") + labs(title="Continuous setting") +
  geom_line(aes(x=1-roc_res2$specificities,y=roc_res2$sensitivities),colour="red") +
  geom_line(aes(x=1-roc_res3$specificities,y=roc_res3$sensitivities),colour=grey(0.25,0.75)) +
  geom_line(aes(x=1-roc_res4$specificities,y=roc_res4$sensitivities),colour="grey") +
  theme_classic() + ylab("True positive rate") + xlab("False positive rate") + geom_abline(slope=1, intercept=0,lty=2)


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

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
auc_p1 = ggplot(NULL,aes(x=lab,y=c(auc1,auc2,auc3,auc4))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + ylim(0.4,1) + geom_hline(yintercept=0.5,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab("Model performance") + xlab("")


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

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
beta_p1 = ggplot(NULL,aes(x=lab,y=c(beta1,beta2,beta3,beta4))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + geom_hline(yintercept=0.0,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab(expression("Estimated "*beta[2])) + xlab("") + 
  geom_text(aes(x=c(1:4,4),y=c(-0.1,0.1,0.1,-0.1,0.1)),label="X",size=4)

p1 = roc_p1 + auc_p1 + beta_p1


#### split
## ROC curve
p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res5 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res6 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res7 = roc(ans~p)

p=c(); ans=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  p = c(p, -log10(results$gwas_out$p_nei))
  ans = c(ans, results$ans[,2])
}
roc_res8 = roc(ans~p)

roc_p2 = ggplot(NULL,aes(x=1-roc_res5$specificities,y=roc_res5$sensitivities)) + 
  geom_line(colour="blue") + labs(title="Split setting") +
  geom_line(aes(x=1-roc_res6$specificities,y=roc_res6$sensitivities),colour="red") +
  geom_line(aes(x=1-roc_res7$specificities,y=roc_res7$sensitivities),colour=grey(0.25,0.75)) +
  geom_line(aes(x=1-roc_res8$specificities,y=roc_res8$sensitivities),colour="grey") +
  theme_classic() + ylab("True positive rate") + xlab("False positive rate") + geom_abline(slope=1, intercept=0,lty=2)


## auc
auc5=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  auc5 = c(auc5,results$AUC[2])
}

auc6=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  auc6 = c(auc6,results$AUC[2])
}

auc7=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  auc7 = c(auc7,results$AUC[2])
}

auc8=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  auc8 = c(auc8,results$AUC[2])
}

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
auc_p2 = ggplot(NULL,aes(x=lab,y=c(auc5,auc6,auc7,auc8))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + ylim(0.4,1) + geom_hline(yintercept=0.5,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab("Model performance") + xlab("")

## beta
beta5=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_NFDS_LMsplit.rds"))
  beta5 = c(beta5,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta6=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_PFDS_LMsplit.rds"))
  beta6 = c(beta6,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta7=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_OD_LMsplit.rds"))
  beta7 = c(beta7,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

beta8=c()
for(i in 1:30) {
  results = readRDS(paste0("output/",i,"_STVS_LMsplit.rds"))
  beta8 = c(beta8,median(results$gwas_out$beta_nei[results$ans[,2]==1]))
}

lab = c(rep("NFDS",30), rep("PFDS",30), rep("OD",30), rep("STVS",30))
lab = factor(lab,levels=c("NFDS","PFDS","OD","STVS"))
beta_p2 = ggplot(NULL,aes(x=lab,y=c(beta5,beta6,beta7,beta8))) + geom_boxplot(outlier.shape = NA,fill=c("skyblue","red",grey(0.25,0.75),"grey")) +
  theme_classic() + geom_hline(yintercept=0.0,lty=2) + geom_jitter(alpha=0.5,pch=16) +
  ylab(expression("Estimated "*beta[2])) + xlab("") + 
  geom_text(aes(x=c(1:4,4),y=c(-0.1,0.1,0.1,-0.1,0.1)),label="X",size=4)

p2 = roc_p2 + auc_p2 + beta_p2
p = (p2 / p1) + plot_annotation(tag_levels="a")
ggsave(p,filename="beta2LMdomi.pdf",width=10,height=6)

