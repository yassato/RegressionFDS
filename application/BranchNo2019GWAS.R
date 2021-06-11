#####################
#2019 Branch No GWAS#
#####################

library(rNeighborGWAS)
source("nei_lmm.R")
source("nei_lm.R")

#load data
geno_d = readRDS(file="sub_snpMAF5.rds")
geno_d[geno_d==0]=-1 #replace 0 into -1

position = readRDS(file="positionsMAF5.rds")

pheno_d = read.csv("./Survey2019CHZ4GWAS_BranchNo.csv",header=T)
pheno_d$gwasID = paste0("X",pheno_d$gwasID)

#reshape pheno.data
pheno_d = subset(pheno_d, gwasID!="X7329")
pheno_d = subset(pheno_d, gwasID!="X3")
naID = which(is.na(pheno_d$InitLeafLen))
pheno_d = pheno_d[-naID,]

n_marker = nrow(geno_d)
n_plants = nrow(pheno_d)

geno = geno_d[,as.character(pheno_d$gwasID)]
geno = t(geno)
n_plants == nrow(geno)

smap = cbind(pheno_d$position_X,pheno_d$position_Y)
scale = sqrt(2)+0.01

rm(geno_d)
gc();gc()

g_nei = nei_coval(geno=geno,smap=smap,scale=scale,grouping=pheno_d$Block)
X = as.matrix(model.matrix(~factor(Block)+scale(InitLeafLen)+Bolting+edge-1,data=pheno_d))
Y = scale(log(pheno_d$BranchNo+1))

res = nei_lmm(geno=geno,g_nei=g_nei,pheno=Y,addcovar=X,n_core=24L,asym=TRUE)
res = data.frame(position,res)
colnames(res) = c("Chr","Position","MAF","beta_self","beta_nei","beta_sxn","P_self","P_nei","P_sxn")
write.csv(res,"210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMMasym.csv")
gc();gc()

res = nei_lm(geno=geno,g_nei=g_nei,pheno=Y,addcovar=X,n_core=24L,asym=TRUE)
res = data.frame(position,res)
colnames(res) = c("Chr","Position","MAF","beta_self","beta_nei","beta_sxn","P_self","P_nei","P_sxn")
write.csv(res,"210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMasym.csv")
gc();gc()

res = nei_lmm(geno=geno,g_nei=g_nei,pheno=Y,addcovar=X,n_core=24L)
res = data.frame(position,res)
colnames(res) = c("Chr","Position","MAF","beta_self","beta_nei","P_self","P_nei")
write.csv(res,"210604_BranchNoNeiGWASbolting_CHZ2019_scaledLMM.csv")
gc();gc()

res = nei_lm(geno=geno,g_nei=g_nei,pheno=Y,addcovar=X,n_core=24L)
res = data.frame(position,res)
colnames(res) = c("Chr","Position","MAF","beta_self","beta_nei","P_self","P_nei")
write.csv(res,"210604_BranchNoNeiGWASbolting_CHZ2019_scaledLM.csv")
gc();gc()

