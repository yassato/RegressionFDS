# script to create toy data
set.seed(1234)

library(vcfR)
library(gaston)
library(rNeighborGWAS)

simy = function(setting=c("split","continuous"),b1,b2) {
  d0 = read.vcf("output/1_NFDS.vcf.gz")
  d = select.snps(d0, d0@snps$maf>0.01)
  
  ann = read.vcfR("output/1_NFDS.vcf.gz")
  ann = ann[which(d0@snps$maf>0.01),]
  
  geno = gaston2neiGWAS(d)
  g = geno$geno[sample(1:nrow(geno$geno),900,replace=FALSE),]
  beta1 = rep(0,nrow(d@snps))
  n_causal1 = length(grep("MT=3",ann@fix[,"INFO"]))
  beta1[grep("MT=3",ann@fix[,"INFO"])] = b1 # b1 = 0.0: no polygenic effects
  
  if(setting=="split") {
    smap = cbind(1:nrow(g),1:nrow(g))
    group = rep((1:(nrow(g)/10)),each=10)
    g_nei = nei_coval(g,smap,scale=20,grouping=group)
    
    g_domi = g
    g_domi[g_domi==0] = 1 # commenting out this line assumes additive effects
    g_nei_domi = nei_coval(g_domi,smap,scale=20,grouping=group)
  } else {
    smap = cbind(rep(1:(nrow(g)/30),each=nrow(g)/30),
                 rep(1:(nrow(g)/30),nrow(g)/30)
    ) #30 x 30 map
    g_nei = nei_coval(g,smap,scale=1.42)
    
    g_domi = g
    g_domi[g_domi==0] = 1 # commenting out this line assumes additive effects
    g_nei_domi = nei_coval(g_domi,smap,scale=1.42)
  }
  
  beta2 = rep(0,nrow(d@snps))
  n_causal2 = length(grep("MT=2",ann@fix[,"INFO"]))
  beta2[grep("MT=2",ann@fix[,"INFO"])][1] = b2 #
  
  y = as.numeric((g%*%beta1 + n_causal1*0.1) + (g_nei_domi%*%beta2 + 1.5))
  y = y + rnorm(length(y),0,0.75*sd(y))
  
  if(setting=="split") {
    if(b1==0) {
      ID = rownames(g)
      gi = rep(NA,length(y))
      gi[g[,which(beta2!=0)]==1] = "A"
      gi[g[,which(beta2!=0)]==0] = "A"
      gi[g[,which(beta2!=0)]==-1] = "a"
      pheno = data.frame(ID,y,gi,group)
    } else {
      ID = rownames(g)
      pheno = data.frame(ID,y,group)
    }
  }
  
  if(setting=="continuous") {
    if(b1==0) {
      ID = rownames(g)
      gi = rep(NA,length(y))
      gi[g[,which(beta2!=0)]==1] = "A"
      gi[g[,which(beta2!=0)]==0] = "A"
      gi[g[,which(beta2!=0)]==-1] = "a"
      X = smap[,1]
      Y = smap[,2]
      pheno = data.frame(ID,y,gi,X,Y)
    } else {
      ID = rownames(g)
      X = smap[,1]
      Y = smap[,2]
      pheno = data.frame(ID,y,X,Y)
    }
  }
  return(pheno)
}

d1 = simy("split",b1=0,b2=-0.2); saveRDS(d1,file="./output/toy1.rds")
d2 = simy("continuous",b1=0,b2=-0.2); saveRDS(d2,file="./output/toy2.rds")
d3 = simy("split",b1=0.1,b2=-0.2); saveRDS(d3,file="./output/toy3.rds")
d4 = simy("continuous",b1=0.1,b2=-0.2); saveRDS(d4,file="./output/toy4.rds")

# toy4
d4 = readRDS("./output/toy4.rds")

g = gaston::read.vcf("output/1_NFDS.vcf.gz")
g = gaston::select.snps(g, g@snps$maf>0.01)
g = rNeighborGWAS::gaston2neiGWAS(g)
geno = g$geno[d4$ID,]
geno[geno==0] = 1
smap = cbind(d4$X,d4$Y)
g_nei = rNeighborGWAS::nei_coval(geno,smap,scale=sqrt(2)+0.01)
res = rNeighborGWAS::nei_lmm(geno,g_nei,d4$y)

x = data.frame(g$gmap,res$p_nei)
colnames(x) = c("chr","pos","p")
gaston::manhattan(x,las=1,thinning=FALSE)
abline(h=-log10(0.05/nrow(g_nei)),lty=2,col="grey")


# toy3
d3 = readRDS("./output/toy3.rds")

g = gaston::read.vcf("output/1_NFDS.vcf.gz")
g = gaston::select.snps(g, g@snps$maf>0.01)
g = rNeighborGWAS::gaston2neiGWAS(g)
geno = g$geno[d3$ID,]
geno[geno==0] = 1
smap = cbind(runif(nrow(d3),0,1),runif(nrow(d3),0,1)) # random dummy map
g_nei = rNeighborGWAS::nei_coval(geno,smap,scale=10^6,grouping=d3$group)
res = rNeighborGWAS::nei_lmm(geno,g_nei,d3$y,addcovar=model.matrix(~d3$group))

x = data.frame(g$gmap,res$p_nei)
colnames(x) = c("chr","pos","p")
gaston::manhattan(x,las=1,thinning=FALSE)
abline(h=-log10(0.05/nrow(g_nei)),lty=2,col="grey")

