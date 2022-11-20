# GWAS for SLiM output .vcf file
set.seed(1234)

library(gaston)
library(vcfR)
library(rNeighborGWAS)
library(pROC)

runGWAS = function(type=c("NFDS","PFDS","OD","STVS"), setting=c("split","continuous"), iter=30) {
  for(i in 1:iter) {
    f_name = paste0(i,"_",type)
    f_path = paste0("output/",f_name,".vcf.gz")
    
    d0 = read.vcf(f_path)
    d = select.snps(d0, d0@snps$maf>0.01)
    
    ann = read.vcfR(f_path)
    ann = ann[which(d0@snps$maf>0.01),]
    
    geno = gaston2neiGWAS(d)
    
    g = geno$geno[sample(1:nrow(geno$geno),900,replace=FALSE),]
    beta1 = rep(0,nrow(d@snps))
    n_causal1 = length(grep("MT=3",ann@fix[,"INFO"]))
    beta1[grep("MT=3",ann@fix[,"INFO"])] = 0.1
    
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
    
    if(type=="NFDS") {
      beta2 = rep(0,nrow(d@snps))
      n_causal2 = length(grep("MT=2",ann@fix[,"INFO"]))
      beta2[grep("MT=2",ann@fix[,"INFO"])] = - 0.1 #NFDS
      
      y = as.numeric((g%*%beta1 + n_causal1*0.1) + (g_nei_domi%*%beta2 + 1.5))
      y = y + rnorm(length(y),0,0.75*sd(y)) 
    } else if(type=="PFDS") {
      beta2 = rep(0,nrow(d@snps))
      n_causal2 = length(grep("MT=2",ann@fix[,"INFO"]))
      beta2[grep("MT=2",ann@fix[,"INFO"])] = 0.1 #PFDS
      
      y = as.numeric((g%*%beta1 + n_causal1*0.1) + (g_nei_domi%*%beta2 + 1.5))
      y = y + rnorm(length(y),0,0.75*sd(y)) 
    } else if(type=="OD") {
      beta2 = rep(0,nrow(d@snps))
      n_causal2 = length(grep("MT=2",ann@fix[,"INFO"]))
      beta2 = g[,grep("MT=2",ann@fix[,"INFO"])] #OD
      beta2[beta2==-1] = 1 - 0.1
      beta2[beta2==1] = 1 + 0.1
      beta2[beta2==0] = 1 + 0.1*2
      beta2 = apply(beta2,1,sum)/ncol(beta2)
      
      y = as.numeric((g%*%beta1 + n_causal1*0.1) + beta2)
      y = y + rnorm(length(y),0,0.75*sd(y)) 
    } else {
      beta2 = rep(0,nrow(d@snps))
      n_causal2 = length(grep("MT=2",ann@fix[,"INFO"]))
      beta2[grep("MT=2",ann@fix[,"INFO"])] = sample(c(-0.1,0.1),length(grep("MT=2",ann@fix[,"INFO"])),replace=TRUE)
      
      y = as.numeric((g%*%beta1+n_causal1*0.1) + (g_domi%*%beta2+n_causal2*0.1))
      y = y + rnorm(length(y),0,0.75*sd(y)) 
    }

    # LMM
    res = nei_lmm(g_domi,g_nei_domi,y) #replacing nei_lmm() by nei_lm() gives LM results

    ans1 = rep(0,nrow(d@snps))
    ans1[grep("MT=3",ann@fix[,"INFO"])] = 1
    p1 = -log10(res$p_self)
    roc1 = roc(ans1~p1)$auc

    ans2 = rep(0,nrow(d@snps))
    ans2[grep("MT=2",ann@fix[,"INFO"])] = 1
    p2 = -log10(res$p_nei)
    roc2 = roc(ans2~p2)$auc

    results = list()
    results[[1]] = c(roc1, roc2)
    results[[2]] = c(n_causal1, n_causal2)
    results[[3]] = cbind(ans1, ans2)
    results[[4]] = res
    names(results) = c("AUC", "n_causal", "ans", "gwas_out")
    print(results$AUC)
    saveRDS(results,file=paste0("output/", f_name, "_LMM", setting, ".rds"))
  }
  return(message("Done. Output files are in ./output"))
}

runGWAS(type="NFDS",setting="continuous",iter=30)
runGWAS(type="NFDS",setting="split",iter=30)

runGWAS(type="PFDS",setting="continuous",iter=30)
runGWAS(type="PFDS",setting="split",iter=30)

runGWAS(type="OD",setting="continuous",iter=30)
runGWAS(type="OD",setting="split",iter=30)

runGWAS(type="STVS",setting="continuous",iter=30)
runGWAS(type="STVS",setting="split",iter=30)
