# load library and functions
library(ggplot2); library(patchwork)

coord = function(chr, pos) {
  if(length(pos)!=length(chr)) stop("chr and pos length differ")
  chr <- as.factor(chr)
  coord <- 0
  M <- 0
  for (i in 1:nlevels(chr)) {
    w <- (chr == levels(chr)[i])
    pos.c <- pos[w]
    coord[w] <- M + pos.c
    mx <- max(pos.c)
    M <- M + mx
  }
  coord <- coord/M
  return(coord)
}

######
# LMM
gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLMM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

x = coord(gwas_out$Chr,gwas_out$Position)
man1 = ggplot(data.frame(gwas_out,x),aes(x=x,y=-log10(gwas_out$P_self))) + geom_point(colour=cols) +
  theme_classic() + ylab(expression(-log[10](p))) + xlab("Chromosomes") + 
  ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) + 
  geom_hline(yintercept=-log10(p_adj),lty=2,col=grey(0.5,0.5)) + ggtitle(expression("self effects "*beta[1]))

x = coord(gwas_out$Chr,gwas_out$Position)
man2 = ggplot(data.frame(gwas_out,x),aes(x=x,y=-log10(gwas_out$P_nei))) + geom_point(colour=cols) +
  theme_classic() + ylab(expression(-log[10](p))) + xlab("Chromosomes") + 
  ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) + 
  geom_hline(yintercept=-log10(p_adj),lty=2,col=grey(0.5,0.5)) + ggtitle(expression("genotype similarity effects "*beta[2]))

h1 = ggplot(subset(gwas_out,P_nei<0.0001),aes(x=beta_nei)) + geom_histogram() + 
  theme_classic() + ylab("No. pf SNPs") + xlab(expression("Estimated "*beta[2]))

gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLMMasym.csv", header=TRUE)
man3 = ggplot(data.frame(gwas_out,x),aes(x=x,y=-log10(gwas_out$P_sxn))) + geom_point(colour=cols) +
  theme_classic() + ylab(expression(-log[10](p))) + xlab("Chromosomes") + 
  ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) + 
  geom_hline(yintercept=-log10(p_adj),lty=2,col=grey(0.5,0.5)) + ggtitle(expression("asymmetric effects "*beta[12]))

mh = ((man1 / man2 / man3) | h1) + plot_layout(widths=c(3,1)) + plot_annotation(tag_levels="a")
ggsave(mh,filename="ManhattanLMM.png",width=12,height=6,dpi=600)

#####
# QQ-plot LMM
gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLMM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

q1 = ggplot(data=gwas_out, mapping=aes(x=-log(ppoints(length(P_self)),10),y=-log(sort(P_self,decreasing=FALSE),10))) +
  geom_point(colour="grey") +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_classic() + ggtitle(expression("self effects "*beta[1])) +
  xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

q2 = ggplot(data=gwas_out, mapping=aes(x=-log(ppoints(length(P_nei)),10),y=-log(sort(P_nei,decreasing=FALSE),10))) +
  geom_point(colour="grey") +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_classic() + ggtitle(expression("genotype similarity effects "*beta[2])) +
  xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLMMasym.csv", header=TRUE)
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

q3 = ggplot(data=gwas_out, mapping=aes(x=-log(ppoints(length(P_sxn)),10),y=-log(sort(P_sxn,decreasing=FALSE),10))) +
  geom_point(colour="grey") +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_classic() + ggtitle(expression("asymmetric effects "*beta[12])) +
  xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

q = q1 + q2 + q3 + plot_annotation(tag_levels="a")
ggsave(q,filename="QQplotLMM.png",width=9,height=3,dpi=600)

#####
# LM
gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

x = coord(gwas_out$Chr,gwas_out$Position)
man1 = ggplot(data.frame(gwas_out,x),aes(x=x,y=-log10(gwas_out$P_self))) + geom_point(colour=cols) +
  theme_classic() + ylab(expression(-log[10](p))) + xlab("Chromosomes") + 
  ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) + 
  geom_hline(yintercept=-log10(p_adj),lty=2,col=grey(0.5,0.5)) + ggtitle(expression("self effects "*beta[1]))

x = coord(gwas_out$Chr,gwas_out$Position)
man2 = ggplot(data.frame(gwas_out,x),aes(x=x,y=-log10(gwas_out$P_nei))) + geom_point(colour=cols) +
  theme_classic() + ylab(expression(-log[10](p))) + xlab("Chromosomes") + 
  ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) + 
  geom_hline(yintercept=-log10(p_adj),lty=2,col=grey(0.5,0.5)) + ggtitle(expression("genotype similarity effects "*beta[2]))

h1 = ggplot(subset(gwas_out,P_nei<0.0001),aes(x=beta_nei)) + geom_histogram() + 
  theme_classic() + ylab("No. pf SNPs") + xlab(expression("Estimated "*beta[2]))

gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLMasym.csv", header=TRUE)
man3 = ggplot(data.frame(gwas_out,x),aes(x=x,y=-log10(gwas_out$P_sxn))) + geom_point(colour=cols) +
  theme_classic() + ylab(expression(-log[10](p))) + xlab("Chromosomes") + 
  ggplot2::theme(axis.ticks.x=ggplot2::element_blank(),axis.text.x=ggplot2::element_blank()) + 
  geom_hline(yintercept=-log10(p_adj),lty=2,col=grey(0.5,0.5)) + ggtitle(expression("asymmetric effects "*beta[12]))

mh = ((man1 / man2 / man3) | h1) + plot_layout(widths=c(3,1)) + plot_annotation(tag_levels="a")
ggsave(mh,filename="ManhattanLM.png",width=12,height=6,dpi=600)

#####
# QQ-plot LMM
gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLM.csv", header=TRUE)
p_adj = 0.05/(nrow(gwas_out))
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

q1 = ggplot(data=gwas_out, mapping=aes(x=-log(ppoints(length(P_self)),10),y=-log(sort(P_self,decreasing=FALSE),10))) +
  geom_point(colour="grey") +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_classic() + ggtitle(expression("self effects "*beta[1])) +
  xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

q2 = ggplot(data=gwas_out, mapping=aes(x=-log(ppoints(length(P_nei)),10),y=-log(sort(P_nei,decreasing=FALSE),10))) +
  geom_point(colour="grey") +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_classic() + ggtitle(expression("genotype similarity effects "*beta[2])) +
  xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

gwas_out = read.csv("BranchNoNeiGWASbolting_CHZ2019_scaledLMasym.csv", header=TRUE)
chr_rep = cumsum(table(gwas_out$Chr))
cols = c(rgb(1,0,0, 2*gwas_out$MAF[1:chr_rep[1]]), rgb(0,1,0, 2*gwas_out$MAF[(chr_rep[1]+1):(chr_rep[2])]), rgb(0,0,1, 2*gwas_out$MAF[(chr_rep[2]+1):(chr_rep[3])]), rgb(0,0,0, 2*gwas_out$MAF[(chr_rep[3]+1):(chr_rep[4])]), rgb(1,0,1, 2*gwas_out$MAF[(chr_rep[4]+1):(chr_rep[5])]))

q3 = ggplot(data=gwas_out, mapping=aes(x=-log(ppoints(length(P_sxn)),10),y=-log(sort(P_sxn,decreasing=FALSE),10))) +
  geom_point(colour="grey") +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_classic() + ggtitle(expression("asymmetric effects "*beta[12])) +
  xlab(expression("Expected "*-log[10](p))) + ylab(expression("Observed "*-log[10](p)))

q = q1 + q2 + q3 + plot_annotation(tag_levels="a")
ggsave(q,filename="QQplotLM.png",width=9,height=3,dpi=600)
