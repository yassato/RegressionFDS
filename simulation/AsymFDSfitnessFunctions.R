# Figure presentation for fitness functions
library(ggplot2)
library(patchwork)

##########
#Main Figure: dominant case

W_AA = function(x,b2,b12) (b12+b2)*(4*x-2*x^2-1)+1
W_Aa = function(x,b2,b12) (b12+b2)*(4*x-2*x^2-1)+1
W_aa = function(x,b2,b12) (b12-b2)*(4*x-2*x^2-1)+1

W_A = function(x,b2,b12) x*W_AA(x,b2,b12)+(1-x)*W_Aa(x,b2,b12)
W_a = function(x,b2,b12) x*W_Aa(x,b2,b12)+(1-x)*W_aa(x,b2,b12)
W_m = function(x,b2,b12) x*W_A(x,b2,b12)+(1-x)*W_a(x,b2,b12)

wplot = function(b2,b12,pch) {
  p = ggplot(NULL,aes(x=0.293,y=0.293*W_A(0.293,b2,b12)+(1-0.293)*W_a(0.293,b2,b12))) + 
    geom_point(pch=pch,size=3) + ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=W_A, col=grey(0.0,1.0),args=list(b2,b12)) +
    geom_function(fun=W_a, col=grey(0.0,0.5),args=list(b2,b12)) + 
    geom_function(fun=W_m, lty=2,args=list(b2,b12)) + 
    ylab("Absolute fitness") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(format(b12,nsmall=1))))
  return(p)
} 

pa = wplot(b2=-0.2,b12=0.0,pch=16)
pa = pa + geom_text(aes(x=0.025,y=1.5),label="----- (black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="----- (gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="- - -  (dashed): mean",size=3,hjust=0)
pb = wplot(b2=-0.2,b12=-0.3,pch=16)
pc = wplot(b2=-0.2,b12=0.3,pch=16)
pd = wplot(b2=0.2,b12=0.0,pch=1)
pe = wplot(b2=0.2,b12=-0.3,pch=1)
pf = wplot(b2=0.2,b12=0.3,pch=1)

p = (pa | pb | pc) / (pd | pe | pf) + plot_annotation(tag_levels = "a")
ggsave(p,filename="AsymFDSdomi.pdf",width=10,height=6)


############
# Supp Figure: additive case
W_AA = function(x,b2,b12) (b12+b2)*(2*x-1)+1
W_Aa = function(x,b2,b12) 1
W_aa = function(x,b2,b12) (b12-b2)*(2*x-1)+1

W_A = function(x,b2,b12) x*W_AA(x,b2,b12)+(1-x)*W_Aa(x,b2,b12)
W_a = function(x,b2,b12) x*W_Aa(x,b2,b12)+(1-x)*W_aa(x,b2,b12)
W_m = function(x,b2,b12) x*W_A(x,b2,b12)+(1-x)*W_a(x,b2,b12)

wplot = function(b2,b12,pch=c(1,16)) {
  p2 = 0.5*(b12-b2)/b12
  w2 = p2*W_A(p2,b2,b12)+(1-p2)*W_a(p2,b2,b12)
  p = ggplot(NULL,aes(x=0.5,y=0.5*W_A(0.5,b2,b12)+(1-0.5)*W_a(0.5,b2,b12))) + 
    geom_point(pch=pch[1],size=3) + ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=W_A, col=grey(0.0,1.0),args=list(b2,b12)) +
    geom_function(fun=W_a, col=grey(0.0,0.5),args=list(b2,b12)) + 
    geom_function(fun=W_m, lty=2,args=list(b2,b12)) + 
    ylab("Absolute fitness") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(format(b12,nsmall=1)))) + 
    geom_point(aes(x=p2,y=w2),pch=pch[2],size=3)
  return(p)
} 

pa = wplot(b2=-0.2,b12=0,pch=c(16,1))
pa = pa + geom_text(aes(x=0.025,y=1.5),label="----- (black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="----- (gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="- - -  (dashed): mean",size=3,hjust=0)
pb = wplot(b2=-0.2,b12=-0.3,pch=c(16,1))
pc = wplot(b2=-0.2,b12=0.3,pch=c(16,1))
pd = wplot(b2=0.2,b12=0,pch=c(1,16))
pe = wplot(b2=0.2,b12=-0.3,pch=c(1,16))
pf = wplot(b2=0.2,b12=0.3,pch=c(1,16))

p = (pa | pb | pc) / (pd | pe | pf) + plot_annotation(tag_levels = "a")
ggsave(p,filename="AsymFDSadd.pdf",width=10,height=6)


###########
# Supp Figure: inbred case

wplot = function(b2,b12,pch) {
  p = ggplot(NULL,aes(x=0.5,y=1)) + 
    geom_point(pch=pch,size=3) + ylim(0.5,1.5) + xlim(0,1) + theme_classic() +
    geom_function(fun=function(x,b2,b12) (b12+b2)*(2*x-1)+1, col=grey(0.0,1.0),args=list(b2,b12)) +
    geom_function(fun=function(x,b2,b12) (b12-b2)*(2*x-1)+1, col=grey(0.0,0.5),args=list(b2,b12)) + 
    geom_function(fun=function(x,b2,b12) 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+1, lty=2,args=list(b2,b12)) + 
    ylab("Absolute fitness") + xlab("Frequency of A alleles") +
    labs(title=bquote(beta[2]*" = "*.(b2)*"; "*beta[12]*" = "*.(format(b12,nsmall=1))))
  return(p)
} 

pa = wplot(b2=-0.2,b12=0.0,pch=16)
pa = pa + geom_text(aes(x=0.025,y=1.5),label="----- (black): A alleles",size=3,hjust=0) + 
  geom_text(aes(x=0.025,y=1.4),label="----- (gray): a alleles",colour=grey(0.5,1),size=3,hjust=0) +
  geom_text(aes(x=0.025,y=1.3),label="- - -  (dashed): mean",size=3,hjust=0)
pb = wplot(b2=-0.2,b12=-0.3,pch=16)
pc = wplot(b2=-0.2,b12=0.3,pch=16)
pd = wplot(b2=0.2,b12=0.0,pch=1)
pe = wplot(b2=0.2,b12=-0.3,pch=1)
pf = wplot(b2=0.2,b12=0.3,pch=1)

p = (pa | pb | pc) /  (pd | pe | pf) + plot_annotation(tag_levels = "a")
ggsave(p,filename="AsymFDSinbred.pdf",width=10,height=6)

