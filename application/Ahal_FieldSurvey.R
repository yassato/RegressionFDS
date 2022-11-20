# Single-marker analysis using Arabidopsis halleri 

######
# data analysis

# load library
library(lme4)
library(ggplot2)

# field data collected at the Omoide River (OM) site, Hyogo, Japan.
# original data are available via Dryad (https://doi.org/10.5061/dryad.53k2d)
OMcensus = read.csv("OMcensusAll2.csv",header=T)

OMcensus$Trichome = as.character(OMcensus$Trichome)
OMcensus$Trichome = as.character(OMcensus$Trichome)

OMcensus$Trichome[OMcensus$Trichome=="h"] = 1 # hairy morph as dominant
OMcensus$Trichome[OMcensus$Trichome=="g"] = -1 # glabrous morph as recessive
OMcensus$Trichome = as.numeric(OMcensus$Trichome)

# calc. morph similarity
xixj = c()
for(i in 1:nrow(OMcensus)) {
  sub = OMcensus[-i,] #Note: if individual i is eliminated, 'sim' varies due to the limited number of individuals
  sub = subset(sub, Year==OMcensus[i,"Year"]&PatchID==OMcensus[i,"PatchID"])
  xixj = c(xixj,OMcensus$Trichome[i]*mean(sub$Trichome))
}

OMcensus = data.frame(OMcensus,xixj)

# Poisson GLMMs of the flower number
# Tri. x Morph similarity is not significant at P<0.05
summary(glmer(Flower~Trichome*xixj+Total+(1|Year/PatchID), data=OMcensus, family=poisson, 
              offset=log(MaxLeaf+1), control=glmerControl(optimizer="bobyqa")))
summary(glmer(Flower~Trichome+xixj+Total+(1|Year/PatchID), data=OMcensus, family=poisson,
              offset=log(MaxLeaf+1), control=glmerControl(optimizer="bobyqa")))

#####
# Figure

# observed
d2 = aggregate(log((OMcensus$Flower/OMcensus$MaxLeaf)+0.1)~Year+factor(Trichome)+OMcensus$gProp,OMcensus,mean)
pch = rep(NA,length(d2$`factor(Trichome)`))
pch[d2$`factor(Trichome)`==-1] = 1; pch[d2$`factor(Trichome)`==1] = 16
hFreq = 1 - d2$`OMcensus$gProp`

# predicted (see Appendix S2)
# b0 = -2.076348; b1 = 0.146011; b2 = -0.163216; b12 = -0.087591
# f_star = 0.5-(b1/(2*b2))

plt_Ah = function(b0,b1,b2,b12) {
  f_star = 0.5-(b1/(2*b2))
  p1 = ggplot(NULL,aes(x=hFreq,y=d2$`log((OMcensus$Flower/OMcensus$MaxLeaf) + 0.1)`)) + theme_classic() +
    geom_point(pch=pch,col="grey") + xlab("phenotype-level frequency of hairy morph") + ylab("Flower production") + 
    geom_function(aes(x=1,y=1),fun=function(x,b0,b1,b2,b12) { (b12+b2)*(2*x-1)+b0+b1 }, args=list(b0,b1,b2,b12)) + 
    geom_function(aes(x=1,y=1),fun=function(x,b0,b1,b2,b12) { (b12-b2)*(2*x-1)+b0-b1 }, args=list(b0,b1,b2,b12), colour=grey(0.0,0.33)) +
    geom_function(aes(x=1,y=1),fun=function(x,b0,b1,b2,b12) { 2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+b0-b1+2*b1*x }, args=list(b0,b1,b2,b12), lty=2) +
    geom_point(aes(x=f_star,y=(b12+b2)*(2*f_star-1)+b0+b1),pch=16,size=3) + labs(title=substitute(paste(italic('Arabidopsis halleri'))))
  return(p1)
}

p1 = plt_Ah(b0 = -2.076348, b1 = 0.146011, b2 = -0.163216, b12 = -0.087591)
 
saveRDS(p1,file="Ahal_plot.rds")
