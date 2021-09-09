# Single-marker analysis using Arabidopsis halleri 

# load library
library(lme4)

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

# Figure
svg("Ahal_plot.svg",width=5,height=5)

# observed
d2 = aggregate(log((OMcensus$Flower/OMcensus$MaxLeaf)+0.1)~Year+factor(Trichome)+OMcensus$gProp,OMcensus,mean)
pch = rep(NA,length(d2$`factor(Trichome)`))
pch[d2$`factor(Trichome)`==-1] = 1; pch[d2$`factor(Trichome)`==1] = 16
hFreq = 1 - d2$`OMcensus$gProp`
plot(d2$`log((OMcensus$Flower/OMcensus$MaxLeaf) + 0.1)`~hFreq,pch=pch,las=1,col="grey",
     ylab="log(no. of flowers / mm)",xlab="phenotype-level frequency of hairy morph",cex.lab=1.2)

# predicted (see Appendix S3)
b0 = -2.076348; b1 = 0.146011; b2 = -0.163216; b12 = -0.087591
curve((b12+b2)*(2*x-1)+b0+b1,ylim=c(-2.5,-1.5),add=T,
      las=1,ylab="Absolute fitness",xlab="Frequency of A alleles",
      cex.axis=1.2,cex.lab=1.2,col=grey(0.0,1.0),lwd=1.5)
curve((b12-b2)*(2*x-1)+b0-b1,add=T,col=grey(0.0,0.33),lwd=1.5)
curve(2*b2*x*(2*x-1)+(b12-b2)*(2*x-1)+b0-b1+2*b1*x,add=T,lty=2,lwd=1.5)
f_star = 0.5-(b1/(2*b2))
points(f_star,(b12+b2)*(2*f_star-1)+b0+b1,pch=16,cex=1.5)

dev.off()
