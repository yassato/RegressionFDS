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
summary(glmer(Flower~factor(Trichome)*xixj+Total+(1|Year/PatchID), data=OMcensus, family=poisson, 
              offset=log(MaxLeaf+1), control=glmerControl(optimizer="bobyqa")))
summary(glmer(Flower~factor(Trichome)+xixj+Total+(1|Year/PatchID), data=OMcensus, family=poisson,
              offset=log(MaxLeaf+1), control=glmerControl(optimizer="bobyqa")))

