# Single-marker analysis on a damselfly, Ischnura elegans

# load library
library(lme4)

# data available upon request for Y. Takahashi (takahashi.yum@gmail.com)
d = read.csv("fitness_data_forSato.csv", header=TRUE)

d$morph[d$morph=="i"] = -1 #infuscans-type gynomorph as ressesive
d$morph[d$morph=="a"] = +1 #andromorph as dominant
d$morph = as.numeric(d$morph)

# exchange frequency conditions into morph similarities
f_sim = c(NA, nrow(d))
f_sim[(d$morph==1)&(d$frequency=="0.8")] = 0.6
f_sim[(d$morph==1)&(d$frequency=="0.5")] = 0.0
f_sim[(d$morph==1)&(d$frequency=="0.2")] = -0.6

f_sim[(d$morph==-1)&(d$frequency=="0.8")] = -0.6
f_sim[(d$morph==-1)&(d$frequency=="0.5")] = 0.0
f_sim[(d$morph==-1)&(d$frequency=="0.2")] = 0.6

d = data.frame(d,f_sim)

# Poisson GLMMs
summary(glmer(matureeggnumber~factor(morph)*f_sim+density+(1|cageID/experimentID),data=d,family=poisson)) #morph x f_sim not signif.
summary(glmer(matureeggnumber~factor(morph)+f_sim+density+(1|cageID/experimentID),data=d,family=poisson))
