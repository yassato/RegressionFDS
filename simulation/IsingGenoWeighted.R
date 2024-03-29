###########################
# Examples of Ising model #
###########################

# case 1: continuous population

#define functions
XiXj1 = function(id,dmat,range) {
  xi = Xi[id]
  if(xi==0) { xi = 1 } else { xi = xi }
  j = as.numeric(which((dmat[,id]<=range)&(dmat[,id]>0)))
  xj = Xi[j]
  xj[xj==0] = 1
  return(xi*sum(xj))
}

prob = function(id,group) {
  Xik = Xi[group==group[id]]
  f_AA = length(which(Xik==1))/length((Xik))
  f_Aa = length(which(Xik==0))/length((Xik))
  f_aa = length(which(Xik==-1))/length((Xik))
  
  if(Xi[id]==-1) {
    p_AA=f_AA+0.5*f_Aa; p_Aa=f_aa+0.5*f_Aa; p_aa=0
  } else if(Xi[id]==0) {
    p_AA=0.5*f_AA+0.25*f_Aa; p_Aa=0.5*f_AA+0.5*f_Aa+0.5*f_aa; p_aa=0.25*f_Aa+0.5*f_aa
  } else if(Xi[id]==1) {
    p_AA=0; p_Aa=f_AA+0.5*f_Aa; p_aa=0.5*f_Aa+f_aa
  } else {
    message(cat("error at", id))
  }
  return(c(p_AA, p_Aa, p_aa))
}

#set coefficients
J = -1.5 #= beta_2
h = 0.0001 #= beta_1

#set a field
rect = 100
N = rect*rect
Xi = sample(c(-1,1),N,replace=T)

smap = cbind(rep(1:rect,each=rect),
             rep(1:rect,rect)
)
dmat = as.matrix(dist(smap))

group = rep(1,rect*rect)


xixj = c()
for(i in 1:N) xixj = c(xixj, XiXj1(i,dmat=dmat,range=sqrt(8)))
Ei_t0 = (J*xixj+h*Xi)

#MCMC by Gibbs sampling
for(j in 1:100) {
  Xi_t0 = Xi

  perm = sample(1:N,N)
  for(i in perm) {
    Xi[i] = sample(c(-1,0,1),1,prob=prob(i,group))
    xi = Xi[i]
    if(xi==0) { xi = 1 }
    Ei = (J*XiXj1(i,dmat=dmat,range=sqrt(8))+h*xi) #set the energy positive to maximize itself
    
    #Metropolis algorithm
    if(Ei_t0[i]<=Ei) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else if(runif(1,0,1)<exp(Ei-Ei_t0[i])) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else {
      Xi[i] = Xi_t0[i]
    }
  }
  E = mean(Ei_t0)
  image(matrix(Xi,rect,rect),main=paste(j,E),col=c("black","grey","white"))
}


svg(file=paste0("IsingJ",J,"_h",h,"cnt.svg"), width=5, height=5)
image(t(matrix(Xi,rect,rect)),main=paste(j,E),col=c("black","grey","white"),las=1)
dev.off()


# case 2: split population

XiXj2 = function(id,group) {
  xi = Xi[id]
  if(xi==0) { xi = 1 } else { xi = xi }
  xj = Xi[-id][group[-id]==group[id]]
  xj[xj==0] = 1
  return(xi*sum(xj))
}

freq = function(Xi,rect) {
  f=c()
  for(i in 1:rect) {
    f_AA = length(which((Xi==-1)&(group==i)))/rect
    f_Aa = length(which((Xi==0)&(group==i)))/rect
    f_aa = length(which((Xi==1)&(group==i)))/rect
    f = cbind(f, c(f_AA,f_Aa,f_aa))
  }
  return(f)
}

prob = function(id,group) {
  Xik = Xi[group==group[id]]
  f_AA = length(which(Xik==1))/length((Xik))
  f_Aa = length(which(Xik==0))/length((Xik))
  f_aa = length(which(Xik==-1))/length((Xik))
  
  if(Xi[id]==-1) {
    p_AA=f_AA+0.5*f_Aa; p_Aa=f_aa+0.5*f_Aa; p_aa=0
  } else if(Xi[id]==0) {
    p_AA=0.5*f_AA+0.25*f_Aa; p_Aa=0.5*f_AA+0.5*f_Aa+0.5*f_aa; p_aa=0.25*f_Aa+0.5*f_aa
  } else if(Xi[id]==1) {
    p_AA=0; p_Aa=f_AA+0.5*f_Aa; p_aa=0.5*f_Aa+f_aa
  } else {
    message(cat("error at", id))
  }
  return(c(p_AA, p_Aa, p_aa))
}

#set coefficients
J = -0.15 #= beta_2
h = 0.0001 #= beta_1

#set a field
rect = 100
N = rect*rect
Xi = sample(c(-1,1),N,replace=T)

smap = cbind(rep(1:rect,each=rect),
             rep(1:rect,rect)
)
dmat = as.matrix(dist(smap))

group = rep(1:rect,each=rect)

xixj = c()
for(i in 1:N) xixj = c(xixj, XiXj2(i,group=group))
Ei_t0 = (J*xixj+h*Xi)

#MCMC by Gibbs sampling
for(j in 1:100) {
  Xi_t0 = Xi

  perm = sample(1:N,N)
  for(i in perm) {
    Xi[i] = sample(c(-1,0,1),1,prob=prob(i,group))
    xi = Xi[i]
    if(xi==0) { xi = 1 }
    Ei = (J*XiXj2(i,group=group)+h*xi) #set the energy positive to maximize itself
    
    #Metropolis algorithm
    if(Ei_t0[i]<=Ei) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else if(runif(1,0,1)<exp(Ei-Ei_t0[i])) {
      Xi[i] = Xi[i]; Ei_t0[i] = Ei
    } else {
      Xi[i] = Xi_t0[i]
    }
  }
  E = mean(Ei_t0)
  f = freq(Xi=Xi,rect=rect); barplot(f, main=paste(j,E), col=c("black","grey","white"))
}


svg(file=paste0("IsingJ",J,"_h",h,"splt.svg"), width=5, height=5)
f = freq(Xi=Xi,rect=rect); barplot(f, main=paste(j,E), col=c("black","grey","white"),las=1)
dev.off()

