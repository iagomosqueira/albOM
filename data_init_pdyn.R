######################################################
# initial popn for IOTC ALB ABC ######################
######################################################
# R. Hillary & I. Mosqueira, 2023 ####################
######################################################

library(Rcpp)

sourceCpp("utilities/init_pdyn.cpp")
logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

# biological and fishery guff

R0 <- 1e+6
del <- 0.5 # eqm female SSB depletion
hh <- 0.75 # steepness
ns <- 4 # number of seasons
na <- 20 # number of ages
psi <- 0.5 # fraction of females at birth
# growth-by-sex
k <- c(0.2,0.2)
linf <- c(100,100)
l0 <- c(10,10)
sdla <- c(0.1,0.1)
# maturity-at-length
ml50 <- 60
ml95 <- 80
# number of fisheries
nf <- 2
# selectivity-at-length
sl50 <- c(45,45)
sl95 <- c(65,65)
# natural mortality (by season)
M <- 0.1
# weight-at-length
alw <- 2e-6
blw <- 3
# eqm catch by season and fishery
C <- array(dim=c(ns,nf))
C[] <- 10000

######################################
# generate biology and fleet objects #
######################################

# mean length-at-age

mula <- array(dim=c(na,ns,2))
ages <- 1:na
for(g in 1:2) {
  for(s in 1:ns) {

    aa <- ages+(s-1)*0.25
    mula[,s,g] <- l0[g]+(linf[g]-l0[g])*(1-exp(-k[g]*aa))

  }
}

# maturity-at-age, weight-at-age and selectivity-at-age

nl <- 20
mata <- wta <- array(dim=c(na,ns,2))
sela <- array(dim=c(na,ns,2,nf))

for(g in 1:2) {
  for(s in 1:4) {
    for(a in 1:na) {

      lmin <- max(0,mula[a,s,g]*(1-sdla[g]*1.96))
      lmax <- mula[a,s,g]*(1+sdla[g]*1.96)
      lref <- seq(lmin,lmax,length=nl)
      dl <- dlnorm(lref,log(mula[a,s,g]),sdla[g])
      dl <- dl/sum(dl)
      mlref <- 1/(1+19^{-(lref-ml50)/(ml95-ml50)})
      wlref <- alw*lref^blw
      mata[a,s,g] <- sum(mlref*dl)
      wta[a,s,g] <- sum(wlref,dl)

      for(f in 1:nf) {
      
        slref <- 1/(1+19^{-(lref-sl50[f])/(sl95[f]-sl50[f])})
        sela[a,s,g,f] <- sum(slref*dl)

      }
    }
  }
}

# SPR ratio at exploited eqm given steepness and depletion
# (based on derivation: del = (4*hh*rho+hh-1)/(5*hh-1))

rhotarg <- (del*(5*hh-1)+1-hh)/(4*hh)

# catch fraction by season 

pctarg <- C/sum(C)

#####################################
# OK now we can estimate initial Fs #
#####################################

hinit <- array(0.025,dim=c(ns,nf))
srec <- 2
psi <- 0.5
res <- initpdyn(c(ns,na,nf),srec,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),as.vector(hinit))

# target vector (rho+pc)

targv <- logit(c(rhotarg,pctarg))

# wrapper objective function to solve (minimise)

fnscale <- 1 # scalar to give objective function some bite

objfn.init <- function(theta) {

  hxinit <- 1/(1+exp(-theta))
  resx <- initpdyn(c(ns,na,nf),srec,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),hxinit) 

  px <- resx$C/sum(resx$C)
  tmpv <- c(resx$rho,as.vector(px))
  objv <- logit(tmpv)

  return(fnscale*(sum((objv-targv)^2)))

}

theta <- logit(as.vector(hinit))
system.time(zz <- optim(theta,objfn.init,method=("L-BFGS-B"),control=list(trace=1)))

hinit <- array(ilogit(zz$par),dim=c(ns,nf))
resinit <- initpdyn(c(ns,na,nf),srec,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),as.vector(hinit)) 

######################
# MSY estimation bit #
######################

sourceCpp("utilities/msy_pdyn.cpp")

res <- msypdyn(c(ns,na,nf),srec,R0,hh,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),hinit)

# relative H-split for MSY calcs

ph <- as.vector(hinit[]/sum(hinit))

# MSY wrapper

msyfn <- function(H) {

  hx <- H*ph
  resx <- msypdyn(c(ns,na,nf),srec,R0,hh,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),hx)

  return(sum(resx$C))

}

msy <- optimise(msyfn,interval=c(0,0.5),maximum=TRUE)
Hmsy <- msy$maximum
Cmsy <- msy$objective
resmsy <- msypdyn(c(ns,na,nf),srec,R0,hh,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),Hmsy*ph)
Bmsy <- resmsy$Bmsy
spr0 <- resinit$spr0
B0 <- R0*spr0
alp <- 4*hh/(spr0*(1-hh))
bet <- (5*hh-1)/(B0*(1-hh))
Bratio <- Bmsy/B0
Bratio
Rratio <- resmsy$Rmsy/R0
Rratio
(4*hh*Bratio)/(hh*(5*Bratio-1)+1-Bratio) # B-H invariant check - should be same as Rratio

# set up initial numbers-at-age for input to population dynamics

Rinit <- R0*(4*hh*del)/(hh*(5*del-1)+1-del)
Ninit <- array(resinit$N,dim=c(na,ns,2))
Ninit[] <- Ninit[]*Rinit
nvec <- as.vector(Ninit)

# expected catch @ hinit

zinit <- msypdyn(c(ns,na,nf),srec,R0,hh,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),as.vector(hinit))
Cinit <- array(zinit$C,dim=c(ns,nf))

# main population stuff

sourceCpp("utilities/pdyn.cpp")

ny <- 10
epsr <- rep(0,ny-1)
Cb <- array(dim=c(ny,ns,nf))
for(y in 1:ny) Cb[y,,] <- Cinit
cvec <- as.vector(Cb)

resp <- pdyn(c(ny,ns,na,nf),srec,R0,hh,psi,epsr,spr0,M,as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec)

N <- array(resp$N,dim=c(ny,na,ns,2))
S <- array(resp$S,dim=c(ny,ns))
H <- array(resp$H,dim=c(ny,ns,nf))

####################################
# generating predicted LF and CPUE #
####################################

sourceCpp("utilities/pdyn_lfcpue.cpp")

# p(l | a, s, g)

lbins <- seq(0,120,by=6)
nbins <- length(lbins)-1
mulbins <- 0.5*(lbins[-1]+lbins[-length(lbins)])

pla <- array(dim=c(nbins,na,ns,2))

for(g in 1:2) {
  for(s in 1:ns) {
    for(a in 1:na) {

      dx <- dnorm(log(mulbins),log(mula[a,s,g]),sdla[g],FALSE)
      dx <- dx/sum(dx)
      pla[,a,s,g] <- dx 

    }
  }
}

# fishery for CPUE generation

fref <- 1

resp2 <- pdynlfcpue(c(ny,ns,na,nl,nf),srec,R0,hh,psi,epsr,spr0,M,as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fref)

N <- array(resp2$N,dim=c(ny,na,ns,2))
S <- array(resp2$S,dim=c(ny,ns))
H <- array(resp2$H,dim=c(ny,ns,nf))
LF <- array(resp2$LF,dim=c(ny,nbins,ns,nf))
I <- array(resp2$I,dim=c(ny,ns))

# exploring "equivalence" between K-L and multinomial ESS

nits <- 1000
neff <- seq(50,200,by=10)
KL <- array(dim=c(length(neff),nits))
ptrue <- LF[1,,1,1]
for(i in 1:length(neff)) {

  nsim <- rmultinom(nits,neff[i],ptrue)
  psim <- apply(nsim,2,function(x){x <- x/sum(x)})

  for(n in 1:nits) {

    xobs <- psim[,n]
    xobs[xobs == 0] <- 1e-6
    xtmp <- xobs*log(xobs/ptrue)
    KL[i,n] <- sum(xtmp[!is.nan(xtmp)])

  }
}

klsumm <- apply(KL,1,quantile,c(0.025,0.5,0.975))
kmax <- max(klsumm)
kmin <- 0
plot(neff,klsumm[2,],ylim=c(kmin,kmax),type='l',xlab='ESS',ylab="K-L divergence")
lines(neff,klsumm[1,],lty=2)
lines(neff,klsumm[3,],lty=2)

# save it

save.image("data/init_pdyn.rda")

