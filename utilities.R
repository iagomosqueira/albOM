# utilities.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/v2/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(FLCore)

# mc.output {{{

mc.output <- function(x, iter=length(x)) {

  # DIMS
  its <- seq(iter)
  dmns <- list(age=0:14, year=2000:2020, unit=c('F', 'M'), season=1:4, 
    area=1:6, iter=its)
  di <- unlist(lapply(dmns, length))

  # stock.n [15,21,2,4,1,500] - x$N
  stock.n <- FLQuant(NA, dimnames=dmns[-5], units='1000')
  for(i in its) {
    stock.n[,,,,,i] <- c(aperm(x[[i]]$N, c(2,1,4,3)) / 1000)
  }

  # m [15,21,2,4,1,500] - N
  m <- FLQuant(dimnames=dmns[-5])
  for(i in its) {
    m[,,,,,i] <- x[[i]]$M
  }

  # hr - H
  hr <- FLQuant(dimnames=dmns[-c(1, 3)], quant="age", units="hr")
  for(i in its) {
    hr[,,,,,i] <- x[[i]]$H
  }

  # - fisheries

  # catch.sel - sela
  catch.sel <- FLQuant(dimnames=dmns, units="")
  for(i in its)
    catch.sel[,,,,,i] <- aperm(array(x[[i]]$sela, dim=c(15,4,2,6,21)),
      c(1,5,3,2,4))

  # C = N * H * S * W
  landings.n <- expand(stock.n, area=dmns$area) * catch.sel *
    expand(hr, age=dmns$age, unit=dmns$unit)
 
  # CHECK: apply(can * expand(wt(sbio)[, ac(2000:2020)], area=1:6), c(2,6), sum)

  # - indices

  # index.hat - Ihat
  index.hat <- FLQuant(dimnames=c(list(age='all'), dmns[c(2, 4, 6)]))
  for(i in its)
    index.hat[,,,,,i] <- x[[i]]$Ihat

  # index.q - lnq
  lnq <- unlist(lapply(x, '[[', 'lnq'))
  index.q <- FLQuant(exp(lnq), dimnames=dmns[c(4, 6)])

  # - FLPar
  
  # srpars
  B0 <- unlist(lapply(x, '[[', 'B0'))
  R0 <- unlist(lapply(x, '[[', 'R0')) / 1000
  h <- unlist(lapply(x, '[[', 'h'))
  srpars <-  FLPar(v=B0, R0=R0, s=h)
  
  # refpts
  cmsy <- unlist(lapply(x, '[[', 'Cmsy'))
  hrmsy <- FLPar(unlist(lapply(x, '[[', 'hmsy')),
    dimnames=list(params='HRMSY', season=1:4, iter=its))
  refpts <- FLPar(B0=B0, MSY=cmsy,
    HRMSY=c(seasonMeans(as(hrmsy, 'FLQuant'))))

  # rho
  rho <- unlist(lapply(x, '[[', 'rho'))

  # - metrics

  list(stock.n=stock.n, m=m, hr=hr, catch.sel=catch.sel, landings.n=landings.n,
    index.q=index.q, index.hat=index.hat,
    srpars=srpars, refpts=refpts, rho=rho 
  )
}
# }}}

# hr & harvest {{{

setMethod("hr", signature(object="FLombf"),
  function(object) {
    if(length(biols(object)) == 1)
      return(attr(object, 'hr')[[1]])
    else
      return(attr(object, 'hr'))
  })

setMethod("harvest", signature(object="FLombf", catch="missing"),
  function(object) {
    if(length(biols(object)) == 1)
      return(attr(object, 'harvest')[[1]])
    else
      return(attr(object, 'harvest'))
  })
# }}}

# metrics {{{

.annual <- list(
  R=function(x) seasonSums(unitSums(rec(biol(x)))),
  B=function(x) unitSums(tsb(biol(x)))[,,,1],
  HR=function(x) areaSums(seasonMeans(hr(x))),
  C=function(x) seasonMeans(unitSums(Reduce('+', catch(fisheries(x)))))
)

.seasonal <- list(
  R=function(x) unitSums(rec(biol(x))),
  B=function(x) unitSums(tsb(biol(x))),
  HR=function(x) areaSums(hr(x)),
  C=function(x) unitSums(Reduce('+', catch(fisheries(x))))
)
# }}}

# FUNCTIONS {{{

logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

# VECTORIZE
aget.sel.age <- function(nf=6,nselg=5,selidx=c(1,2,3,4,5,5),selpars) {

  seltmp <- rep(NA,20)
  sela <- array(dim=c(na,ns,2,nf))
  
  # mula, sdla [a,s,u], 

  almin <- array(pmax(0, mula - sdla * 1.96), dim=dim(mula))
  almax <- array(pmax(0, mula + sdla * 1.96), dim=dim(mula))
  # n, a, s, u
  alref <- aperm(array(c(mapply(function(x,y) seq(x, y, length=20),
    x=almin, y=almax, SIMPLIFY=TRUE)), dim=c(20, 15, 4, 2)),
    c(2,3,4,1))

  adl <- array(dnorm(c(alref), rep(mula, 20), rep(sdla, 20)),
    dim=c(dim(sdla), 20))
  adl <- adl / as.numeric(array(apply(adl, 1:3, sum, simplify=FALSE),
    dim=dim(adl)))

  bb <-   lapply(1:5, function(f) ifelse(lref < selpars[f,1],
    2^((lref - selpars[f,1]) / selpars[f,2] ^ 2),
    2^(-((lref - selpars[f,1]) / selpars[f,3])^2)))

        for(f in 1:nf) {
          fref <- selidx[f]
          for(l in 1:20) seltmp[l] <- ifelse(lref[l] < selpars[fref,1],2^{-((lref[l]-selpars[fref,1])/selpars[fref,2])^2},2^{-((lref[l]-selpars[fref,1])/selpars[fref,3])^2})

          sela[a,s,g,f] <- sum(seltmp*dl)
        }

  return(sela)
}


get.sel.age <- function(nf=6,nselg=5,selidx=c(1,2,3,4,5,5),selpars)
{

  seltmp <- rep(NA,20)
  sela <- array(dim=c(na,ns,2,nf))
  for(g in 1:2) {
    for(s in 1:4) {
      for(a in 1:na) {

        lmin <- max(0,mula[a,s,g]-sdla[a,s,g]*1.96)
        lmax <- mula[a,s,g]+sdla[a,s,g]*1.96
        lref <- seq(lmin,lmax,length=20)
        dl <- dnorm(lref,mula[a,s,g],sdla[a,s,g])
        dl <- dl/sum(dl)

        for(f in 1:nf) {
        
          fref <- selidx[f]
          for(l in 1:20) seltmp[l] <- ifelse(lref[l] < selpars[fref,1],2^{-((lref[l]-selpars[fref,1])/selpars[fref,2])^2},2^{-((lref[l]-selpars[fref,1])/selpars[fref,3])^2})

          sela[a,s,g,f] <- sum(seltmp*dl)
        }
      }
    }
  }

  return(sela)
}

objfn.init <- function(theta,targv,sela) {

  hxinit <- 1 / (1 + exp(-theta))
  resx <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
  as.vector(wta), as.vector(sela), hxinit) 

  cx <- resx$C
  px <- cx/sum(cx)
  tmpv <- c(resx$rho, as.vector(px))
  objv <- logit(tmpv)

  return(fnscale * (sum((objv - targv) ^ 2)))

}

msyfn <- function(H,ph,sela) {

  hx <- H * ph
  resx <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),hx)
    
  return(sum(resx$C))
}
# }}}

# mcmc.abc {{{
mcmc.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar)
  acp <- rep(0,ngibbs)

  # get initial guess discrepancy

  xx <- sim(R0,dep,h,M,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        wtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- parvecold

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc2.abc {{{
mcmc2.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+2)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
    
  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        } 

        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
            
      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      }  

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc2a.abc {{{
mcmc2a.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+2)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    } 

    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  if(length(yfmsy) == 1) {

    sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean)  

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 
       
      if(length(yfmsy) == 1) {

        sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc3.abc {{{
mcmc3.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+3)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmarold,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]  

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmarold,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # conditional posterior for sigmaR

    res.tmp <- sum(0.5*epsrx^2)
    atmp <- alpR+length(epsrx)/2
    btmp <- betR+res.tmp
    sigmarold <- sqrt(1/rgamma(1,atmp,btmp))

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold,sigmarold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc3a.abc {{{
mcmc3a.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+3)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  if(length(yfmsy) == 1) {

    sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmarold,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      if(length(yfmsy) == 1) {

        sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmarold,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # conditional posterior for sigmaR

    res.tmp <- sum(0.5*epsrx^2)
    atmp <- alpR+length(epsrx)/2
    btmp <- betR+res.tmp
    sigmarold <- sqrt(1/rgamma(1,atmp,btmp))

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold,sigmarold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc4.abc {{{
mcmc4.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+2)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    } 

    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue.mcmc <- matrix(nrow=nits,ncol=prod(dim(resq)))

  dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
  dcpue <- sum(dcpue.loo)

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yof,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  if(length(yof) == 1) {

    zof <- max(hmsyrat-1,0)
    dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
    sprior <- sprior+dof

  } else {

    zof <- hmsyrat-1
    zof[zof < 0] <- 0
    dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
    sprior <- sprior+sum(dof) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
      dcpue <- sum(dcpue.loo)

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yof,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean)  

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      if(length(yof) == 1) {

        zof <- max(hmsyrat-1,0)
        dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
        sprior <- sprior+dof

      } else {

        zof <- hmsyrat-1
        zof[zof < 0] <- 0
        dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
        sprior <- sprior+sum(dof) 

      } 
        

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) {

      theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold)
      dcpue.mcmc[(n-burn)/thin,] <- as.vector(dcpue.loo)

    }

  }

  return(list(pars=theta.mcmc,cpuelogl=dcpue.mcmc,acp=acp))
} 
# }}}

# mcmc5.abc {{{
mcmc5.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+3)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue.mcmc <- matrix(nrow=nits,ncol=prod(dim(resq))) 
  dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
  dcpue <- sum(dcpue.loo)

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yof,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  if(length(yof) == 1) {

    zof <- max(hmsyrat-1,0)
    dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
    sprior <- sprior+dof

  } else {

    zof <- hmsyrat-1
    zof[zof < 0] <- 0
    dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
    sprior <- sprior+sum(dof) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmarold,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
      dcpue <- sum(dcpue.loo) 

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]  
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yof,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      if(length(yof) == 1) {

        zof <- max(hmsyrat-1,0)
        dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
        sprior <- sprior+dof

      } else {

        zof <- hmsyrat-1
        zof[zof < 0] <- 0
        dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
        sprior <- sprior+sum(dof) 

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmarold,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # conditional posterior for sigmaR

    res.tmp <- sum(0.5*epsrx^2)
    atmp <- alpR+length(epsrx)/2
    btmp <- betR+res.tmp
    sigmarold <- sqrt(1/rgamma(1,atmp,btmp))

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) {

      theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold,sigmarold)
      dcpue.mcmc[(n-burn)/thin,] <- as.vector(dcpue.loo)
    }

  }

  return(list(pars=theta.mcmc,cpuelogl=dcpue.mcmc,acp=acp))
} 
# }}}

# sim {{{
sim <- function(R0=1e6, dep=0.5, h=0.75, M=0.075, selpars, epsr, dms, pctarg,selidx) {

  na <- dms[1]
  ns <- dms[2]
  nf <- dms[3]
  nselg <- dms[4]

  # SPR ratio at exploited eqm given steepness and depletion
  # (based on derivation: dep = (4*h*rho+h-1)/(5*h-1))
  rhotarg <- (dep*(5*h-1)+1-h)/(4*h)

  # create selectivity-at-age
  sela <- get.sel.age(nf,nselg,selidx,selpars) 
    
  # target vector (rho+pc)
  targv <- logit(c(rhotarg,as.vector(pctarg)))

  # wrapper objective function to solve (minimise)

  # Estimate initial Fs
  
  hinit <- array(0.15*pctarg, dim=c(ns, nf)) 

  # scalar to give objective function some bite

  theta <- logit(as.vector(hinit))

  zz <- optim(theta,objfn.init,targv=targv,sela=sela,method=("L-BFGS-B"),control=list(trace=0))

  hinit <- array(ilogit(zz$par), dim=c(ns, nf))
  resinit <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), as.vector(hinit)) 

  # relative H-split for MSY calcs
  ph <- as.vector(hinit[] / sum(hinit))

  msy <- optimise(msyfn,interval=c(0,0.9),ph=ph,sela=sela,maximum=TRUE)
  Hmsy <- msy$maximum
  Cmsy <- msy$objective
  resmsy <- msypdyn(c(ns,na,nf), srec, R0, h, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), Hmsy * ph)
  
  Bmsy <- resmsy$Bmsy
  spr0 <- resinit$spr0
  B0 <- R0*spr0
  alp <- 4*h/(spr0*(1-h))
  bet <- (5*h-1)/(B0*(1-h))

  Bratio <- Bmsy/B0
  Rratio <- resmsy$Rmsy/R0
  hmsyv <- apply(array(Hmsy*ph,dim=c(ns,nf)),1,sum)

  if(!all.equal(Rratio, (4*h*Bratio)/(h*(5*Bratio-1)+1-Bratio)))
    warning("B-H invariant check - should be same as Rratio")

  # set up initial numbers-at-age for input to population dynamics

  Rinit <- R0*(4*h*dep)/(h*(5*dep-1)+1-dep)
  Ninit <- array(resinit$N,dim=c(na,ns,2))
  Ninit[] <- Ninit[]*Rinit
  nvec <- as.vector(Ninit)

  # expected catch @ hinit
  zinit <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),
    as.vector(wta),as.vector(sela),as.vector(hinit))
  Cinit <- array(zinit$C,dim=c(ns,nf))

  # main population stuff

  cvec <- as.vector(C)

  ## generating predicted LF and CPUE 
    
  # fishery for CPUE generation

  resp2 <- pdynlfcpue(c(ny,ns,na,nbins,nf),srec,R0,h,psi,epsr,spr0,M,
    as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fcpue)

  N <- array(resp2$N,dim=c(ny,na,ns,2))
  S <- array(resp2$S,dim=c(ny,ns))
  H <- array(resp2$H,dim=c(ny,ns,nf))
  LFhat <- array(resp2$LF,dim=c(ny,nbins,ns,nf))
  Ihat <- array(resp2$I,dim=c(ny,ns))

  # predicted LF distro for relevant fisheries

  LFhat <- apply(LFhat,c(2,4),sum)
  phat <- LFhat[,flf]
  phat <- apply(phat,2,function(x){x <- x/sum(x)})

  return(list(N=N,S=S,H=H,LF=phat,I=Ihat,Bmsy=Bmsy,Cmsy=Cmsy,Hmsy=hmsyv,B0=B0))

}
# }}}

# get.mcmc.vars {{{

get.mcmc.vars <- function(parsmat) {

  varlist <- list()                    
  nnits <- dim(parsmat)[1]
  for(nn in 1:nnits) {

    R0x <- exp(mcpars[nn,1])
    depx <- ilogit(mcpars[nn,2])
    epsrx <- mcpars[nn,3:(ny+1)]
    selvx <- exp(mcpars[nn,(ny+2):npar])
    selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
    xx <- sim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

    varlist[[nn]] <- list()
    varlist[[nn]][['Rtot']] <- apply(xx$N[,1,srec,],1,sum)
    varlist[[nn]][['SSB']] <- xx$S[,srec-1]
    varlist[[nn]][['dep']] <- xx$S[,srec-1]/xx$B0
    varlist[[nn]][['dbmsy']] <- xx$S[,srec-1]/xx$Bmsy
    varlist[[nn]][['Cmsy']] <- xx$Cmsy
    hmsy <- xx$Hmsy
    hy <- apply(xx$H[,,],c(1,2),sum)
    hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 
    varlist[[nn]][['hmsyrat']] <- hmsyrat
    varlist[[nn]][['Ihat']] <- xx$I
    varlist[[nn]][['LFhat']] <- xx$LF 

    if(nn %% 100 == 0) cat("Iteration",nn,"of",nnits,"\n")
  }

  return(varlist)

}
# }}}

# get.mcmc2.vars {{{

get.mcmc2.vars <- function(mcpars) {

  nnits <- dim(mcpars)[1]
  varlist <- vector("list", nnits)
  for(nn in 1:nnits) {

    # R0, pars[,1]
    R0x <- exp(mcpars[nn,1])
    # dep, pars[,2]
    depx <- ilogit(mcpars[nn,2])
    epsrx <- mcpars[nn,3:(ny+1)]
    selvx <- exp(mcpars[nn,(ny+2):npar])
    selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
    hx <- mcpars[nn,npar+1]
    Mx <- mcpars[nn,npar+2]
    xx <- sim(R0x,depx,hx,Mx,selparsx,epsrx,dms,pctarg,selidx)

    # N Rtot SSB dep dbmsy Cmsy hmsyrat H Ihat LFhat B0 R0 M h sela
    varlist[[nn]] <- setNames(nm=list("N", "Rtot", "SSB", "dep", "dbmsy",
      "Cmsy", "hmsy", "hmsyrat", "H", "Hy", "Ihat", "LFhat", "B0", "R0",
      "M", "h", "sela"))

    varlist[[nn]][['N']] <- xx$N
    varlist[[nn]][['Rtot']] <- apply(xx$N[,1,srec,],1,sum)
    varlist[[nn]][['SSB']] <- xx$S[,srec-1]
    varlist[[nn]][['dep']] <- xx$S[,srec-1]/xx$B0
    varlist[[nn]][['dbmsy']] <- xx$S[,srec-1]/xx$Bmsy
    varlist[[nn]][['Cmsy']] <- xx$Cmsy
    # HMSY per season[s,f]
    hmsy <- xx$Hmsy
    varlist[[nn]][['hmsy']] <- hmsy
    # ADD over f
    hy <- apply(xx$H[,,],c(1,2),sum)
    varlist[[nn]][['hy']] <- hy
    hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 
    # Annual hrmsy
    varlist[[nn]][['hmsyrat']] <- hmsyrat 
    varlist[[nn]][['H']] <- xx$H
    varlist[[nn]][['Ihat']] <- xx$I
    resq <- log(I[,,fcpue] / xx$I)
    varlist[[nn]][['lnq']] <- apply(resq, 2, mean)
    varlist[[nn]][['LFhat']] <- xx$LF 
    varlist[[nn]][['B0']] <- xx$B0
    varlist[[nn]][['R0']] <- R0x
    varlist[[nn]][['M']] <- Mx
    varlist[[nn]][['h']] <- hx
    varlist[[nn]][['rho']] <- c(rho(FLQuant(epsrx)))
    varlist[[nn]][['sela']] <- get.sel.age(nf,nselg,selidx,selpars) 

    if(nn %% 100 == 0) cat("Iteration",nn,"of",nnits,"\n")
  }

  return(varlist)

}
# }}}

# plot.mcmc.vars {{{
plot.mcmc.vars <- function(varlist,ptype) {

  nnits <- length(varlist)

  if(ptype == 'dep') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dep
    vmin <- 0
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmax <- max(vq[3,])
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='SSB depletion',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

  if(ptype == 'bmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dbmsy
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmin <- 0 
    vmax <- max(vq) 
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='Bmsy ratio',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

  if(ptype == 'hmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$hmsyrat
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmin <- 0
    vmax <- max(vq) 
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab=expression(H[y]/H[msy]),col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  } 

  if(ptype == 'rec') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$Rtot
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmin <- 0
    vmax <- max(vq) 
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='Recruitment',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

}
# }}}

# {{{ plot.mcmc.cpue
plot.mcmc.cpue <- function(varlist) {

  nnits <- length(varlist)

  vv <- array(dim=c(nnits,ny,ns))
   for(nn in 1:nnits) {
     
     tmpv <- varlist[[nn]]$Ihat
     iobs <- I[,,fcpue]
     if(seasonq) {

       if(qtrend) {

        resq <- log(iobs/(tmpv*qt))

       } else {

         resq <- log(iobs/tmpv) 

       }

       lnq <- apply(resq,2,mean)
       if(qtrend) {

         vv[nn,,] <- t(apply(tmpv*qt,1,function(x,lnq){x <- x*exp(lnq)},lnq))*rlnorm(ny*ns,0,sdcpue)

       } else {

         vv[nn,,] <- t(apply(tmpv,1,function(x,lnq){x <- x*exp(lnq)},lnq))*rlnorm(ny*ns,0,sdcpue) 

       }

     } else {

       if(qtrend) {

        resq <- log(iobs/(tmpv*qt))

       } else {

         resq <- log(iobs/tmpv) 

       } 

       lnq <- mean(resq)
       if(qtrend) {

         vv[nn,,] <- tmpv*qt*exp(lnq)*rlnorm(ny*ns,0,sdcpue)  

       } else {

         vv[nn,,] <- tmpv*exp(lnq)*rlnorm(ny*ns,0,sdcpue) 

       }

     }
   }   
   
   vq <- apply(vv,c(2,3),quantile,c(0.025,0.5,0.975))

   # ggplot the sumbitch

   vdf <- expand.grid(year=yrs,season=1:ns,obs=NA,hat=NA,lq=NA,uq=NA)
   vdf$obs <- as.vector(iobs)
   vdf$hat <- as.vector(vq[2,,])
   vdf$lq <- as.vector(vq[1,,]) 
   vdf$uq <- as.vector(vq[3,,]) 
   ggplot(vdf)+geom_line(aes(x=year,y=hat),colour='blue')+geom_line(aes(x=year,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=year,y=uq),colour='blue',linetype='dashed')+geom_point(aes(x=year,y=obs),colour='magenta')+facet_wrap(~season)+ylab("CPUE")+theme_bw()

}
# }}}

# {{{ plot.mcmc.lf
plot.mcmc.lf <- function(varlist) {

  nnits <- length(varlist)
  vv <- array(dim=c(nnits,nbins,nselg)) 
    for(nn in 1:nnits) vv[nn,,] <- varlist[[nn]]$LF
    vq <- apply(vv,c(2,3),quantile,c(0.025,0.5,0.975))
    vdf <- expand.grid(length=mulbins,fishery=1:nselg,obs=NA,hat=NA,lq=NA,uq=NA) 
    vdf$obs <- as.vector(pobs)
    vdf$hat <- as.vector(vq[2,,])
    vdf$lq <- as.vector(vq[1,,]) 
    vdf$uq <- as.vector(vq[3,,]) 
    ggplot(vdf)+geom_line(aes(x=length,y=hat),colour='blue')+geom_line(aes(x=length,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=length,y=uq),colour='blue',linetype='dashed')+geom_point(aes(x=length,y=obs),colour='magenta')+facet_wrap(~fishery)+ylab("Length frequency")+theme_bw()

}
# }}}

# {{{ plot.mcmc.sel 
plot.mcmc.sel <- function(mcpars) {

  nnits <- dim(mcpars)[1]
  selparsx <- exp(mcpars[,paridx[[3]]])
  msel <- array(dim=c(nnits,nbins,nselg))

  for(nn in 1:nnits) {
    for(f in 1:nselg) {
      
      sx50 <- selparsx[nn,1+(f-1)]
      sxL <- selparsx[nn,nselg+f]
      sxR <- selparsx[nn,2*nselg+f]

      for(l in 1:nbins) {
  
        lref <- mulbins[l]
        msel[nn,l,f] <- ifelse(lref<sx50,2^{-(lref-sx50)^2/(sxL^2)},2^{-(lref-sx50)^2/(sxR^2)})
        
      }
    }
  }

  vq <- apply(msel,c(2,3),quantile,c(0.025,0.5,0.975))
  vdf <- expand.grid(length=mulbins,fishery=1:nselg,med=NA,lq=NA,uq=NA) 
  vdf$med <- as.vector(vq[2,,])
  vdf$lq <- as.vector(vq[1,,])
  vdf$uq <- as.vector(vq[3,,]) 
  ggplot(vdf)+geom_line(aes(x=length,y=med),colour='blue')+geom_line(aes(x=length,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=length,y=uq),colour='blue',linetype='dashed')+facet_wrap(~fishery)+ylab("Size selectivity")+theme_bw()

}

# }}}

# get.mcmc.vars {{{

get.mcmc.vars <- function(varlist,vtype) {

  nnits <- length(varlist)

  if(vtype == 'dep') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dep
    vmin <- 0
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  if(vtype == 'bmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dbmsy
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  if(vtype == 'hmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$hmsyrat
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  if(vtype == 'rec') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$Rtot
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  return(round(vq,2))

} # }}}

# cpue.mase {{{

get.cpue.mase <- function(k) {
  
  iobs <- I[,,fcpue]
  xhat <- mcvars[[k]]$Ihat
  resq <- log(iobs/xhat) 
  if(seasonq) {

    lnq <- apply(resq,2,mean) 
    ihat <- t(apply(xhat,1,function(x,lnq){x <- x*exp(lnq)},lnq))

  } else {

    lnq <- mean(resq)
    ihat <- iobs*exp(lnq)

  }

  # naive forecast error by season

  idiff <- log(iobs[-1,]/iobs[-dim(iobs)[1],])

  # model forecast error

  ediff <- log(iobs[-1,]/ihat[-1,])
  mudiff <- apply(abs(idiff),2,mean)
  zdiff <- apply(abs(ediff),1,function(x,mudiff){x <- x/mudiff},mudiff)
  mases <- apply(zdiff,1,mean)

  return(mases)

}
# }}}
