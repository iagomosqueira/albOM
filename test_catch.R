# test_catch.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albOM/test_catch.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

load("output/base.rda")

run <- mget(load("model/abc_run5b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# --- COMPARE input and SS3 catch ~ computed from mcvars

# iter 1
i <- 2
aa <- run$mcvars[[i]]

# C = N * hr * sel
dim(aa$N)
dim(aa$H)
dim(aa$sela)

# [a,y,u,s,(f)] [15,21,2,4,6]
nn <- array(aperm(aa$N, c(2,1,4,3)), dim=c(15,21,2,4,6))
hh <- aperm(array(aa$H, dim=c(21,4,6,15,2)), c(4,1,5,2,3))
se <- aperm(array(aa$sela, dim=c(15,4,2,6,21)), c(1,5,3,2,4))
wa <- aperm(array(run$wta, dim=c(15,4,2,21,6)), c(1,4,3,2,5))

ca <- nn * hh * se * wa

# INPUT to abc
apply(run$C, 1, sum)

# FROM mcvars output
setNames(apply(ca, 2, sum), nm=2000:2020)

foo <- function(i) {
  aa <- run$mcvars[[i]]
  nn <- array(aperm(aa$N, c(2,1,4,3)), dim=c(15,21,2,4,6))
  hh <- aperm(array(aa$H, dim=c(21,4,6,15,2)), c(4,1,5,2,3))
  se <- aperm(array(aa$sela, dim=c(15,4,2,6,21)), c(1,5,3,2,4))
  wa <- aperm(array(run$wta, dim=c(15,4,2,21,6)), c(1,4,3,2,5))
  # C = N * H * S * W
  ca <- nn * hh * se * wa

  setNames(apply(ca, 2, sum), nm=2000:2020)
}

cab <- lapply(seq(500), foo)
caf  <- divide(apply(Reduce('+', lapply(cas, catch)), c(2,6), sum), 6)

seasonSums(unitSums(iter(catch(om)[[1]], 1)))
cab[[1]]


x <- run$mcvars

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
