# data.R - LOAD SS3 model and add ABC output
# abc_tuna/om/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

source("utilities.R")

source("output_base.R")

# LOAD SS3 base run

load("output/base.rda")

fy <- 2045


# --- LOAD abc run

run <- mget(load("model/abc_run6.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out <- mc.output(run$mcvars, run$C)

#  CREATE om

its <- dims(out$m)$iter

# - FLBiol (stock.n, m)

bio <- propagate(window(sbio, start=2000), its)

n(bio) <- out$stock.n
m(bio) <- out$m

sr(bio) <- predictModel(model=bevholtss3()$model, params=out$srpars)

# spwn, mid Q4
spwn(bio)[,,,4] <- 0.5

# FIX mat BUG: CHECK readFLSss3 +ss3om 
mat(bio)[,,'F',1:4] <- mat(bio)[,,'F',1]

# - FLFisheries (catch.n, catch.sel)

# FLCatch(es)
cas <- Map(function(x, y) FLCatch(landings.n=x, landings.wt=wt(bio),
  catch.sel=y, discards.n=x %=% 0, discards.wt=wt(bio)),
  x=divide(out$catch.n, 5), y=divide(out$catch.sel %*% (out$catch.n %=% 1), 5))

# FLFisheries
fis <- FLFisheries(lapply(cas, function(x)
  FLFishery(effort=unitSums(catch(cas[[1]])) %=% 0, ALB=x)))

names(fis) <- c(paste0("LL", 1:4), "PS", "Other")

om <- FLombf(biols=FLBiols(ALB=bio), fisheries=fis,
  refpts=FLPars(ALB=out$refpts))

# EXTEND

om <- fwdWindow(om, end=fy)

# ADD hr & hrbar
attr(om, 'hr') <- FLQuants(ALB=expand(n(biol(om)), area=1:6) %=% as.numeric(NA))
attr(om, 'hrbar') <- FLQuants(ALB=expand(quantSums(n(biol(om))), area=1:6)
  %=% as.numeric(NA))

# SRR deviances

deviances(om) <- rlnormar1(n=500, meanlog=0, sdlog=0.355, rho=out$rho,
  years=2020:2045)

# - BUILD oem

# idx: FLIndexBiomass by season, with sel.pattern by sex

NW <- FLIndexBiomass(
  index=window(out$index.hat %*% out$index.q, end=fy),
  index.q=expand(out$index.q, year=2000:fy),
  sel.pattern=expand(out$sel, year=2000:fy),
  catch.wt=wt(biol(om)),
  range=c(startf=0.5, endf=0.5))

# stk: no units
oem <- FLoem(observations=list(ALB=list(idx=FLIndices(NW=NW),
  stk=simplify(stock(om)[[1]], 'unit'))),
  method=sampling.oem)

# TODO: verify(oem, om)


# -- SAVE
save(om, oem, file='output/om6b.rda', compress='xz')

# --- TESTS

# - C = 0

ctrl <- fwdControl(year=2021:fy, quant="catch", value=0)

system.time(
  tes2 <- fwdabc.om(iter(om, 1:10), ctrl, pcbar=pcbar, pla=pla)
)

system.time(
tesf0 <- fwdabc.om(om, ctrl, pcbar=pcbar, pla=pla)
)

mets <- list(rec=function(x) unitSums(rec(x)), B=function(x) unitSums(tsb(x)))
plot(biol(tes2$om), metrics=mets) + geom_vline(xintercept=ISOdate(2011,1,1))

mets <- list(rec=function(x) seasonSums(unitSums(rec(x))),
  B=function(x) unitSums(tsb(x))[,,,1])
plot(biol(tes2$om), metrics=mets) + geom_vline(xintercept=2021)


