# data.R - LOAD SS3 model and add ABC output
# abc_tuna/om/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

source("utilities.R")

# source("output_base.R")

# LOAD SS3 base run

load("output/base.rda")

fy <- 2045


# --- LOAD abc run

run <- mget(load("model/abc_run5b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out <- mc.output(run$mcvars)
out <- mc.output(run$mcvars, 50)

# CHECK

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

# FLCatch(es), RESCALE landings.n & areaSums to drop area names
cas <- Map(function(x, y)
  FLCatch(landings.n=areaSums(x), landings.wt=wt(bio),
  catch.sel=areaMeans(y), discards.n=x %=% 0, discards.wt=wt(bio)), 
  x=divide(out$landings.n, 5),
  y=divide(out$catch.sel %*% (out$landings.n %=% 1), 5))

# FLFisheries
fis <- FLFisheries(lapply(cas, function(x)
  FLFishery(effort=unitSums(catch(cas[[1]])) %=% 0, ALB=x)))

names(fis) <- c(paste0("LL", 1:4), "PS", "Other")

om <- FLombf(biols=FLBiols(ALB=bio), fisheries=fis,
  refpts=FLPars(ALB=out$refpts))

# EXTEND
om <- fwdWindow(om, end=fy)

# ADD hr [1, y, 1, s, f, it]
attr(om, 'hr') <- FLQuants(ALB=window(out$hr, end=2045))

attr(om, 'harvest') <- FLQuants(ALB=Reduce(abind, Map(function(i, j) i %*% j,
  i=divide(hr(om), 5),
  j=lapply(fisheries(om), function(y) catch.sel(y[[1]])))))
units(attr(om, 'harvest')$ALB) <- "hr"

# BUG: show(FLQuant) slow because of mean/MAD calculation

# SRR deviances
deviances(om) <- rlnormar1(n=500, meanlog=0, sdlog=run$sigrpost, rho=out$rho,
  years=2020:2045)

# PLOT
plot(metrics(om, metrics=.annual))

# - BUILD oem

# idx: FLIndexBiomass by season, with sel.pattern by sex

sp <- expand(out$catch.sel[,,,,1], year=2000:fy)
dimnames(sp) <- list(area='unique')

NW <- FLIndexBiomass(
  index=window(out$index.hat %*% out$index.q, end=fy),
  index.q=expand(out$index.q, year=2000:fy),
  catch.wt=expand(unitMeans(wt(biol(om))), year=2000:fy),
  sel.pattern=sp)

# stk: no units
oem <- FLoem(observations=list(ALB=list(idx=FLIndices(NW=NW),
  stk=simplify(stock(om)[[1]], 'unit'))),
  method=sampling.oem)

# SAVE
save(om, oem, file='output/om5b-test50.rda', compress='xz')
