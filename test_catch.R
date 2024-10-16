# test_catch.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albOM/test_catch.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLCore)

load("output/base.rda")

run <- mget(load("model/abc_run5b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# --- COMPARE input and SS3 catch ~ computed from mcvars

# iter
i <- 2
aa <- run$mcvars[[i]]

# [a,y,u,s,(f)] [15,21,2,4,6]

# N [21,15,4,2] - nn [15,21,2,4,(6)]
nn <- array(aperm(aa$N, c(2,1,4,3)), dim=c(15,21,2,4,6))

# H [21,4,6] - hh [(15),21,(2),4,6]
hh <- aperm(array(aa$H, dim=c(21,4,6,15,2)), c(4,1,5,2,3))

# sela [15,4,2,6] - se [(15),21,2,4,6]
se <- aperm(array(aa$sela, dim=c(15,4,2,6,21)), c(1,5,3,2,4))

# wta [15,4,2] - wa [15,(21),2,4,(6)]
wa <- aperm(array(run$wta, dim=c(15,4,2,21,6)), c(1,4,3,2,5))

# CA
ca <- nn * hh * se * wa

# FROM mcvars output
setNames(apply(ca, 2, sum), nm=2000:2020)

# INPUT to abc
apply(run$C, 1, sum)


# RUN along iters
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

# 50 iters
cab <- lapply(seq(50), foo)

# LOADED om
library(mse)

load('output/om5b-test50.rda')

omca <- apply(catch(om)$ALB, c(2,6), sum)[, ac(2000:2020)]

cab[[1]]
iter(omca, 1)


