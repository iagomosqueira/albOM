# test_om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albOM/test_om.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: verify(oem, om)


library(mse)

load("output/om5b-test50.rda")

pcbar <- as.matrix(fread('boot/data/pcbar.dat'))
load('boot/initial/data/pla.rda')

fy <- 2045

# C = 0

ctrl <- fwdControl(year=2021:fy, quant="catch", value=1e-8)

system.time(
tesf0 <- fwdabc.om(om, ctrl, pcbar=pcbar, pla=pla)$om
)

plot(metrics(tesf0, metrics=.annual)) +
  geom_vline(xintercept=2021, linetype=2, alpha=0.5)

# C

ctrl <- fwdControl(year=2021:fy, quant="catch", value=3500)
ctrl <- fwdControl(year=2021:fy, quant="catch", value=35000)

system.time(
tesmsy <- fwdabc.om(om, ctrl, pcbar=pcbar, pla=pla)$om
)

plot(metrics(tesmsy, metrics=.annual)) +
  geom_vline(xintercept=2021, linetype=2, alpha=0.5)

metrics(tesmsy, metrics=.annual)

# Cs

tescs <- parallel::mclapply(seq(3e4, 5e4, length=3), function(x) {
  fwdabc.om(om, fwdControl(year=2021:fy, quant="catch", value=x),
    pcbar=pcbar, pla=pla)$om
  }, mc.cores=3)

taf.png('~/catch_fwd.png', width=3200)
(plot(metrics(tescs[[1]], metrics=.annual)) + ggtitle("C=30000")) +
(plot(metrics(tescs[[2]], metrics=.annual)) + ggtitle("C=40000")) +
(plot(metrics(tescs[[3]], metrics=.annual)) + ggtitle("C=50000"))
dev.off()

plot(metrics(tescs[[3]], metrics=.seasonal)) + ggtitle("C=30000")

plot(hr(tescs[[2]]))
