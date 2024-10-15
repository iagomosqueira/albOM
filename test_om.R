# test_om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albOM/test_om.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: verify(oem, om)


library(mse)

load("output/om6b.rda")

fy <- 2045

# C = 0

ctrl <- fwdControl(year=2021:fy, quant="catch", value=1200)

system.time(
tesf0 <- fwdabc.om(iter(om, 1:50), ctrl, pcbar=pcbar, pla=pla)$om
)

plot(metrics(tesf0, metrics=.annual)) +
  geom_vline(xintercept=2021, linetype=2, alpha=0.5)

# C = MSY

ctrl <- fwdControl(year=2021:fy, quant="catch", value=refpts(om)$MSY)

system.time(
tesmsy <- fwdabc.om(om, ctrl, pcbar=pcbar, pla=pla)$om
)

plot(metrics(tesmsy$om, metrics=.annual)) +
  geom_vline(xintercept=2021, linetype=2, alpha=0.5)
