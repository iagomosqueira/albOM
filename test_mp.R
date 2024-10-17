# test_om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albOM/test_om.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: verify(oem, om)


library(mse)

load("output/om5b.rda")

it <- 1

om <- iter(om, it)
attr(om, 'harvest')$ALB <- iter(attr(om, 'harvest')$ALB, it)
attr(om, 'hr')$ALB <- iter(attr(om, 'hr')$ALB, it)
oem <- iter(oem, it)

pcbar <- as.matrix(fread('boot/data/pcbar.dat'))
load('boot/initial/data/pla.rda')

fy <- 2045

# TEST fwdabc.om {{{

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
ctrl <- fwdControl(year=2021:fy, quant="catch",
  value=rnorm(100, 35000, 10000))

system.time(
tesmsy <- fwdabc.om(om, ctrl, pcbar=pcbar, pla=pla)$om
)

plot(metrics(tesmsy, metrics=.annual)) +
  geom_vline(xintercept=2021, linetype=2, alpha=0.5) +
  ylim(0, NA)

metrics(tesmsy, metrics=.annual)$C
metrics(tesmsy, metrics=.annual)$HR / refts(om)

for(y in 2021:fy)
  tesmsy <- fwdabc.om(tesmsy, fwdControl(year=y, quant="catch", value=3500), 
    pcbar=pcbar, pla=pla)$om




# TODO: SSB

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

# }}}

# - TEST MP {{{

source('../albMSE/utilities.R', echo=TRUE)

trace("goFish", browser, exit=browser, signature = c("FLombf"))
untrace("goFish", signature = c("FLombf"))

projection(om) <- mseCtrl(method=fwdabc.om, args=list(pcbar=pcbar, pla=pla))

load_all('~/Projects/FLR/code/mse/mse')
method(oem) <- sampling.oem

# ctrl

mpctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=perfect.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=rep(23000, dims(om)$iter)))))

# mp
tes <- mp(om, oem=oem, ctrl=mpctrl, args=list(iy=2020, fy=2030, frq=1))

index(observations(oem(tes))$ALB$idx[[1]])[, ac(2020:2029)]

plot(metrics(om(tes), metrics=.annual)) +
  ylim(0, NA)

metrics(om(tes), metrics=.annual)$C
metrics(om(tes), metrics=.annual)$B


# fwd @ catch=3500
ctrl <- fwdControl(year=2021, quant="catch", value=23000)
tfwd <- fwdabc.om(om, ctrl, pcbar=pcbar, pla=pla)$om

plot(index(survey(stock(tfwd)[[1]], observations(oem)$ALB$idx[[1]])))

plot(metrics(tfwd, metrics=.annual)) +
  geom_vline(xintercept=2021, linetype=2, alpha=0.5) +
  ylim(0, NA)

# cpue MP

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind),
  # HCR
  hcr = mseCtrl(method=cpue.hcr, args=list(target=1))))

tes <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2020, fy=2025))

plot(om, tes)



