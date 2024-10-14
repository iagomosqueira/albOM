# summarise the OM outputs

library(icesTAF)
mkdir("report")

library(ggplot2)
library(patchwork)

source("utilities.R")

# abc_run4 {{{

load("model/abc_run4.rda")
mkdir("report/run4")
qtrend <- FALSE

# plots
taf.png("report/run4/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run4/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run4/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run4/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run4/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run4a {{{

load("model/abc_run4a.rda")
mkdir("report/run4a")

qtrend <- FALSE

# plots
taf.png("report/run4a/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run4a/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run4a/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run4a/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run4a/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run4b {{{

load("model/abc_run4b.rda")
mkdir("report/run4b")

qtrend <- FALSE

# plots
taf.png("report/run4b/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run4b/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run4b/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run4b/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run4b/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run5 {{{

load("model/abc_run5.rda")
mkdir("report/run5")
qtrend <- FALSE

# plots
taf.png("report/run5/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run5/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run5/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run5/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run5/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run5a {{{

load("model/abc_run5a.rda")
mkdir("report/run5a")

qtrend <- FALSE

# plots
taf.png("report/run5a/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run5a/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run5a/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run5a/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run5a/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run5b {{{

load("model/abc_run5b.rda")
mkdir("report/run5b")

qtrend <- FALSE

# plots
taf.png("report/run5b/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run5b/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run5b/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run5b/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run5b/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run6 {{{

load("model/abc_run6.rda")
mkdir("report/run6")

# plots
taf.png("report/run6/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run6/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run6/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run6/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run6/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}

# abc_run6a {{{

load("model/abc_run6a.rda")
mkdir("report/run6a")

# plots
taf.png("report/run6a/timeseries.png")
par(mfrow = c(2,2))
plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
dev.off()

taf.png("report/run6a/cpue.png")
plot.mcmc.cpue(mcvars)
dev.off()

taf.png("report/run6a/lf.png")
plot.mcmc.lf(mcvars)
dev.off()

taf.png("report/run6a/sel.png")
plot.mcmc.sel(mcpars)
dev.off()

taf.png("report/run6a/mase.png")
ggplot(melt(cpuemase, measure.vars=c('LL1', 'LL2', 'LL3','LL4')),
  aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  geom_hline(yintercept=1, linetype=2) +
  xlab("CPUE index") + ylab("MASE") +
  theme(legend.position="none")
dev.off()

# }}}






# sigmaR stuff

sigrpost <- mcpars[,npar+3]
dpost <- density(sigrpost)
dprior <- density(psigmaR)
p1 <- dpost$y
p1 <- p1/sum(p1)
x1 <- dpost$x
p2 <- dprior$y
p2 <- p2/sum(p2)
x2 <- dprior$x
pmax <- max(c(max(p1),max(p2)))
xmin <- min(c(min(x1),min(x2)))
xmax <- max(c(max(x1),max(x2)))
plot(x1,p1,type='l',xlim=c(xmin,xmax),ylim=c(0,pmax),xlab=expression(sigma[R]),ylab='density',main='OM scenario R2a')
lines(x2,p2,lty=2,col='purple')
legend("topright",lty=c(1,2),col=c("black","purple"),legend=c("Posterior","Prior"),bty='n')

round(quantile(psigmaR,c(0.025,0.5,0.975)),2)
round(quantile(sigrpost,c(0.025,0.5,0.975)),2)

