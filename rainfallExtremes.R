# Calculate trends in rainfall extremes
## January 2020

# functions
# Set functions
AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}


## import
dat <- read.table("monthlyprecipdata.csv", header=T, sep=",")

## calculate yearly sums
precip.yr.sum <- xtabs(dat$Monthly.Precipitation.Total..millimetres. ~ dat$Year)
precip.yr.sum <- precip.yr.sum[-length(precip.yr.sum)]
year.vec <- as.numeric(names(precip.yr.sum))
plot(year.vec, as.numeric(precip.yr.sum), type="l", pch=19, xlab="year", ylab="annual precipitation (mm)")
fit.yr <- lm(precip.yr.sum ~ year.vec)
abline(fit.yr, lty=2, lwd=2, col="red")
abline(h=mean(as.numeric(precip.yr.sum)),lty=2, lwd=3)
linreg.ER(year.vec, as.numeric(precip.yr.sum))

## calculate seasonal variation
year.vec <- as.numeric(unique(dat$Year))
precip.sd <- precip.cv <- precip.mn <- rep(0,length(year.vec))
for (y in 1:length(year.vec)) {
  dat.sub <- subset(dat, Year == year.vec[y])
  precip.sd[y] <- sd(dat.sub$Monthly.Precipitation.Total..millimetres.)
  precip.mn[y] <-  mean(dat.sub$Monthly.Precipitation.Total..millimetres.)
  precip.cv[y] <- precip.sd[y]/precip.mn[y]
}
par(mfrow=c(1,3))
plot(year.vec, precip.sd, xlab="year", ylab="precipitation SD", pch=19, cex=0.7, type="b")
fit.precipsd <- lm(precip.sd ~ year.vec)
abline(fit.precipsd, lty=2, col="red")
abline(h=mean(precip.sd), lty=2, lwd=2)
linreg.ER(year.vec, precip.sd)

plot(year.vec, precip.cv, xlab="year", ylab="precipitation CV", pch=19, cex=0.7, type="b")
fit.precipcv <- lm(precip.cv ~ year.vec)
abline(fit.precipcv, lty=2, col="red")
abline(h=mean(precip.cv), lty=2, lwd=2)
linreg.ER(year.vec, precip.cv)

plot(precip.mn, precip.cv, xlab="mean precipitation", ylab="precipitation CV", pch=19, cex=0.7)
fit.precipmncv <- lm(precip.cv ~ precip.mn)
abline(fit.precipmncv, lty=2, col="red")
abline(h=mean(precip.cv), lty=2, lwd=2)
linreg.ER(precip.mn, precip.cv)
par(mfrow=c(1,1))

plot(year.vec, precip.cv, xlab="year", ylab="monthly precipitation CV", pch=19, cex=0.7, type="b")
fit.precipcv <- lm(precip.cv ~ year.vec)
abline(fit.precipcv, lty=2, col="red")
abline(h=mean(precip.cv), lty=2, lwd=2)
linreg.ER(year.vec, precip.cv)

# frequency of extremes
n.sds <- 1
yearly.sd <- sd(as.numeric(precip.yr.sum))
yearly.mn <- mean(as.numeric(precip.yr.sum))
upp.lim <- yearly.mn + (n.sds*yearly.sd)
low.lim <- yearly.mn - (n.sds*yearly.sd)
times.upp <- which(as.numeric(precip.yr.sum) > upp.lim)
times.low <- which(as.numeric(precip.yr.sum) < low.lim)
time.diffs.upp <- diff(times.upp)
time.diffs.low <- diff(times.low)
hist(time.diffs.upp)
hist(time.diffs.low)

par(mfrow=c(2,1))
plot(year.vec[times.upp[-1]], time.diffs.upp, pch=19, xlab="year", ylab="years between high extremes")
times.upp.fit <- lm(time.diffs.upp ~ year.vec[times.upp[-1]])
abline(times.upp.fit, lty=2, col="red")
linreg.ER(year.vec[times.upp[-1]], time.diffs.upp)
plot(year.vec[times.low[-1]], time.diffs.low, pch=19, xlab="year", ylab="years between low extremes")
times.low.fit <- lm(time.diffs.low ~ year.vec[times.low[-1]])
abline(times.low.fit, lty=2, col="red")
linreg.ER(year.vec[times.low[-1]], time.diffs.low)
par(mfrow=c(1,1))

