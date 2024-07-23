# analyzing tick data

dat <- readRDS("ticks2.RDS")

# look at the distribution of number of ticks in batch
table(dat$PoolSize) # access names element in list using $

Z <- dat$Result0   # test 
M <- dat$PoolSize  # pool size

# construct "sufficient statistics"
T1 <- sum(Z)
T2 <- sum((1-Z)*M)
T3 <- sum(Z*(M==2))

# sample from Beta
n.smp <- 10000
P <- rbeta( n.smp, T1+1, T2+1 )
W.log <-  T3*log(2-P) 

# normalize weights
W.log <- W.log - mean(W.log) - log(n.smp)
W <- exp(W.log)/sum(exp(W.log))

# plot weights
# plot(W,pch=20,cex=0.2) # pch=20 full circle, cex=size of dots
# abline( h=1/n.smp, lwd=3,col="red") # uniform weight

# reuse code from previous example to make histogram
# 
nbins <- 200 # number of bins
histogram.breaks <- seq(0,0.5,length=nbins+1)
iidx <- cut( P, histogram.breaks )
split( W, iidx ) %>%
  sapply(., sum) -> histogram.hight 
histogram.hight <- histogram.hight * nbins  # normalization 

plot(histogram.breaks[-(nbins+1)],histogram.hight,type="s",
     xlab="probability",ylab="histogram",
     sub="histogram of posterior distribution")

# quantiles of posterior distribution
# sort sampled values and associated weights
idx <- sort.list(P)
P.sort <- P[idx]
W.sort <- W[idx]

# calculate cumulative sum (cdf)
Wcdf <- cumsum(W.sort)

sum.stat <- list(q10=max(P.sort[Wcdf < 0.1]),
     q50=max(P.sort[Wcdf < 0.5]),
     q90=max(P.sort[Wcdf < 0.9]) )

sapply(sum.stat, round, 4)