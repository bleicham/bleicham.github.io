# analyzing tick data using MCMC

dat <- readRDS("ticks.RDS")

Z <- dat$Result0   # test 
M <- dat$PoolSize  # pool size

# lets assume that the prior is a beta 
# with parameter 
a <- 1
b <- 2

# set up variables
n.smp <- 10000 # number of samples
sigma <- 0.1   # step size

logitP <- rep(0, n.smp)
LL <- rep(0,n.smp)

# calculate loglikelihood
pp.new <- exp( logitP[1])/(1+ exp( logitP[1] ) )
qq.new <- exp( M*log(1-pp.new) ) # calculate probability of escaping infection
LL[1] <- sum(Z*log(1-qq.new) ) + sum((1-Z)*log(qq.new)) + (a-1)*log(pp.new) + (b-1)*log(1-pp.new)

# run the Metropolis sampler
for ( k in 2:n.smp ){
  # generate proposal
  theta.new <- logitP[k-1] + sigma*rnorm(1)
  pp.new <- exp( theta.new )/(1 + exp( theta.new ) )
  qq.new <- exp( M*log(1-pp.new) )
  LL.new <- sum(Z*log(1-qq.new) ) + sum((1-Z)*log(qq.new)) + (a-1)*log(pp.new) + (b-1)*log(1-pp.new)
  
  # calculate ratio
  R <- min( exp( LL.new - LL[k-1] ), 1 )
  U <- runif(1)
  
  # move chain
  if ( U < R ){
    # go to proposal
    LL[k] <- LL.new
    logitP[k] <- theta.new
  } else {
    # stay put
    LL[k] <- LL[k-1]
    logitP[k] <- logitP[k-1]
  }
}

# index after burn in
burn.in <- 500
idx <- burn.in:n.smp

# histogram
P <- exp(logitP)/( 1 + exp(logitP ) )
hist(P[idx], nclass=100, xlab="propensity",main="CCHV")

# summary statistics for inference
sum.stat <- list(
  mean = mean(P[idx]),
  qq = quantile(P[idx],c(0.1,0.5,0.9) )
)

sapply(sum.stat, round, 3)
  
# trace of MCMC
plot( logitP, xlab="index", cex=0.2, pch=20 )
