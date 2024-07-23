# data
N <- 50
Z <- 8
m <- 3

# parameter of prior
a <- 2
b <- 3

# number of samples
n.smp <- 10000
P <- rep( 0, n.smp ) # value of the parameter
LL <- rep(0, n.smp ) # loglikelihood

P[1] <- 0.2
pp <- P[1]
LL[1] <- (Z+a-1)*log(pp) + (m*(N-Z) + b - 1)*log(1-pp) + Z*log(m-m*pp+pp*pp)

for ( k in 2:n.smp ){
  pp.new <- rbeta(1,a,b)
  # loglikelihood of proposal
  LL.new <- (Z+a-1)*log(pp.new) + (m*(N-Z) + b - 1)*log(1-pp.new) + Z*log(m-m*pp.new+pp.new*pp.new)
  logprior <- (a-1)*( log(P[k-1]) - log(pp.new) ) + (b-1)*( log(1-P[k-1]) - log(1-pp.new) )
  R <- exp( min(LL.new - LL[k-1] + logprior, 0) ) # ratio
  U <- runif(1) # random uniform
  # move
  if ( U < R ){
    # go to proposal
    P[k] <- pp.new
    LL[k] <- LL.new
  } else {
    # stay where you are
    P[k] <- P[k-1]
    LL[k] <- LL[k-1]
  }
}

# plot Markov chain
plot(P, xlab="iteration", ylab="P",
     pch=20, cex=0.2)

# histogram of sample
hist(P,nclass=100, xlab="infection probability", ylab="frequency",xlim=c(0,0.3),
     main="Histogram of p")

# auto correlation function
acf(P, main="ACF of MCMC for p")