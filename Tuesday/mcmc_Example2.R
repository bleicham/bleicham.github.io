
# set random seed.  Helps reproduce results and ddebugging
set.seed(249)


# generate synthetic data

# set the parameters for simulating the data
p <- 0.15
q <- 0.4
n <- 40

# generate the sample, using the auxiliary variables A
N <- sample(seq(2,6,by=1),n,replace=TRUE)
A <- rbinom(n,1,q)
Z <- rep(0,n)
for ( k in 1:n ){
  if ( A[k] == 1 ) Z[k] <- rbinom(1,N[k]-1,p)
}

#
X <- Z==0
#
n.smp <- 30000 # number of steps
sdd <- 0.075 # step size

# set up variables
logitP <- rep(0,n.smp)
logitQ <- rep(0,n.smp)
LL <- rep(0,n.smp)

pnew <- exp(logitP[1])/(1+exp(logitP[1]))
qnew <- exp(logitQ[1])/(1+exp(logitQ[1]))

# calculate likelihood 
pn <- exp((N-1)*log(1-pnew))
LL[1] <- sum(X*log((1-qnew)+qnew*pn)) + sum((1-X)*(log(qnew)+Z*log(pnew)+(N-1-Z)*log(1-pnew)))

# run sampler
for ( k in 2:n.smp ){
  
  # make proposal
  logit.pnew <- logitP[k-1] + sdd*rnorm(1)
  logit.qnew <- logitQ[k-1] + sdd*rnorm(1)
  
  # transform into probability
  pnew <- exp(logit.pnew)/(1+exp(logit.pnew))
  qnew <- exp(logit.qnew)/(1+exp(logit.qnew))
  
  # calculate likelihood 
  pn <- exp((N-1)*log(1-pnew))
  LL.new <- sum(X*log((1-qnew)+qnew*pn)) + sum((1-X)*(log(qnew)+Z*log(pnew)+(N-1-Z)*log(1-pnew)))

  # calculate ratio
  R <- exp( min(LL.new-LL[k-1], 0) )
  U <- runif(1)
  #make update
  if ( U < R ){
    # proposal is the new value
    logitP[k] <- logit.pnew
    logitQ[k] <- logit.qnew
    LL[k] <- LL.new
  } else {
    # stay 
    logitP[k] <- logitP[k-1]
    logitQ[k] <- logitQ[k-1]
    LL[k] <- LL[k-1]
  }
}
  

P <- exp(logitP)/(1+exp(logitP))
Q <- exp(logitQ)/(1+exp(logitQ))

# contour plot
jd <- kde2d(P,Q,n=100)
contour(jd,xlab="p",ylab="q",
        lwd=2, nlevels = 15,
        sub="joint posterior distribtuion",
        col=hcl.colors(15, "Spectral"))

points(P,Q,pch=".")

# trace of estimates
plot(logitP,pch=20,cex=0.2)
plot(logitQ,pch=20,cex=0.2)


