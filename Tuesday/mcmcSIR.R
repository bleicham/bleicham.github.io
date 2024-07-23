# set random seed.  Helps reproduce results and ddebugging
set.seed(912)

# read dataset
flu.data <- read.table("Maryland_incidence.csv", 
                       sep=",", header=TRUE )
flu.data[,1] <- as.Date(flu.data[,1],format = "%d/%m/%Y")
pop.maryland <- scan("pop_Maryland_1920.csv",
                     sep=",", skip=1, what=0)
names(pop.maryland) <- scan("pop_Maryland_1920.csv", n=5, what="", sep=",")

# 
n.days <- 22        # consider only the first three weeks of data
SS <- rep(0,n.days) # susceptible
II <- rep(0,n.days) # infected

# S depends only on population size and incidence, and is observed
Z <-flu.data[1:n.days,2]
SS <- pop.maryland[1] - cumsum(Z)
II[1] <- Z[1]

# set-up variables
n.iter <- 20000

logit.q  <- rep(0, n.iter)
logit.pi <- rep(0,n.iter)
VV <- matrix(0,n.iter,n.days)
LL <- rep(0,n.iter)

# initialize logit pi
logit.pi[1] <- -5
logit.q[1] <- -2

qq <- exp(logit.q[1])/(1+exp(logit.q[1]))

# calculate loglikelihood to initialize the chain
for ( j in 2:n.days ){
  VV[1,j] <- rbinom(1,II[j-1],qq)
  II[j] <- II[j-1]-VV[1,j]+Z[j]
}

ppi <- exp(logit.pi[1])/(1+exp(logit.pi[1]) )
pp  <- 1 - exp( II[-n.days] * log( 1 - ppi ) )

# loglikelihood
ll <- sum( Z[-1]*log(pp) + ( SS[-n.days]-Z[-1] )*log(1-pp) )

# need to include loglikelihood of recovered individuals
ll.V <- sum( dbinom(VV[1,-1],II[-n.days],qq,log = TRUE) )
LL[1] <- ll + ll.V

# step size
sdp <- 0.05

# mcmc

for ( k in 2:n.iter ){
  #. Metropolis proposal
  logit.pi.new <- logit.pi[k-1] + rnorm(1,mean=0,sd=sdp)
  pi.new <- exp( logit.pi.new )/( 1 + exp( logit.pi.new ) )
  
  logit.q.new <- logit.q[k-1] + sdp*rnorm(1)
  qq.new <- exp(logit.q.new)/(1+exp(logit.q.new) )
  # sample recovery (data augmentation step)
  
  for ( j in 2:n.days ){
    VV[k,j] <- rbinom(1,II[j-1],qq.new)
    II[j] <- II[j-1]-VV[k,j]+Z[j]
  }
  

  # Calculate likelihood of full data
  pp.new <- 1-exp( II[-n.days]*log(1-pi.new) )
  LL.new <- sum( Z[-1]*log(pp.new) + (SS[-n.days]-Z[-1])*log(1-pp.new) )
  ll.V <- sum( dbinom(VV[1,-1],II[-n.days],qq.new,log = TRUE) )
  
  LL.new <- LL.new + ll.V
  
  # make move
  R <- exp( min(LL.new-LL[k-1],0) )
  U <- runif(1)
  if ( U < R ){
    # accepts
    logit.pi[k] <- logit.pi.new
    logit.q[k] <- logit.q.new
    LL[k] <- LL.new
  } else {
    # reject 
    logit.pi[k] <- logit.pi[k-1]
    logit.q[k] <- logit.q[k-1]
    LL[k] <- LL[k-1]
  }
}

burn.in <- 5000
iidx <- burn.in:n.iter
jd <- kde2d(logit.pi[iidx],logit.q[iidx],n=100)
contour(jd,xlab="logit p",ylab="logit q",
        lwd=2, nlevels = 15,
        sub="joint posterior distribtuion",
        col=hcl.colors(15, "Spectral"))
points(logit.pi[iidx],logit.q[iidx],pch=".")

plot(logit.pi,pch=20,cex=0.2,
     xlab="iteration",ylab="logit of infection",
     sub="Trace of parameter")
plot(logit.q,pch=20,cex=0.2,
     xlab="iteration",ylab="logit of infection",
     sub="Trace of parameter")

burn.in <- 1000
XX <- matrix(0,n.iter-burn.in, n.days-1)
V  <- rep(0,n.days)
II <- rep(0,n.days)
II[1] <- Z[1]

for ( k in 1:(n.iter-burn.in) ){
  pi <- exp(logit.pi[k+burn.in])/(1+exp(logit.pi[k+burn.in]))
  for ( j in 2:n.days ){
    pp <- 1-exp(II[j-1]*log(1-pi))
    V[j] <- rbinom(1,II[j-1],qq)
    XX[k,j-1] <- rbinom(1,SS[j-1],pp)
    II[j] <- II[j-1]-V[j]+Z[j]
    if ( II[j] == 0 ) break
  }
}

boxplot(as.data.frame(XX),xlab="",ylim=c(0,75))
points(1:21,Z[-1],pch=20,col=2)
