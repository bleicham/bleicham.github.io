setwd("~/Lecture/Lecture=current/Introduction")

flu.data <- read.table("Maryland_incidence.csv", 
                       sep=",", header=TRUE )
flu.data[,1] <- as.Date(flu.data[,1],format = "%d/%m/%Y")
pop.maryland <- scan("pop_Maryland_1920.csv",
                     sep=",", skip=1, what=0)
names(pop.maryland) <- scan("pop_Maryland_1920.csv", n=5, what="", sep=",")

tot.infected <- apply(flu.data[,-1],2,sum)

n.iter <- 10000
logit.pi <- rep( 0, n.iter )
LL <- rep( 0,n.iter )
AA <- rep(0, n.iter )

#. initialize variables
sdd <- 0.1
logit.pi[1] <- 0

ssize <- 5

mu <- pop.maryland * exp( logit.pi[1])/(1+exp( logit.pi[1]) )
LL[1] <- sum( dnbinom(tot.infected, ssize, mu=mu, log=TRUE) )

for ( k in 2:n.iter ){
  
  #. propose a new value
  new.logit.pi <- logit.pi[k-1] + rnorm(1,mean=0,sd=sdd )
  
  #. evaluate the loglikelihood
  new.mu <- pop.maryland * exp( new.logit.pi )/(1+exp( new.logit.pi ) )
  new.LL <- sum( dnbinom(tot.infected, ssize, mu=new.mu, log=TRUE) )
  
  RR <- exp( new.LL - LL[k-1] )
  UU <- runif(1)
  
  if ( UU < RR ){
    AA[k] <- 1
    LL[k] <- new.LL
    logit.pi[k] <- new.logit.pi
  } else {
    LL[k] <- LL[k-1]
    logit.pi[k] <- logit.pi[k-1]
  }
}

#. transform into a distribution for varrho
ppi <- exp(logit.pi)/(1 + exp(logit.pi) )
varrho <- -log( 1-ppi )/ppi

plot(varrho[1:200],pch=20,cex=0.5,
     xlab="index",ylab="estimate")

#. make histogram of value after burn in
hist( varrho[1001:n.iter], nclass=500, xlim=c(1,1.2),
      xlab="basic reproductive number",
      main="histogram of estimated reproductive number")

acf(varrho)

MovingAve <- function(x,n.lag){
  n.x <- length(x)
  M <- suppressWarnings( matrix(x,n.x+1,n.lag,byrow = FALSE) )
  m <- apply(M,1,mean)
  return(m[1:(n.x-n.lag+1)])
}

mA <- MovingAve(AA,50)
plot(mA)


