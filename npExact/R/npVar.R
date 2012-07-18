######################################################################################
###         Testing Variance of Singe Sample - Upper and Lower Confidence Bounds    ##
###                             R Code                                              ##
######################################################################################

##----------------------------------------------------------------------------------
## Reference: K. H. Schlag (2008), A New Method for Constructing Exact Tests
## without Making any Assumptions, Universitat Pompeu Fabra working paper Working 1109
##----------------------------------------------------------------------------------

## This is a program that implements a nonrandomized one-sided test
## for the variance of a single sample of independent observations.
## Null hypothesis: H0: VarY >= w , if type="upper", finetuning test for H1: VarY <= w1
## Null hypothesis: H0: VarY <= w , if type="lower", finetuning test for H1: VarY >= w1

## This test can be used to derive (1-alpha) UPPER CONFIDENCE BOUNDS, type="upper":
## Find largest w such that not rejects the null

## This test can be used to derive (1-alpha) LOWER CONFIDENCE BOUNDS, type="lower":
## Find smallest w such that not rejects the null.

## The test transforms the data into a new sample that has mean equal to
## variance of orginal process and then applies the A test.

## example
## z <- rnorm(101,5, 3)
## npVar(z, low=-4, up=12, w=4, w1=0.1, type="upper", alpha=0.05)
## npVar(z, low=-4, up=12, w=2, w1=4, type="lower", alpha=0.05)

npVar <- function(test1,low, up, w, w1, type="upper",
                        alpha=0.025, T=50000)
{
  if (min(test1-low,up-test1)<0)
    stop("values have to lie in [low,up]")

  if (w > (up-low)/4)
    stop("w too large")

  if (w1>=w & type=="upper")
    stop("need w1 < w when type is upper")

  if (w1<=w & type=="lower")
    stop("need w1 > w when type is lower")

  ## Computation of sample mean and variance for output
  m1 <- length(test1)
  m <- floor(m1/2)
  x <- (test1 - low)/(up - low)  ## Normalization so that x in [0,1]
  w0 <- w / (up - low)^2  ## normalized threshold
  p <- 2*w0  ## threshold for mean comparison below
  p1 <- 2*w1 / (up - low)^2  ## threshold for finetuning and computing typeII error

  mean1 <- sum(test1/m1)
  var1 <- sum((test1-mean1)^2/(m1-1))

  if (type=="upper")
    {
      it <- as.numeric(min_value(n=m, p=1-p, p1=1-p1, alpha=alpha))
      if (it[2]>=0.99)
        stop("decrease w1 so that typeII is below 1")
    }
  else
    {
      it <- as.numeric(min_value(n=m, p=p, p1=p1, alpha=alpha))
      if (it[2]>=0.99)
        stop("increase w1 so that typeII is below 1")
    }

  theta <- it[1]
  pseudoalpha <- alpha*theta ##size adjusted downward for deterministic test

  cat(paste("Sample N = ",m1 , ',   avg = ',
            round(mean1,3),',   var = ',
            round(var1,3),',  theta =  ',
            round(theta,3),',  alpha =  ',
            alpha, sep=""))
  cat("\n")

  ##Monte Carlo simulations
  rj <- 0 ##Probability of rejection under size pseudoalpha (PMP test)

  for(t in 1:T)
    {
      c1 <- sample(x)
      ##transformation into sample in [0,1] that has mean equal to 1/2 + Var(X)
      y <- (c1[c(1:m)*2]-c1[(c(1:m)*2-1)])^2

      ##Random transformation of [0,1] data into {0,p,1} data, later only use {0,1}
      s1 <- 0 ##number of 0 values in transformed data
      s2 <- 0 ##number of 1 values in transformed data
      q <- runif(m)
      s2 <- length(which(y - p > (q*(1-p))))
      s1 <- length(which(y< q*p))
      ##Evaluation of randomized binomial test, see if the number of s2 relative
      ##to (s1+s2) is significantly below p
      h1 <- 0
      k <- switch(type, upper=0:s2, lower=s2:(s1+s2))
      h1 <- sum((p^k)*((1-p)^(s1+s2-k))*choose((s1+s2),k))
      if(h1 <= pseudoalpha) ##(reject with probability 1)
        {rj <- rj + 1/T}
      else
        {
          h <- (p^s2)*((1-p)^s1)*choose(s1+s2,s2)
          if (h1 <= pseudoalpha + h) ##(reject with positive probability)
            {    rj = rj + ((pseudoalpha - h1 + h) / h) /T}
        }

    }
  if (rj >= theta)
    {
      cat(paste('H0: VarY ',
                switch(type, upper=">=", lower="<="),
                w,' REJECTED  (rj=  ',
                round(rj,4),'), typeII =',
                round(it[2],4),sep=""))
    }
  else
    {
      cat(paste('H0: VarY ',
                switch(type, upper=">=", lower="<="),
                w,' NOT rejected  (rj=  ',
                round(rj,4),'), typeII =',
                round(it[2],4), sep="") )
    }
  cat("\n")

}
