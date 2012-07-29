###################################################
###         Testing Variance of Singe Sample     ##
###         Upper and Lower Confidence Bounds    ##
###                   R Code                     ##
###################################################

##-------------------------------------------------
## Reference: K. H. Schlag (2008), A New Method for Constructing Exact
## Tests without Making any Assumptions, Universitat Pompeu Fabra working
## paper Working 1109
##-------------------------------------------------

## This is a program that implements a nonrandomized one-sided test for
## the variance of a single sample of independent observations.

## Null hypothesis: H0: VarY >= w, if type="upper",
##                      finetuning test for H1: VarY <= w1
## Null hypothesis: H0: VarY <= w, if type="lower",
##                      finetuning test for H1: VarY >= w1

## This test can be used to derive (1 - alpha) UPPER CONFIDENCE BOUNDS,
## type="upper": Find largest w such that not rejects the null

## This test can be used to derive (1 - alpha) LOWER CONFIDENCE BOUNDS,
## type="lower": Find smallest w such that not rejects the null.

## The test transforms the data into a new sample that has mean equal to
## variance of orginal process and then applies the A test.

## example
## npVar(runif(50), low = 0, up = 1, w = 0.1, w1 = 0.01)

npVar <- function(x, low, up, w, w1,
                  type="upper",
                  alpha=0.025, iterations=500)
{
  DNAME <- deparse(substitute(x))

  if (min(x) < low | max(x) > up)
    stop("Some values are out of bounds!")

  if (w > (up - low)/4)
    stop("Parameter w is too large.")

  if (w1 >= w & type == "upper")
    stop("A w1 < w is needed when type is upper.")

  if (w1 <= w & type == "lower")
    stop("A w1 > w is needed when type is lower.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  ## Computation of sample mean and variance for output
  m <- floor(length(x) / 2)
  x <- (x - low)/(up - low)  ## Normalization so that x in [0,1]
  w0 <- w / (up - low)^2  ## normalized threshold
  p <- 2 * w0  ## threshold for mean comparison below
  p1 <- 2 * w1 / (up - low)^2  ## threshold for finetuning and computing typeII error

  ## if (type=="upper")
  ##   {
  ##     it <- as.numeric(min_value(n=m, p=1-p, p1=1-p1, alpha=alpha))
  ##     if (it[2]>=0.99)
  ##       stop("decrease w1 so that typeII is below 1")
  ##   }
  ## else
  ##   {
  ##     it <- as.numeric(min_value(n=m, p=p, p1=p1, alpha=alpha))
  ##     if (it[2]>=0.99)
  ##       stop("increase w1 so that typeII is below 1")
  ##   }

  ## theta <- it[1]
  theta <- 0.4 ## just as long the code above doesn't work
  pseudoalpha <- alpha * theta ## size adjusted downward for deterministic test

  ## Monte Carlo simulations
  rj <- 0 ## Probability of rejection under size pseudoalpha (PMP test)

  for(t in 1:iterations)
    {
      c1 <- sample(x)

      ## transformation into sample in [0,1] that has mean equal to 1/2 +
      ## Var(X)

      y <- (c1[c(1:m)*2] - c1[(c(1:m)*2 - 1)])^2

      ## Random transformation of [0,1] data into {0,p,1} data,
      ## later only use {0,1}
      q <- runif(m)

      ## number of 0 values in transformed data
      s2 <- sum(y - p > (q*(1 - p)))
      ## number of 1 values in transformed data
      s1 <- sum(y < q*p)

      ## Evaluation of randomized binomial test, see if the number of s2
      ## relative to (s1+s2) is significantly below p

      k <- switch(type,
                  upper = 0:s2,
                  lower = s2:(s1+s2))
      h1 <- sum((p^k)*((1-p)^(s1 + s2 - k))*choose((s1 + s2), k))
      if(h1 <= pseudoalpha) ## reject with probability 1
        {
          rj <- rj + 1/iterations
        }
      else
        {
          h <- (p^s2)*((1 - p)^s1) * choose(s1 + s2, s2)
          if (h1 <= pseudoalpha + h) ##(reject with positive probability)
            {
              rj = rj + ((pseudoalpha - h1 + h) / h) /iterations
            }
        }

    }

  sample.est <- var(x)
  names(sample.est) <- "variance"
  null.value <- w1
  names(null.value) <- "variance"
  bounds <- paste("[", low, ", ", up, "]", sep = "")
  rejection <- ifelse(rj >= theta, TRUE, FALSE)
  alternative <- ifelse(type == "upper", "less", "greater")

  structure(list(method = "Nonparametric Variance Test",
                 data.name = DNAME,
                 alternative = alternative,
                 estimate = sample.est,
                 probrej = rj,
                 rejection = rejection,
                 alpha = alpha,
                 theta = theta,
                 iterations = iterations,
                 pseudoalpha = pseudoalpha,
                 bounds = bounds,
                 null.value = null.value),
            class = "nphtest")
}
