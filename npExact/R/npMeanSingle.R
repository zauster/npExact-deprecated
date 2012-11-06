
## npMeanSingle
## a function to test the mean of a variable

## x ... the data vector
## mu ... the value to be tested, the supposed true value of the mean
## lower, upper ... logical bounds of x, thus x in [lower, upper]
## iterations ... number of iterations
## alpha ... the level of the test
## alternative ... a character string specifying the alternative
## hypothesis, must be one of 'two.sided', 'greater' or 'less'.

## example
## npMeanSingle(runif(20), mu = 0.3)

npMeanSingle <- function(x, mu,
                         lower = 0, upper = 1,
                         iterations = 500, alpha = 0.05,
                         alternative = "greater")
{
  DNAME <- deparse(substitute(x))
  sample.est <- mean(x)
  x <- as.vector(x)

  ## warnings
  if (min(x) < lower | max(x) > upper)
    stop("Some values are out of bounds!")

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  if(mu <= lower | mu >= upper )
    stop("Please supply a sensible value for mu")

  if(alternative != "greater" & alternative != "less" & alternative !=
  "two.sided")
    stop("Please specify the alternative you want to test. Possible value are: 'greater' (default), 'less' or 'two.sided'")

  ## standardize variables
  x <- (x - lower)/(upper - lower)
  p <- (mu - lower)/(upper - lower)

  n <- length(x)
  xp <- x - p

  if(alternative == "two.sided")
    {
      theta.par <- optim(c(0.4, p/4),
                         optimizeTheta, alpha = alpha,
                         mu0 = p, N = n)$par
      theta <- theta.par[1]
      ## mu1 <- mu + theta.par[2] ## not needed right now
      pseudoalpha <- (alpha/2) * theta

      rj.greater <- mean(replicate(iterations,
                                   transBinomTest(x, p, xp, n,
                                                  pseudoalpha)))

      x <- 1 - x
      p <- 1 - p
      xp <- x - p
      theta.par <- optim(c(0.4, p/4),
                         optimizeTheta, alpha = alpha,
                         mu0 = p, N = n)$par
      theta <- theta.par[1]
      ## mu1 <- mu + theta.par[2] ## not needed right now

      rj.less <- mean(replicate(iterations,
                                transBinomTest(x, p, xp, n,
                                               pseudoalpha)))

      rj <- rj.greater + rj.less
    }
  else
    {
      if(alternative == "less")
        {
          x <- 1 - x
          p <- 1 - p
        }

      theta.par <- optim(c(0.4, p/4),
                         optimizeTheta, alpha = alpha,
                         mu0 = p, N = n)$par
      theta <- theta.par[1]
      ## mu1 <- mu + theta.par[2] ## not needed right now
      pseudoalpha <- alpha * theta

      rj <- mean(replicate(iterations,
                           transBinomTest(x, p, xp, n,
                                          pseudoalpha)))
    }

  method <- "Nonparametric Single Mean Test"
  names(sample.est) <- "mean"
  null.value <- mu
  names(null.value) <- "mean"
  rejection <- ifelse(rj >= theta, TRUE, FALSE)
  bounds <- paste("[", lower, ", ", upper, "]", sep = "")

  structure(list(method = method,
                   data.name = DNAME,
                   alternative = alternative,
                   estimate = sample.est,
                   probrej = rj,
                   rejection = rejection,
                   alpha = alpha,
                   theta = theta,
                   mu.alternative = mu + theta.par[2],
                   iterations = iterations,
                   pseudoalpha = pseudoalpha,
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
}


## transBinomTest
## executes a binomial test on the (randomly) transformed data

## x ... data vector
## p ... transformed mu
## xp ... simply x - p
## n ... length of x
## pseudoalpha ... theta times alpha, the (new) level of the test

transBinomTest <- function(x, p, xp, n, pseudoalpha)
  {
    q <- runif(n)
    zeros <- sum(x < (q * p))  ## counts how often values of x < q*p
    ones <- sum(xp > (q * (1 - p))) ## counts how often xp > q*(1-p)
    res.binomtest <- 1 - pbinom(ones - 1, ones + zeros, p)
    res <- 0
    ## or
    if (res.binomtest <= pseudoalpha)
      {
        res <- 1
      }
    else
      {
        h2 <- (p^ones) * ((1 - p)^zeros) * choose(zeros + ones, ones)
        if (res.binomtest <= (pseudoalpha + h2))
          {
            res <- ((pseudoalpha - res.binomtest + h2)/h2)
          }
      }
    res
  }

## w
## helper function, to ease the reading of the code

w <- function(x)
  {
    as.numeric(x >= 0)
  }

## g1fun
## helper function, to calculate the type II error in
## npMeanSingleTypeIIerror

g1fun <- function(k, N, z, alpha)
  {
    summationterm1 <- alpha - (1 - pbinom(k - 1, N, z))
    summationterm2 <- alpha - (1 - pbinom(k, N, z))
    term3 <- summationterm2/dbinom(k, N, z)
    res <- w(summationterm1) + (1 - w(summationterm1)) *
  (w(summationterm2)) * term3
    res
  }

## g2fun
## helper function, to calculate the type II error in
## npMeanSingleTypeIIerror

g2fun <- function(alpha, mu, N, z)
  {
    k <- 0:(N - 1)
    term2 <- w(alpha - z^N) + (1 - w(alpha - z^N))*(alpha/z^N)
    res <- sum(dbinom(k, N, mu) * g1fun(k, N, z, alpha)) + mu^N*term2
    res
  }

## npMeanSingleTypeIIError
## calculates the type II error for the npMeanSingle function

## alpha ... level of the test
## theta ... parameter for the random testing
## mu.alternative ... value where the type II error is evaluated
## mu0 ... hypothesized "true" value of the mean (transformed)
## N ... length of x

npMeanSingleTypeIIError <- function(alpha, theta,
                                    mu.alternative, mu0, N)
  {
    res <- (1 - g2fun(alpha * theta,
                      mu.alternative, N, mu0))/(1 - theta)
    res
  }

## optimizeTheta
## finds the value of theta and d (difference of mu.alternative and
## mu0), such that the type II error is 0.5 at this value

## par ... the two parameters to be optimized
## alpha ... level of the test
## mu0 ... hypothesized "true" value of the mean
## N ... length of x
optimizeTheta <- function(par, alpha, mu0, N)
  {
    res <- (npMeanSingleTypeIIError(alpha, par[1],
                                    mu0 + par[2], mu0, N) - 0.5)^2
    res
  }
