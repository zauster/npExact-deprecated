
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
                         iterations = 5000, alpha = 0.05,
                         alternative = "greater")
{
  DNAME <- deparse(substitute(x))
  sample.est <- mean(x)
  x <- as.vector(x)

  null.hypothesis <- paste("E(", DNAME, ") ",
                           ifelse(alternative == "greater", "<= ",
                                  ifelse(alternative == "less", ">= ",
                                         "= ")),
                           mu, sep = "")
  alt.hypothesis <- paste("E(", DNAME, ") ",
                          ifelse(alternative == "greater", "> ",
                                 ifelse(alternative == "less", "< ",
                                        "!= ")),
                          mu, sep = "")

  ## warnings
  if (min(x) < lower | max(x) > upper)
    stop("Some values are out of bounds!")

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  if(mu <= lower | mu >= upper )
    stop("Please supply a sensible value for mu.")

  if(lower >= upper)
    stop("Please supply sensible values for the bounds.")

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
      optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                               c(0, 1), p = p, N = n,
                               alpha = alpha)
      theta <- minTypeIIError(optimaltypeII[[1]],
                              p = p, N = n, alpha = alpha)

      ## print(theta) ## for debugging only
      ## theta.par <- optim(c(0.4, ifelse(p + p/4 > 1, (1 - p)/2, p/4)),
      ##                    optimizeTheta, alpha = alpha,
      ##                    mu0 = p, N = n)$par
      ## theta <- theta.par[1]
      ## mu1 <- mu + theta.par[2] ## not needed right now
      pseudoalpha <- (alpha/2) * theta$theta

      rej.greater <- mean(replicate(iterations,
                                   transBinomTest(x, p, xp, n,
                                                  pseudoalpha)))

      x <- 1 - x
      p <- 1 - p
      xp <- x - p
      optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                               c(0, 1), p = p, N = n,
                               alpha = alpha)
      theta <- minTypeIIError(optimaltypeII[[1]],
                              p = p, N = n, alpha = alpha)

      ## theta.par <- optim(c(0.4, ifelse(p + p/4 > 1, (1 - p)/2, p/4)),
      ##                    optimizeTheta, alpha = alpha,
      ##                    mu0 = p, N = n)$par
      ## theta <- theta.par[1]
      ## mu1 <- mu + theta.par[2] ## not needed right now
      ## print(theta)
      pseudoalpha <- (alpha/2) * theta$theta

      rej.less <- mean(replicate(iterations,
                                transBinomTest(x, p, xp, n,
                                               pseudoalpha)))

      rej <- rej.greater + rej.less
    }
  else ## alternative == "greater" => default
    {
      if(alternative == "less")
        {
          ## print("less")
          x <- 1 - x
          p <- 1 - p
          xp <- x - p
        }

      optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                               c(0, 1), p = p, N = n,
                               alpha = alpha)
      theta <- minTypeIIError(optimaltypeII[[1]],
                              p = p, N = n, alpha = alpha)
      ## theta.par <- optim(c(0.4, ifelse(p + p/4 > 1, (1 - p)/2, p/4)),
      ##                    optimizeTheta, alpha = alpha,
      ##                    mu0 = p, N = n)$par
      ## theta <- ifelse(theta.par[1] < 0.1, 0.1, theta.par[1])
      ## mu1 <- mu + theta.par[2] ## not needed right now
      pseudoalpha <- alpha * theta$theta

      rej <- mean(replicate(iterations,
                           transBinomTest(x, p, xp, n,
                                          pseudoalpha)))
    }
  ## print(p + theta.par[2])

  method <- "Nonparametric Single Mean Test"
  names(sample.est) <- "mean"
  null.value <- mu
  names(null.value) <- "mean"
  rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
  bounds <- paste("[", round(lower, digits = 3), ", ",
                  round(upper, digits = 3), "]", sep = "")

  structure(list(method = method,
                 data.name = DNAME,
                 alternative = alternative,
                 null.hypothesis = null.hypothesis,
                 alt.hypothesis = alt.hypothesis,
                 estimate = sample.est,
                 probrej = rej,
                 rejection = rejection,
                 alpha = alpha,
                 theta = theta$theta,
                 mu.alternative = optimaltypeII$root,
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

    ## Code below tests H0: p_true < p and returns p-value
    res.binomtest <- 1 - pbinom(ones - 1, ones + zeros, p)

    res <- 0
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

#### Functions to calculate the typeII error ####

## w
## helper function, to ease the reading of the code

## w <- function(x)
##   {
##     as.numeric(x >= 0)
##   }

## ## g1fun
## ## helper function, to calculate the type II error in
## ## npMeanSingleTypeIIerror

## g1fun <- function(k, N, z, alpha)
##   {
##     summationterm1 <- alpha - (1 - pbinom(k - 1, N, z))
##     summationterm2 <- alpha - (1 - pbinom(k, N, z))
##     term3 <- summationterm2/dbinom(k, N, z)
##     res <- w(summationterm1) + (1 - w(summationterm1)) *
##   (w(summationterm2)) * term3
##     res
##   }

## ## g2fun
## ## helper function, to calculate the type II error in
## ## npMeanSingleTypeIIerror

## g2fun <- function(alpha, mu, N, z)
##   {
##     k <- 0:(N - 1)
##     term2 <- w(alpha - z^N) + (1 - w(alpha - z^N))*(alpha/z^N)
##     res <- sum(dbinom(k, N, mu) * g1fun(k, N, z, alpha)) + mu^N*term2
##     res
##   }

## ## npMeanSingleTypeIIError
## ## calculates the type II error for the npMeanSingle function

## ## alpha ... level of the test
## ## theta ... parameter for the random testing
## ## mu.alternative ... value where the type II error is evaluated
## ## mu0 ... hypothesized "true" value of the mean (transformed)
## ## N ... length of x

## npMeanSingleTypeIIError <- function(alpha, theta,
##                                     mu.alternative, mu0, N)
##   {
##     res <- (1 - g2fun(alpha * theta,
##                       mu.alternative, N, mu0))/(1 - theta)
##     ifelse(res >= 1, 1, res)
##   }

## ## optimizeTheta
## ## finds the value of theta and d (difference of mu.alternative and
## ## mu0), such that the type II error is 0.5 at this value

## ## par ... the two parameters to be optimized, theta and the
## ## difference of the 'true' mean to the hypothesized mean
## ## alpha ... level of the test
## ## mu0 ... hypothesized "true" value of the mean
## ## N ... length of x
## optimizeTheta <- function(par, alpha, mu0, N)
##   {
##     ## for debugging purposes:
##     ## print(paste("theta: ", par[1], " mu.alt: ", mu0 + par[2]))
##     ## par[1] <- ifelse(par[1] < 0, par[1] <- 0.8,
##     ##                  ifelse(par[1] > 1, par[1] <- 0.9, par[1] <- par[1]))
##     ## par[2] <- ifelse(par[2] < 0, par[2] <- 0.1,
##     ##                  ifelse(par[2] > 1, par[2] <- 0.9, par[2] <- par[2]))
##     ## print(paste("theta: ", par[1], " mu.alt: ", mu0 + par[2]))


##     ## real function content:
##     res <- abs(npMeanSingleTypeIIError(alpha, par[1],
##                                     mu0 + par[2], mu0, N) - 0.5) ##^2
##     res
##   }

## ## Problems:
## ## If mu0 near 1 -> optimization procedure becomes unstable!
## optim(c(0.4, .05), optimizeTheta, alpha = .05, mu0 = .4, N = 50)

## ## If mu0 near 1, mu.alt -> 1 and theta -> 0 (even -> negative)!

## ## alpha <- 1:5/25
## ## theta <- 5:20/25
## ## mu0 <- 1:19/40
## ## mu.alt <- mu0 + .05

## ## for(i in alpha)
## ##   {
## ##     for(j in theta)
## ##       {
## ##         for(k in 1:19)
## ##           {
## ##             print(paste("TII:", round(npMeanSingleTypeIIError(i, j,
## ##                                                         mu.alt[k], mu0[k],
## ##                                                         50),
## ##                                       digits = 3),
## ##                         " alpha: ", i, " theta: ", j,
## ##                         " mualt: ", mu.alt[k], " mu0: ", mu0[k]))
## ##           }
## ##       }
## ##   }
