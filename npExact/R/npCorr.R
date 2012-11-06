###########################################################################
##                      Testing for Correlation                          ##
##                               R Code                                  ##
##                         disclaimer as usual                           ##
##                          by Karl H. Schlag                            ##
##                 July 3, 2008 - Copyrights reserved                    ##
###########################################################################

## Reference: UPF working paper 1097

## This is a program for testing the null hypothesis that the covariance
## of two random variables x_1 and x_2 is less than 0 against the
## alternative hypothesis that it is strictly above 0.

## So we are testing H_0: Cov(x_1,x_2) <= 0 against H_1: Cov(x_1,x_2) > 0
## based on a sample of n matched pairs, drawn independently from
## (x_1,x_2) where it is known ex-ante that x_1 and x_2 realize outcomes
## in [a1,b2] and [a2,b2] respectively

## TEST: Fix a threshold parameter theta.
## Transform the data randomly into a binary data.
## Apply the Fisher Tocher one-sided test.
## Repeat and record the probability of rejection rj.
## Reject null if rj is greater than theta.

## For H0: Cov(x_1,x_2) >= 0 simply replace X2 by -X2

## Note: cannot be reversed unless Var(x_2) known (conditional approach)

## Initial Values
## a1 = 0; #specify the ex-ante known lower bound of x_1
## b1 = 1000; #specify the ex-ante known upper bound of x_1
## a2 = 0; #specify the ex-ante known lower bound of x_2
## b2 = 1; #specify the ex-ante known upper bound of x_2
## T = 10000; #number of iterations
## theta = 0.2; #scale parameter (to be chosen ex-ante)
## alpha = 0.1; #size of the nonrandomized test
## pseudoalpha = alpha * theta; #size of the randomized test

npCorr <- function(x1, x2, a1 = 0, b1 = 1, a2 = 0, b2 = 1,
                       iterations = 5000, theta = 0.2, alpha = 0.1,
                       pseudoalpha = theta * alpha,
                       alternative = "greater",
                       conditional = FALSE)
  {
    ## What this test does:
    ## Fix a threshold parameter theta.
    ## Transform the data randomly into a binary data.
    ## Apply the Fisher Tocher one-sided test.
    ## Repeat and record the probability of rejection rj.
    ## Reject null if rj is greater than theta.
    DNAME <- paste(deparse(substitute(x1)), "and",
                   deparse(substitute(x2)))

    x1 <- as.vector(x1)
    x2 <- as.vector(x2)

    ## Warnings
    if(min(x1) < a1 | max(x1) > b1 | (min(x2) < a2 | max(x2) > b2) &
    conditional == FALSE)
      stop("Some values are out of bounds!")

    if(length(x1) != length(x2))
      stop("Unequal length of input vectors!")

    if(iterations < 500)
      warning("Low number of iterations. Results may be inaccurate.")

    if(alternative == "two.sided")
      stop("Currently not supported. Please test for greater and less at alpha/2.")

    if(alpha >= 1 | alpha <= 0)
      stop("Please supply a sensible value for alpha.")

    n <- length(x1)
    rj <- 0

    ## Normalization of x1 and x2
    x1 <- (x1 - a1)/(b1 - a1)

    if(conditional == TRUE)
      {
        a2 <- min(x2)
        b2 <- max(x2)
      }
    x2 <- (x2 - a2)/(b2 - a2)

    if(alternative == "less")
      {
        x2 <- 1 - x2
      }

    for(t in 1:iterations)
      {
        bin.x1 <- as.numeric(runif(n) < x1) ## create a random
        ## transformation of x1 in [0, 1] to {0, 1}
        bin.x2 <- as.numeric(runif(n) < x2)

        a <- sum(bin.x1 == 1 & bin.x2 == 0) ## #{i: (x1, x2) = (1,0)}
        b <- sum(bin.x1 == 1 & bin.x2 == 1) ## #{i: (x1, x2) = (1,1)}
        n1 <- sum(bin.x2 == 0) ## #{i: (x1, x2) = (*, 0)}
        n2 <- sum(bin.x2 == 1) ## #{i: (x1, x2) = (*, 1)}

        if(n1 * n2 > 0)
          {
            r <- a + b
            h1 <- 0
            h2 <- 0
            if(a >= 1 + max(c(0, r - n2)))
              {
                for(k in max(c(0, r - n2)):(a - 1))
                  {
          ##           h1 <- h1 + (factorial(n1)/(factorial(k) *
          ## factorial(n1 - k))) * (factorial(n2)/(factorial(r - k) *
          ## factorial(n2 - (r -k)))) / (factorial(n1 + n2)/(factorial(r)
          ## * factorial(n1 + n2 - r)))
                    h1 <- h1 + choose(n1, k) * choose(n2, r - k) /
          choose(n1 + n2, r)
                  }
              }

            if(h1 <= pseudoalpha)
              {
          ##       h <- (factorial(n1)/(factorial(a) *
          ## factorial(n1-a))) * (factorial(n2)/(factorial(r-a) *
          ## factorial(n2-(r-a)))) / (factorial(n1+n2)/(factorial(r) *
          ## factorial(n1+n2-r)))
                h <- choose(n1, a) * choose(n2, r - a) / choose(n1 +
              n2, r)

                h <- h + h1
                if(h <= pseudoalpha)
                  {
                    rj <- rj + 1/iterations
                  }
                else
                  {
                    rj <- rj + (pseudoalpha - h1)/(iterations * (h - h1))
                  }
              }
          }
        else
          {
            rj <- rj + pseudoalpha/iterations
          }
      }

    method <- paste(ifelse(conditional == TRUE,
                           "Conditional", "Unconditional"),
                    "Nonparametric Correlation Test")
    sample.est <- cov(x1, x2)
    names(sample.est) <- "covariance"
    null.value <- 0
    names(null.value) <- "covariance"
    rejection <- ifelse(rj >= theta, TRUE, FALSE)
    bounds <- paste("[", a1, ", ", b1,
                    "] and [", a2, ", ", b2, "]",
                    sep = "")

    structure(list(method = method,
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

w <- function(x)
  {
    as.numeric(x > 0)
  }

gfun <- function(a, r, m, n)
  {
    dhyper(a, m, n, r)
  }

hfun <- function(b, r, n1, n2)
  {
    k <- max(0, r - n2):b
    ## res <- (choose(n1, k) * choose(n2, r - k))/choose(n1 + n2, r)
    ## print(sum(res))

    ## res <- phyper(b, n1, n2, r)
    res <- sum(dhyper(k, n1, n2, r))
    res
  }

qfun <- function(n1, n2, s1, s2, alpha)
  {
    term1 <- (alpha - hfun(s1, s1 + s2, n1, n2)) > 0
    term2 <- (alpha - hfun(s1 - 1, s1 + s2, n1, n2)) > 0
    term3 <- (alpha - hfun(s1 - 1, s1 + s2, n1, n2))/gfun(s1, s1 + s2,
                                                          n1, n2)
    res <- term1 + term2*(1 - term1)*term3
    res
  }

phi1fun <- function(n1, n2, y1, y2, alpha)
  {
    res <- 0
    ## res1 <- 0
    for(s1 in 0:n1)
      {
        for(s2 in 0:n2)
          {
  ##           res1 <- res1 + choose(n1, s1) * y1^s1 * (1 - y1)^(n1 - s1) *
  ## choose(n2, s2) * y2^s2 * (1 - y2)^(n2 - s2) * qfun(n1, n2, s1, s2,
  ## alpha)
            res <- res + dbinom(s1, n1, y1) * dbinom(s2, n2, y2) *
  qfun(n1, n2, s1, s2, alpha)
          }
      }
    res
  }

phi2fun <- function(n, p, q0, q1, alpha)
  {
    res <- 0
    for(j in 0:n)
      {
  ##       res <- res + choose(n, j) * p^j * (1 - p)^(n - j) * phi1fun(n
  ## - j, j, q0, q1, alpha)
        res <- res + dbinom(j, n, p) * phi1fun(n - j, j, q0, q1,
  alpha)

      }
    res
  }

typeIIerror <- function(x, n, theta, alpha, c) ## x = (p, q0)
  {
    p <- x[1]
    q0 <- x[2]
    z <- min(q0, 1 - c/(p*(1 - p)))
    print(paste("p = ", p, ", q0 = ", q0))
    print(paste("z = ", z))
    print(paste("theta = ", theta))
    res <- (1 - phi2fun(n, p, z,
                        z + c/(p*(1 - p)),
                        theta*alpha))/(1 - theta)
    res
  }





## npCorrOptimalTheta <- function(n, alpha, c = 0.12)
##   {
##     pbounds <- c((1 - sqrt(1 - 4*c))/2,
##                  1 - (1 - sqrt(1 - 4*c))/2)

##     control <- list()
##     control$fnscale <- -1




##     f2 <- function(theta, alpha, c = c)
##       {
##         res <- optim(c((pbounds[1] + pbounds[2])/2, 0.5),
##                      f1, c = c,
##                      theta = theta, alpha = alpha,
##                      control = control,
##                      lower = c(pbounds[1], 0),
##                      upper = c(pbounds[2], 1),
##                      method = "L-BFGS-B")
##         ## print(res)
##         res$value
##       }
##     optim(0.4, f2, alpha = alpha, c = c,
##           lower = 0.01, upper = 0.99,
##           method = "L-BFGS-B")
##           ## method = "Brent")

##   }

## TODO
## check first function: optimization of p and q0
## -> seems not to work very good
## res <- matrix(0, nrow = 20, ncol = 20)
## index <- seq(0, 1, length.out = 20)
## for(i in 1:20)
##   {
##     for(j in 1:20)
##       {
##         res[i, j] <- f1(x = c(index[i], index[j]), 0.4, 0.05, c = 0.12)
##       }
##   }
## other optimization package?
## multidimensional optimization seems to be the problem.

## c <- 0.12
##     pbounds <- c((1 - sqrt(1 - 4*c))/2,
##                  1 - (1 - sqrt(1 - 4*c))/2)
## n <- 20
## j <- 1
## res <- NULL
## for(theta in seq(0, 1, length.out = 10)[-10])
##   {
##     res[j] <- f2(theta, 0.1, 0.12)
##     j <- j + 1
##   }
