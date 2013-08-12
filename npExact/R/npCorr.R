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
## in [lower.x1,upper.x2] and [lower.x2,upper.x2] respectively

## TEST: Fix a threshold parameter theta.
## Transform the data randomly into a binary data.
## Apply the Fisher Tocher one-sided test.
## Repeat and record the probability of rejection rej.
## Reject null if rej is greater than theta.

## For H0: Cov(x_1,x_2) >= 0 simply replace X2 by -X2

## Note: cannot be reversed unless Var(x_2) known (conditional approach)

## Initial Values
## lower.x1 = 0; #specify the ex-ante known lower bound of x_1
## upper.x1 = 1000; #specify the ex-ante known upper bound of x_1
## lower.x2 = 0; #specify the ex-ante known lower bound of x_2
## upper.x2 = 1; #specify the ex-ante known upper bound of x_2
## T = 10000; #number of iterations
## theta = 0.2; #scale parameter (to be chosen ex-ante)
## alpha = 0.1; #size of the nonrandomized test
## pseudoalpha = alpha * theta; #size of the randomized test

npCorr <- function(x1, x2, lower.x1 = 0, upper.x1 = 1,
                   lower.x2 = 0, upper.x2 = 1,
                   iterations = 5000, alpha = 0.05,
                   alternative = "greater",
                   epsilon = 1 * 10^(-6),
                   conditional = FALSE)
  {
    ## What this test does:
    ## Fix a threshold parameter theta.
    ## Transform the data randomly into a binary data.
    ## Apply the Fisher Tocher one-sided test.
    ## Repeat and record the probability of rejection rej.
    ## Reject null if rej is greater than theta.
    names.x1 <- deparse(substitute(x1))
    names.x2 <- deparse(substitute(x2))
    DNAME <- paste(names.x1, "and", names.x2)

    x1 <- as.vector(x1)
    x2 <- as.vector(x2)
    sample.est <- cov(x1, x2)

    ## H_0: Cov(x_1,x_2) <= 0 against H_1: Cov(x_1,x_2) > 0
    null.hypothesis <- paste("Cov(", names.x1, ", ", names.x2, ")",
                             ifelse(alternative == "greater", " <= ",
                                    ifelse(alternative == "less", " >= ",
                                           " = ")),
                             "0", sep = "")
    alt.hypothesis <- paste("Cov(", names.x1, ", ", names.x2, ")",
                             ifelse(alternative == "greater", " > ",
                                    ifelse(alternative == "less", " < ",
                                           " != ")),
                             "0", sep = "")

    ## if(is.null(upper.x1) & is.null(lower.x1))
    ##     {
    ##         upper.x1 <- upper.x2
    ##         lower.x1 <- lower.x2
    ##     }
    ## if(is.null(upper.x2) & is.null(lower.x2))
    ##     {
    ##         upper.x2 <- upper.x1
    ##         lower.x2 <- lower.x1
    ##     }

    ## Warnings
    if(min(x1) < lower.x1 | max(x1) > upper.x1 | (min(x2) < lower.x2 | max(x2) > upper.x2) &
       conditional == FALSE)
      stop("Some values are out of bounds!")

    if(length(x1) != length(x2))
      stop("Unequal length of input vectors!")

    if(iterations < 500)
      warning("Low number of iterations. Results may be inaccurate.")

    ## if(alternative == "two.sided")
    ##   stop("Currently not supported. Please test for greater and less at alpha/2.")

    if(alpha >= 1 | alpha <= 0)
      stop("Please supply a sensible value for alpha.")

    n <- length(x1)
    rej <- 0

    ## hardcoding, until it works
    theta <- NULL
    theta$theta <- 0.4

    ## Normalization of x1 and x2
    x1 <- (x1 - lower.x1)/(upper.x1 - lower.x1)

    if(conditional == TRUE)
      {
        lower.x2 <- min(x2)
        upper.x2 <- max(x2)
      }
    x2 <- (x2 - lower.x2)/(upper.x2 - lower.x2)

    error <- i <- 1
    rejMatrix <- NULL

    if(alternative == "two.sided")
      {
          pseudoalpha <- (alpha/2)* theta$theta
          while(error > epsilon & i <= 20)
            {
                rejMatrix <- cbind(rejMatrix,
                                   replicate(iterations,
                                             randomFischerTocherTest(x1, x2, n,
                                                                     pseudoalpha)))
                rejUpper <- mean(rejMatrix)
                error <- exp(-2 * (iterations * i) * (rejUpper - theta$theta)^2)
                i <- i + 1
            }

          x2 <- 1 - x2
          error <- i <- 1
          rejMatrix <- NULL
          while(error > epsilon & i <= 20)
            {
                rejMatrix <- cbind(rejMatrix,
                                   replicate(iterations,
                                             randomFischerTocherTest(x1, x2, n,
                                                                     pseudoalpha)))
                rejLess <- mean(rejMatrix)
                error <- exp(-2 * (iterations * i) * (rejLess - theta$theta)^2)
                i <- i + 1
            }
          rej <- rejUpper + rejLess
      }
        else
          {
              if(alternative == "less")
                {
                    x2 <- 1 - x2
                }
              pseudoalpha <- alpha * theta$theta
              while(error > epsilon & i <= 20)
                {
                    rejMatrix <- cbind(rejMatrix,
                                       replicate(iterations,
                                                 randomFischerTocherTest(x1, x2,
                                                                         n, pseudoalpha)))
                    rej <- mean(rejMatrix)
                    error <- exp(-2 * (iterations * i) * (rej - theta$theta)^2)
                    i <- i + 1
                }
          }
    if(i == 21)
      warning("The maximum number of iterations (100,000) was reached. Rejection may be very sensible to the choice of the parameters.")

    method <- paste(ifelse(conditional == TRUE,
                           "Conditional", "Unconditional"),
                    "Nonparametric Correlation Test")
    names(sample.est) <- "covariance"
    null.value <- 0
    names(null.value) <- "covariance"
    rejection <- ifelse(rej >= theta, TRUE, FALSE)
    bounds <- paste("[", lower.x1, ", ", upper.x1,
                    "] and [", lower.x2, ", ", upper.x2, "]",
                    sep = "")

    structure(list(method = method,
                   data.name = DNAME,
                   null.hypothesis = null.hypothesis,
                   alt.hypothesis = alt.hypothesis,
                   alternative = alternative,
                   estimate = sample.est,
                   probrej = rej,
                   rejection = rejection,
                   alpha = alpha,
                   theta = theta$theta,
                   iterations = iterations * (i - 1),
                   pseudoalpha = pseudoalpha,
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
  } ## end of npCorr

randomFischerTocherTest <- function(x1, x2, n, pseudoalpha)
  {
    bin.x1 <- as.numeric(runif(n) < x1) ## create a random
    ## transformation of x1 in [0, 1] to {0, 1}
    bin.x2 <- as.numeric(runif(n) < x2)

    a <- sum(bin.x1 == 1 & bin.x2 == 0) ## #{i: (x1, x2) = (1,0)}
    b <- sum(bin.x1 == 1 & bin.x2 == 1) ## #{i: (x1, x2) = (1,1)}
    zeros <- sum(bin.x2 == 0) ## #{i: (x1, x2) = (*, 0)}
    ones <- sum(bin.x2 == 1) ## #{i: (x1, x2) = (*, 1)}
    rej <- 0

    if(zeros * ones > 0)
      {
        r <- a + b
        ## prob <- phyper(b, ones, zeros, r, lower.tail = FALSE)
        k <- max(0, r - zeros):b
        prob <- 1 - sum(dhyper(k, ones, zeros, r))

        ## prob <- hfun(b, r, ones, zeros) ## safer variant than above
        if(prob <= pseudoalpha)
          {
            h2 <- prob + choose(zeros,
                              a) * choose(ones,
                                          r - a) / choose(zeros + ones,
                                                          r)
            rej <- ifelse(prob + h2 <= pseudoalpha, 1,
                          (pseudoalpha - prob)/(h2 - prob))
          }
        }
    else
      {
        rej <- pseudoalpha
      }

    return(rej)
  }

w <- function(x)
  {
    as.numeric(x >= 0)
  }

gfun <- function(a, r, m, n)
  {
    dhyper(a, m, n, r)
  }

hfun <- function(b, r, n1, n2)
  {
    k <- max(0, r - n2):b
    return(sum(dhyper(k, n1, n2, r)))
    ## phyper(b, n1, n2, r) ## would be also possible, however
    ## problems when r > n2!
  }

qfun <- function(n1, n2, s1, s2, alpha)
  {
    term1 <- w(alpha - hfun(s1, s1 + s2, n1, n2))
    term2 <- alpha - hfun(s1 - 1, s1 + s2, n1, n2)
    term3 <- term2 / gfun(s1, s1 + s2, n1, n2)

    ## print(term1)
    ## print(term2)
    ## print(term3)

    res <- term1 + w(term2) * (1 - term1) * term3
    res
  }

phi1fun <- function(n1, n2, y1, y2, alpha)
  {
    res <- 0
    for(s1 in 0:n1)
      {
        for(s2 in 0:n2)
          {
            res <- res + dbinom(s1, n1, y1) * dbinom(s2, n2, y2) *
              qfun(n1, n2, s1, s2, alpha)
          }
      }
    ## for(s1 in 0:n1)
    ##   {
    ##     res <- res + sum(dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) *
    ##           qfun(n1, n2, s1, 0:n2, alpha))
    ## qfun(.., 0:n2) is not the same as qfun(.., 0), .., qfun(.., n2)!
    ##   }
    res
  }

phi2fun <- function(n, p, q0, q1, alpha)
  {
    res <- 0
    for(j in 0:n)
      {
        res <- res + dbinom(j, n, p) * phi1fun(n - j, j, q0, q1,
                                               alpha)
      }
    res
  }

typeIIError <- function(n, p, q0, q1, theta, alpha)
  {
    res <- (1 - phi2fun(n, p, q0, q1,
                        theta * alpha))/(1 - theta)
    res
  }

typeIIError_pq0 <- function(pq0, n, theta, alpha, c = 0.125)
  {
    z <- min(pq0[2], 1 - c/(pq0[1]*(1 - pq0[1])))
    ## print(pq0)
    res <- typeIIError(n, pq0[1], z, z + c/(pq0[1]*(1 - pq0[1])),
                theta, alpha)
    ## res <- ifelse(res >= 1, 1, res)
    res
  }

optimal_pq0<- function(theta, n, alpha, pq0 = c(.5, .5), c = .125)
  {
    control <- list()
    control$fnscale <- -1
    res <- optim(c(pq0[1], pq0[2]), typeIIError_pq0,
                 n = n, theta = theta, alpha = alpha, c = c,
                 ## method = "L-BFGS-B",
                 control = control)
    res$value
  }

## optimize(optimal_pq0, c(0, 1), n = 20, alpha = .1, c = .12)

## optim(.5, optimal_pq0, n = n, alpha = alpha, method = "Brent",
##       lower = 0, upper = 1)

## typeIIError_maxpq0 <- function(theta, n, alpha, c = 0.125)
##   {
##     pbounds <- c((1 - sqrt(1 - 4*c))/2,
##                  1 - (1 - sqrt(1 - 4*c))/2)
##         z <- min(q0, 1 - 0.25/(p*(1 - p)))
##   }

## typeIIError_maxpq<- function(x, n, theta, alpha) ## x = (p, q0)
##   {
##     p <- x[1]
##     q0 <- x[2]
##     z <- min(q0, 1 - 0.25/(p*(1 - p)))
##     ## print(paste("p = ", p, ", q0 = ", q0))
##     ## print(paste("z = ", z))
##     ## print(paste("theta = ", theta))
##     res <- (1 - phi2fun(n, p, z,
##                         z + 0.25/(p*(1 - p)),
##                         theta*alpha))/(1 - theta)
##     res
##   }

## res <- optim(c(0.5, 0.5),
##              typeIIError_maxpq, n = n,
##              theta = theta, alpha = alpha)
##              ## control = control,
##              ## lower = c(pbounds[1], 0),
##              ## upper = c(pbounds[2], 1),
##              ## method = "L-BFGS-B")

## npCorrOptimalTheta <- function(n, alpha, c = 0.12)
##   {
## c <- 1/4
##     pbounds <- c((1 - sqrt(1 - 4*c))/2,
##                  1 - (1 - sqrt(1 - 4*c))/2)
##     control <- list()
##     control$fnscale <- -1
##     f2 <- function(theta, alpha, c = c)
##       {
        ## res <- optim(c((pbounds[1] + pbounds[2])/2, 0.5),
        ##              f1, c = c,
        ##              theta = theta, alpha = alpha,
        ##              control = control,
        ##              lower = c(pbounds[1], 0),
        ##              upper = c(pbounds[2], 1),
        ##              method = "L-BFGS-B")
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
