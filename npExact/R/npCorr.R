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
                       iter = 1000, theta = 0.2, alpha = 0.1,
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

    ## Warnings
    if(min(x1) < a1 | max(x1) > b1 | (min(x2) < a2 | max(x2) > b2) &
    conditional == FALSE)
      stop("Some values are out of bounds!")

    if(length(x1) != length(x2))
      stop("Unequal length of input vectors!")

    if(iter < 500)
      warning("Low number of iterations. Results may be inaccurate.")

    if(alternative == "two.sided")
      stop("Currently not supported. Please test for greater and less at alpha/2.")

    DNAME <- paste(deparse(substitute(x1)), "and",
                   deparse(substitute(x2)))

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

    for(t in 1:iter)
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
                    rj <- rj + 1/iter
                  }
                else
                  {
                    rj <- rj + (pseudoalpha - h1)/(iter * (h - h1))
                  }
              }
          }
        else
          {
            rj <- rj + pseudoalpha/iter
          }
      }

    method <- paste(ifelse(conditional == TRUE,
                           "Conditional", "Unconditional"),
                    "Nonparametric Correlation Test")
    sample.est <- cov(x1, x2)
    names(sample.est) <- "covariance"
    null.value <- 0
    names(null.value) <- "covariance"
    rejection <- ifelse(rj > theta, TRUE, FALSE)
    bounds <- paste("[", a1, ", ", b1, "] and [", a2, ", ", b2, "]",
                    sep = "")

    structure(list(method = method,
                   data.name = DNAME,
                   alternative = alternative,
                   estimate = sample.est,
                   probrej = rj,
                   rejection = rejection,
                   alpha = alpha,
                   theta = theta,
                   iterations = iter,
                   pseudoalpha = pseudoalpha,
                   bounds = bounds,
                   ## bounds = list(c(a1, b1),
                   ##   c(round(a2, digits = 3), round(b2, digits = 3))),
                   null.value = null.value),
              class = "nphtest")
  }
