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

        theta <- calcOptimalTheta(n = n,
                                  ## c = sample.est/((upper.x1 - lower.x1)*(upper.x2 - lower.x2)),
                                  alpha = alpha)

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
        rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
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
                       typeIIerror = theta$typeII,
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
        ## k <- max(0, r - n2):b
        ## print(sum(dhyper(k, n1, n2, r)))
        phyper(ifelse(b >= r - n2, b, r - n2),
               n1, n2, r) ## would be also possible, however
        ## problems when r > n2!
    }

qfun <- function(n1, n2, s1, s2, alpha = 0.05)
    {
        ## term1 <- w(alpha - hfun(s1, s1 + s2, n1, n2))
        ## term2 <- alpha - hfun(s1 - 1, s1 + s2, n1, n2)
        ## term3 <- term2 / gfun(s1, s1 + s2, n1, n2)

        ## ## print(term1)
        ## ## print(term2)
        ## ## print(term3)

        ## res <- term1 + w(term2) * (1 - term1) * term3
        ## res

        term2 <- hfun(s1 - 1, s1 + s2, n1, n2)
        ifelse(alpha < term2, 0,
               ifelse(alpha >= hfun(s1, s1 + s2, n1, n2),
                      1,
                      (alpha - term2)/gfun(s1, s1 + s2, n1, n2)))
    }

phi1fun <- function(n1, n2, y1, y2, alpha = 0.05)
    {
        ## res <- 0
        ## for(s1 in 0:n1)
        ##     {
        ##         for(s2 in 0:n2)
        ##             {
        ##                 res <- res + dbinom(s1, n1, y1) * dbinom(s2, n2, y2) *
        ##                     qfun(n1, n2, s1, s2, alpha)
        ##             }
        ##     }
        ## print(res)

        ## res <- 0
        ## for(s1 in 0:n1)
        ##     {
        ##         res <- res + sum(ifelse(qfun(n1, n2, s1, 0:n2, alpha) > 0,
        ##                dbinom(s1, n1, y1)*dbinom(0:n2, n2, y2)*qfun(n1, n2, s1, 0:s2, alpha),
        ##                0))
        ##     }
        ## print(res)

        ## res <- NULL
        res <- sum(sapply(0:n1,
                      function(s1)
                      sum(ifelse(qfun(n1, n2, s1, 0:n2, alpha) > 0,
                                 dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * qfun(n1, n2, s1, 0:n2, alpha),
                       0))))
        res
    }

typeIIError <- function(n, p, q0, q1, theta, alpha = 0.05)
    {
        ## res <- 0
        ## for(j in 0:n)
        ##     {
        ##         res <- res + dbinom(j, n, p) * phi1fun(n - j, j, q0, q1,
        ##                                                alpha)
        ##     }
        ## print(res)

        res <- sum(sapply(0:n,
                          function(j)
                          dbinom(j, n, p) * phi1fun(n - j, j, q0, q1,
                                                    theta * alpha)))
        return((1 - res)/(1 - theta))
    }

typeIIError_pq0 <- function(pq0, n, theta, alpha, c = 0.125)
    {
        z <- min(pq0[2], 1 - c/(pq0[1]*(1 - pq0[1])))
        res <- typeIIError(n, pq0[1], z, z + c/(pq0[1]*(1 - pq0[1])),
                           theta, alpha)
        ## res <- ifelse(res >= 1, 1, res)
        res
    }

typeIIError_theta <- function(theta, n, pq0, alpha, c = 0.125)
    {
        typeIIError_pq0(pq0, n, theta, alpha, c)
    }

calcOptimalTheta <- function(theta = 0.2,
                             pq0 = c(0.5, 0.5),
                             n, alpha, c = 0.12)
    {
        ## maximization step
        control <- list()
        control$fnscale <- -1
        res <- optim(c(pq0[1], pq0[2]), typeIIError_pq0,
                     n = n, alpha = alpha,
                     theta = theta,
                     ## method = "L-BFGS-B",
                     ## lower = 0, upper = 1,
                     control = control)
        print(res)

        ## minimization step
        res <- optimize(typeIIError_theta, c(0, 1),
                        n = n, pq0 = res$par, alpha = alpha,
                        c = c)
        print(res)

        return(list(typeII = res$objective,
                    theta = res$minimum))
    }

calcOptimalTypeIIError <- function(theta, n, alpha, c = 0.12,
                                   typeIIgoal = 0.5)
    {
        (calcOptimalTheta(theta, n = n, alpha = alpha,
                          c = c)$typeII - typeIIgoal)^2
    }

## uniroot(calcOptimalTypeIIError, c(0, 1), n = 20, alpha = 0.05)
## optimize(calcOptimalTypeIIError, c(0, 1), n = 20, alpha = 0.05)
