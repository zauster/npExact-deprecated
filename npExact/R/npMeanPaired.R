#######################################################################
##   Program Name   npMeanPaired                                      #
##   Purpose        Exact nonparametric test of Schlag (2008)         #
##                  for comparing expected values given matched pairs #
##                  for H0: E(Y1) >= E(Y2)  if Y1,Y2 in [low,up]      #
##                                                                    #
##   Input variables                                                  #
##        required:   y1       independent sample of Y1               #
##                    y2       independent sample of Y2               #
##                    low      exogenously known lower                #
##                             bound for any outcome that can         #
##                             be generated in either sample          #
##                    up       exogenously known upper                #
##                             bound for any outcome that can         #
##                             be generated in either sample          #
##                    d        difference (d= EY2 - EY1)              #
##                             s.t. type II is minimized              #
##                             when H1: E(Y1)+d <= E(Y2)              #
##                             so 0< d <= up-low                      #
##                                                                    #
##        optional (default): times  (30000) number of Monte Carlo    #
##                                   iterations                       #
##                            alpha  (0.05) significance level of test#
##                                                                    #
##   Return Values                                                    #
##                            obs    number of observations in sample #
##                            low    lower bound for data values      #
##                            up     upper bound for data values      #
##                            theta  cutoff level theta               #
##                            alpha  significance level alpha         #
##                            avg1   average of values in y1          #
##                            avg2   average of values in y2          #
##                            P(rj)  probability of rejection in the  #
##                                   under lying randomdomized test   #
##                                                                    #
##   Authors           Christian Pechhacker and Karl Schlag           #
##   Date              04.03.2012                                     #
##                                                                    #
##                                                                    #
## For documentation see                                              #
##   (Schlag, Karl H. 2008, A New Method for Constructing Exact       #
##   Tests without Making any Assumptions, Department of              #
##   Economics and Business Working Paper 1109, Universitat           #
##   Pompeu Fabra)                                                    #
##                                                                    #
#######################################################################


## Examples:

## y1 <- sample(c(1,2), size = 100, replace = TRUE)
## y2 <- sample(c(1,2,3,4), size = 100, replace = TRUE)
## npMeanPaired(y1, y2, low = 0, up = 5, alpha = 0.05)
## npMeanPaired(runif(20), runif(20), low = 0, up = 1,
##              alternative = "greater", iterations = 2000)


npMeanPaired <- function(x1, x2, lower = 0, upper = 1, ## d = 0,
                         alpha = 0.05,
                         alternative = "greater",
                         epsilon = 1 * 10^(-6),
                         iterations = 5000)
    {
        names.x1 <- deparse(substitute(x1))
        names.x2 <- deparse(substitute(x2))

        DNAME <- paste(names.x1, "and", names.x2)

        null.hypothesis <- paste("E(", names.x1, ")",
                                 ifelse(alternative == "less", " >= ",
                                        ifelse(alternative == "greater",
                                               " <= ",
                                               " = ")),
                                 "E(", names.x2, ")", sep = "")

        alt.hypothesis <- paste("E(", names.x1, ")",
                                ifelse(alternative == "less", " < ",
                                       ifelse(alternative == "greater",
                                              " > ",
                                              " != ")),
                                "E(", names.x2, ")", sep = "")

        x1 <- as.vector(x1)
        x2 <- as.vector(x2)

        if(any(is.na(c(x1, x2))) == TRUE)
            {
                warning("Pairs containing NA's were removed completely.")
                complete <- complete.cases(cbind(x1, x2))
                x1 <- x1[complete]
                x2 <- x2[complete]
            }

        ## Warnings
        if(min(x1, x2) < lower | max(x1, x2) > upper)
            stop("Some values are out of bounds!")

        if(length(x1) != length(x2))
            stop("Unequal length of input vectors!")

        if(iterations < 500)
            warning("Low number of iterations. Results may be inaccurate.")

        if(alpha >= 1 | alpha <= 0)
            stop("Please supply a sensible value for alpha.")

        if(alternative != "greater" & alternative != "less" & alternative !=
           "two.sided")
            stop("Please specify the alternative you want to test. Possible value are: 'greater' (default), 'less' or 'two.sided'")

        ## d <- d/(upper - lower)
        ## if(d > 1 | d < 0)
        ##   stop("Please supply a sensible value for d.")

        n <- length(x1)

        sample.est <- c(mean(x1), mean(x2))

#### Normalize vectors to [0,1]
        x1 <- (x1 - lower)/(upper - lower)
        x2 <- (x2 - lower)/(upper - lower)

        error <- 1
        i <- 1
        rejMatrix <- NULL

        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = 0.5, N = n,
                                 alpha = alpha)
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = 0.5, N = n, alpha = alpha)

        if(alternative == "two.sided")
            {
                pseudoalpha <- (alpha/2) * theta$theta
                while(error > epsilon & i <= 20)
                    {
                        rejMatrix <- cbind(rejMatrix,
                                           replicate(iterations,
                                                     McNemarTestRandom(runif(n) < x1,
                                                                       runif(n) < x2,
                                                                       pseudoalpha)))
                        rejUpper <- mean(rejMatrix)
                        error <- exp(-2 * (iterations * i) * (rejUpper - theta$theta)^2)
                        i <- i + 1
                    }
                error <- i <- 1
                rejMatrix <- NULL

                x1 <- 1 - x1
                x2 <- 1 - x2
                while(error > epsilon & i <= 20)
                    {
                        rejMatrix <- cbind(rejMatrix,
                                           replicate(iterations,
                                                     McNemarTestRandom(runif(n) < x1,
                                                                       runif(n) < x2,
                                                                       pseudoalpha)))
                        rejLess <- mean(rejMatrix)
                        error <- exp(-2 * (iterations * i) * (rejLess - theta$theta)^2)
                        i <- i + 1
                    }


                rej <- rejUpper + rejLess
            }
        else
            {
                if(alternative == "greater")
                    {
                        x1 <- 1 - x1
                        x2 <- 1 - x2
                    }
                pseudoalpha <- alpha * theta$theta
                while(error > epsilon & i <= 20)
                    {
                        rejMatrix <- cbind(rejMatrix,
                                           replicate(iterations,
                                                     McNemarTestRandom(runif(n) < x1,
                                                                       runif(n) < x2,
                                                                       pseudoalpha)))
                        rej <- mean(rejMatrix)
                        error <- exp(-2 * (iterations * i) * (rej - theta$theta)^2)
                        i <- i + 1
                    }
            }

        if(i == 21)
            warning("The maximum number of iterations (100,000) was reached. Rejection may be very sensible to the choice of the parameters.")

        names(sample.est) <- c(paste("mean(", names.x1, ")", sep = ""),
                               paste("mean(", names.x2, ")", sep = ""))

        null.value <- 0
        names(null.value) <- "E[x2] - E[x1]" ##"mean difference"
        rejection <- ifelse(rej > theta$theta, TRUE, FALSE)
        bounds <- paste("[", lower, ", ", upper, "]", sep = "")

        structure(list(method = "Nonparametric Mean Test for Matched Pairs",
                       data.name = DNAME,
                       alternative = alternative,
                       null.hypothesis = null.hypothesis,
                       alt.hypothesis = alt.hypothesis,
                       estimate = sample.est,
                       probrej = rej,
                       rejection = rejection,
                       alpha = alpha,
                       theta = theta$theta,
                       d.alternative = (optimaltypeII$root - 0.5)*2*(upper - lower),
                       typeIIerror = theta$typeII,
                       iterations = iterations * (i - 1),
                       pseudoalpha = pseudoalpha,
                       bounds = bounds,
                       null.value = null.value),
                  class = "nphtest")
    } ## end of npMeanPaired


McNemarTestRandom <- function(x1, x2, alpha)
    {
        ## performs the randomized McNemar test
        ## x1, x2 are binary-valued vectors of equal length
        ## n10, n01 ... counts of (1,0) and (0,1) respectively.
        ##              (1, 1) and (0, 0) are dropped

        ## returns either 1 (rejection), p (prob of rejection) or 0 (no
        ## rejection)

        n10 <- sum(x1 > x2)
        n01 <- sum(x1 < x2)

        k <- n01:(n10 + n01)
        prob <- sum(choose(n10 + n01, k)/(2^(n10 + n01)))

        res <- 0
        if(prob <= alpha)
            {
                res <- 1
            }
        else
            {
                h <- choose(n10 + n01, n01)
                if(prob <= alpha)
                    {
                        res <- (alpha - prob + h)/h
                    }
            }
        return(res)
    }
