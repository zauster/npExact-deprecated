#######################################################################
##   Program Name   npMeanPaired                                      #
##   Purpose        Exact nonparametric test of Schlag (2008)         #
##                  for comparing expected values given matched pairs #
##                  for H0: E(Y1) >= E(Y2)  if Y1,Y2 in [low,up]      #
##                                                                    #
##   Authors           Christian Pechhacker and Karl Schlag           #
##   Date              04.03.2012                                     #
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
## npMeanPaired(y1, y2, low = 0, up = 5, alpha = 0.05, alternative = "less")
## npMeanPaired(runif(20), runif(20), low = 0, up = 1,
##              alternative = "greater", iterations = 2000)


npMeanPaired <- function(x1, x2, lower = 0, upper = 1, ## d = 0,
                         alpha = 0.05,
                         alternative = "two.sided",
                         epsilon = 1 * 10^(-6),
                         iterations = 5000,
                         max.iterations = 100000)
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

    if(alpha >= 1 | alpha <= 0)
        stop("Please supply a sensible value for alpha.")

    if(alternative != "greater" & alternative != "less" & alternative !=
       "two.sided")
        stop("Please specify the alternative you want to test. Possible value are: 'greater' (default), 'less' or 'two.sided'")

    n <- length(x1)

    sample.est <- c(mean(x1), mean(x2))

    ## Normalize vectors to [0,1]
    x1 <- (x1 - lower)/(upper - lower)
    x2 <- (x2 - lower)/(upper - lower)

    error <- 1

    ## rejMatrix will store the results of the McNemarTest Monte
    ## Carlo replications
    rejMatrix <- vector(mode = "numeric", length = 0)

    if(alternative == "two.sided")
    {
        ##
        ## first alternative at alpha / 2
        ##

        ## compute the theta that optimizes the type II error
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = 0.5, N = n,
                                 alpha = alpha / 2 - epsilon)
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = 0.5, N = n, alpha = alpha / 2 - epsilon)
        pseudoalpha <- (alpha/2) * theta$theta
        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     McNemarTestRandom(runif(n) < x1,
                                                       runif(n) < x2,
                                                       pseudoalpha)))
            rejUpper <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rejUpper - theta$theta)^2)
        }
        rejectionUpper <- ifelse(rejUpper > theta$theta, TRUE, FALSE)

        ##
        ## other alternative at alpha / 2
        ##
        error <- 1
        rejMatrix <- vector(mode = "numeric", length = 0)
        x1 <- 1 - x1
        x2 <- 1 - x2
        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     McNemarTestRandom(runif(n) < x1,
                                                       runif(n) < x2,
                                                       pseudoalpha)))
            rejLess <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rejLess - theta$theta)^2)
        }
        rejectionLess <- ifelse(rejLess > theta$theta, TRUE, FALSE)

        rej <- rejUpper + rejLess
        rejection <- ifelse(rejectionUpper + rejectionLess >= 1, TRUE, FALSE)
    }
    else
    {
        if(alternative == "greater")
        {
            x1 <- 1 - x1
            x2 <- 1 - x2
        }

        ## compute the theta that optimizes the type II error
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = 0.5, N = n,
                                 alpha = alpha - epsilon)
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = 0.5, N = n, alpha = alpha - epsilon)
        pseudoalpha <- alpha * theta$theta

        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     McNemarTestRandom(runif(n) < x1,
                                                       runif(n) < x2,
                                                       pseudoalpha)))
            rej <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rej - theta$theta)^2)
        }

        rejection <- ifelse(rej > theta$theta, TRUE, FALSE)
    }

    if(!is.null(iterations) & length(rejMatrix) < 1000)
        warning("Low number of iterations. Results may be inaccurate.")

    if(length(rejMatrix) >= max.iterations)
        warning(paste("The maximum number of iterations (",
                      format(max.iterations, scientific = FALSE),
                      ") was reached. Rejection may be very sensible to the choice of the parameters.", sep = ""))

    names(sample.est) <- c(paste("mean(", names.x1, ")", sep = ""),
                           paste("mean(", names.x2, ")", sep = ""))

    null.value <- 0
    names(null.value) <- "E[x2] - E[x1]" ##"mean difference"

    ## if rejection in a two.sided setting, we inform the user of the
    ## side of rejection
    if(rejection == TRUE & alternative == "two.sided")
    {
        alt.hypothesis <- paste("E(", names.x1, ")",
                                ifelse(rejectionUpper == TRUE, " < ", " > "),
                                "E(", names.x2, ")", sep = "")
    }

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
                   d.alternative = (optimaltypeII$root - 0.5) * 2 * (upper - lower),
                   typeIIerror = theta$typeII,
                   iterations = length(rejMatrix),
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
        h <- choose(n10 + n01, n01)/(2^(n10 + n01))
        if(prob <= (alpha + h))
        {
            res <- (alpha - prob + h)/h
        }
    }
    return(res)
}
