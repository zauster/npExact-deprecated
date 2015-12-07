##' A test for the variance of a bounded random variable based on a single
##' sample of iid observations.
##' 
##' This test requires that the user knows upper and lower bounds before
##' gathering the data such that the properties of the data generating process
##' imply that all observations will be within these bounds. The data input
##' consists of a sequence of observations, each being an independent
##' realization of the random variable. No further distributional assumptions
##' are made.
##' 
##' This is a test of the null hypothesis \eqn{H_0: Var(X) \le v} against
##' \eqn{H_1 : Var(X) > v}.
##' 
##' This test randomly matches the data into pairs, then computes for each pair
##' the square of the difference and continues with the resulting sequence with
##' half as many observations as npMeanSingle. See the cited paper for more
##' information.
##' 
##' @param x a (non-empty) numeric vector of data values.
##' @param v the value of the variance to be tested as \eqn{H_0: Var(x) \le v}.
##' @param lower,upper the theoretical lower and upper bounds on the data
##' outcomes known ex-ante before gathering the data.
##' @param alternative a character string describing the alternative
##' hypothesis, can take values "greater", "less" or "two.sided"
##' @param alpha the type I error.
##' @param iterations the number of iterations used, should not be changed if
##' the exact solution should be derived.
##' @param epsilon the tolerance in terms of probability of the Monte Carlo
##' simulations.
##' @param ignoreNA if \code{TRUE}, NA values will be omitted. Default:
##' \code{FALSE}
##' @param max.iterations the maximum number of iterations that should be
##' carried out. This number could be increased to achieve greater accuracy in
##' cases where the difference between the threshold probability and theta is
##' small. Default: \code{10000}
##' @return A list with class "nphtest" containing the following components:
##' 
##' \item{method}{ a character string indicating the name and type of the test
##' that was performed.  } \item{data.name}{ a character string giving the
##' name(s) of the data.  } \item{alternative}{ a character string describing
##' the alternative hypothesis.  } \item{estimate}{ the estimated mean or
##' difference in means depending on whether it was a one-sample test or a
##' two-sample test.  } \item{probrej}{ numerical estimate of the rejection
##' probability of the randomized test, derived by taking an average of
##' \code{iterations} realizations of the rejection probability.  }
##' \item{bounds}{ the lower and upper bounds of the variables.  }
##' \item{null.value}{ the specified hypothesized value of the correlation
##' between the variables.  } \item{alpha}{ the type I error.  } \item{theta}{
##' the parameter that minimizes the type II error.  } \item{pseudoalpha}{
##' \code{theta}*\code{alpha}, this is the level used when calculating the
##' average rejection probability during the iterations.  } \item{rejection}{
##' logical indicator for whether or not the null hypothesis can be rejected.
##' } \item{iterations}{ the number of iterations that were performed.  }
##' @author Karl Schlag and Oliver Reiter
##' @seealso
##' \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
##' @references Karl Schlag (2008).  Exact tests for correlation and for the
##' slope in simple linear regressions without making assumptions. Available at
##' \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1097}.
##' @keywords variance test single sample
##' @examples
##' 
##' ## see if the minority share holder shores have a variance greater
##' ## than 0.05
##' data(mshscores)
##' 
##' scores <- as.vector(as.matrix(mshscores))
##' npVarianceSingle(scores, lower = 0, upper = 1, v = 0.05, ignoreNA = TRUE)
##' 
##' @export npVarianceSingle
npVarianceSingle <- function(x, v, lower = 0, upper = 1,
                             alternative = "two.sided",
                             alpha = 0.05, iterations = 5000,
                             epsilon = 1 * 10^(-6),
                             ignoreNA = FALSE,
                             max.iterations = 100000)
{
    DNAME <- deparse(substitute(x))
    x <- as.vector(x)

    null.hypothesis <- paste("Var(", DNAME, ") ",
                             ifelse(alternative == "greater", "<= ",
                             ifelse(alternative == "less", ">= ",
                                    "= ")),
                             v, sep = "")
    alt.hypothesis <- paste("Var(", DNAME, ") ",
                            ifelse(alternative == "greater", "> ",
                            ifelse(alternative == "less", "< ",
                                   "!= ")),
                            v, sep = "")

    if(ignoreNA == TRUE)
    {
        x <- x[!is.na(x)]
    }
    else if(any(is.na(x)) == TRUE)
    {
        stop("The data contains NA's!")
    }

    if(min(x) < lower | max(x) > upper)
        stop("Some values are out of bounds (or NA)!")

    if(v > 0.25*(upper - lower)^2)
        stop("Hypothesized variance v is too large.")

    if(alpha >= 1 | alpha <= 0)
        stop("Please supply a sensible value for alpha.")

    if(iterations < 500)
        warning("Low number of iterations. Results may be inaccurate.")

    ## Computation of sample mean and variance for output
    sample.est <- var(x)
    m <- floor(length(x) / 2)
    x <- (x - lower)/(upper - lower)  ## Normalization so that x in [0,1]
    p <- 2 * v / (upper - lower)^2  ## normalized threshold

    error <- 1
    rejMatrix <- vector(mode = "numeric", length = 0)

    ## deterministic test
    if(alternative == "two.sided")
    {
        ##
        ## first alternative
        ##
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = p, N = m,
                                 alpha = alpha / 2 - epsilon)
        thetaUpper <- minTypeIIError(optimaltypeII[[1]],
                                p = p, N = m, alpha = alpha / 2 - epsilon)
        pseudoalpha <- alpha * thetaUpper$theta / 2
        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     sampleBinomTestnpVar(x, m, p,
                                                          alternative = "greater",
                                                          alpha = pseudoalpha)))
            rejUpper <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rejUpper - thetaUpper$theta)^2)
        }
        rejectionUpper <- ifelse(rejUpper >= thetaUpper$theta, TRUE, FALSE)

        ##
        ## other alternative
        ##
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = 1 - p, N = m,
                                 alpha = alpha / 2 - epsilon)
        thetaLess <- minTypeIIError(optimaltypeII[[1]],
                                    p = 1 - p, N = m,
                                    alpha = alpha / 2 - epsilon)
        pseudoalpha <- alpha * thetaLess$theta / 2
        error <- 1
        rejMatrix <- vector(mode = "numeric", length = 0)

        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     sampleBinomTestnpVar(x, m, p,
                                                          alternative = "less",
                                                          alpha = pseudoalpha)))
            rejLess <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rejLess - thetaLess$theta)^2)
        }
        rejectionLess <- ifelse(rejLess >= thetaLess$theta, TRUE, FALSE)
        rej <- rejUpper + rejLess
        rejection <- ifelse(rejectionUpper + rejectionLess >= 1, TRUE, FALSE)
    }
    else
    {
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1),
                                 p = ifelse(alternative == "greater",
                                            p, 1 - p),
                                 N = m,
                                 alpha = alpha - epsilon)
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = ifelse(alternative == "greater",
                                           p, 1 - p),
                                N = m, alpha = alpha - epsilon)
        pseudoalpha <- alpha * theta$theta

        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     sampleBinomTestnpVar(x, m, p,
                                                          alternative = alternative,
                                                          alpha = pseudoalpha)))
            rej <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rej - theta$theta)^2)
        }
        rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
    }

    if(!is.null(iterations) & length(rejMatrix) < 1000)
        warning("Low number of iterations. Results may be inaccurate.")

    if(length(rejMatrix) >= max.iterations)
        warning(paste("The maximum number of iterations (",
                      format(max.iterations, scientific = FALSE),
                      ") was reached. Rejection may be very sensible to the choice of the parameters.", sep = ""))

    names(sample.est) <- "variance"
    null.value <- v
    names(null.value) <- "variance"
    bounds <- paste("[", lower, ", ", upper, "]", sep = "")

    ## if rejection in a two.sided setting, we inform the user of the
    ## side of rejection
    if(alternative == "two.sided")
    {
        if(rejection == TRUE)
        {
            alt.hypothesis <- paste("Var(", DNAME, ")",
                                    ifelse(rejectionUpper == TRUE, " > ", " < "),
                                    v, sep = "")
        }
        if(rejectionUpper == TRUE) {
            theta <- thetaUpper
        } else {
            theta <- thetaLess
        }
    }

    structure(list(method = "Nonparametric Variance Test",
                   data.name = DNAME,
                   alternative = alternative,
                   null.hypothesis = null.hypothesis,
                   alt.hypothesis = alt.hypothesis,
                   estimate = sample.est,
                   probrej = rej,
                   rejection = rejection,
                   alpha = alpha,
                   theta = theta$theta,
                   d.alternative = 0.5 * optimaltypeII$root * (upper - lower)^2,
                   typeIIerror = theta$typeII,
                   iterations = length(rejMatrix),
                   pseudoalpha = pseudoalpha,
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
}


sampleBinomTestnpVar <- function(x, m, p, alternative, alpha)
{
    x <- sample(x)

    ## transformation into sample in [0,1] that has mean equal to 1/2 +
    ## Var(X)
    ## subtract the odd from the even indexed values, to the power 2
    x.folded <- (x[c(1:m)*2] - x[(c(1:m)*2 - 1)])^2

    ## Random transformation of [0,1] data into {0,p,1} data,
    ## later only use {0,1}
    q <- runif(m)

    ## number of 0 values in transformed data
    zeros <- sum(x.folded - p > (q*(1 - p)))
    ## number of 1 values in transformed data
    ones <- sum(x.folded < q*p)

    ## Evaluation of randomized binomial test, see if the number of zeros
    ## relative to (ones+zeros) is significantly below p

    k <- switch(alternative,
                less = 0:zeros,
                greater = zeros:(ones+zeros))

    ## prob <- sum((p^k)*((1-p)^(ones + zeros - k))*choose((ones +
    ## zeros), k)) ## inefficient
    prob <- sum(dbinom(k, ones + zeros, p))
    ## prob <- pbinom(zeros, ones + zeros, p,
    ##                lower.tail = ifelse(alternative == "greater",
    ##                  TRUE, FALSE)) ## not exact

    if(prob <= alpha) ## reject with probability 1
    {
        res <- 1
    }
    else
    {
        h <- dbinom(zeros, zeros + ones, p) ## more efficient
        if (prob <= alpha + h) ##(reject with positive probability)
        {
            res <- ((alpha - prob + h) / h)
        }
        else
        {
            res <- 0
        }
    }
    return(res)
}
