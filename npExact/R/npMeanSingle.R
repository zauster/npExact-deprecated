##' A test for the mean of a bounded random variable based on a single sample
##' of iid observations.
##' 
##' This test requires that the user knows upper and lower bounds before
##' gathering the data such that the properties of the data generating process
##' imply that all observations will be within these bounds. The data input
##' consists of a sequence of observations, each being an independent
##' realization of the random variable. No further distributional assumptions
##' are made.
##' 
##' For any \eqn{\mu} that lies between the two bounds, under alternative =
##' "greater", it is a test of the null hypothesis \eqn{H_0 : E(X) \le \mu}
##' against the alternative hypothesis \eqn{H_1 : E(X) > \mu}.
##' 
##' Using the known bounds, the data is transformed to lie in [0, 1] using an
##' affine transformation. Then the data is randomly transformed into a new
##' data set that has values 0, \code{mu} and 1 using a mean preserving
##' transformation. The exact randomized binomial test is then used to
##' calculate the rejection probability of this under new data when level is
##' \code{theta}*\code{alpha}. This random transformation is repeated
##' \code{iterations} times. If the average rejection probability is greater
##' than theta, one can reject the null hypothesis. If however the average
##' rejection probability is too close to theta then the iterations are
##' continued. The values of \code{theta} and a value of \code{mu} in the
##' alternative hypothesis is found in an optimization procedure to maximize
##' the set of parameters in the alternative hypothesis under which the type II
##' error probability is below 0.5. Please see the cited paper below for
##' further information.
##' 
##' @param x a (non-empty) numeric vector of data values.
##' @param mu threshold value for the null hypothesis.
##' @param lower,upper the theoretical lower and upper bounds on the data
##' outcomes known ex-ante before gathering the data.
##' @param iterations the number of iterations used, should not be changed if
##' the exact solution should be derived
##' @param alpha the type I error.
##' @param alternative a character string describing the alternative
##' hypothesis, can take values "greater", "less" or "two.sided".
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
##' between the variables.  } \item{alpha}{ the type I error } \item{theta}{
##' the parameter that minimizes the type II error.  } \item{pseudoalpha}{
##' \code{theta}*\code{alpha}, this is the level used when calculating the
##' average rejection probability during the iterations.  } \item{rejection}{
##' logical indicator for whether or not the null hypothesis can be rejected.
##' } \item{iterations}{ the number of iterations that were performed.  }
##' @author Karl Schlag, Peter Saffert and Oliver Reiter
##' @seealso
##' \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
##' @references Schlag, Karl H. 2008, A New Method for Constructing Exact Tests
##' without Making any Assumptions, Department of Economics and Business
##' Working Paper 1109, Universitat Pompeu Fabra. Available at
##' \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}.
##' @keywords single sample mean test
##' @examples
##' 
##' ## test whether Americans gave more than 5 dollars in a round of
##' ## the Ultimatum game
##' data(bargaining)
##' us_offers <- bargaining$US
##' npMeanSingle(us_offers, mu = 5, lower = 0, upper = 10, alternative =
##' "greater", ignoreNA = TRUE) ## no rejection
##' 
##' ## test if the decrease in pain before and after the surgery is smaller
##' ## than 50
##' data(pain)
##' pain$decrease <- with(pain, before - after)
##' without_pc <- pain[pain$pc == 0, "decrease"]
##' npMeanSingle(without_pc, mu = 50, lower = 0, upper = 100,
##' alternative = "less")
##' 
##' @export npMeanSingle
npMeanSingle <- function(x, mu,
                         lower = 0, upper = 1,
                         alternative = "two.sided",
                         iterations = 5000, alpha = 0.05,
                         epsilon = 1 * 10^(-6),
                         ignoreNA = FALSE,
                         max.iterations = 100000)
{
    method <- "Nonparametric Single Mean Test"
    DNAME <- deparse(substitute(x))

    ## if x is a 1-column data.frame, convert it to a vector
    if(is.data.frame(x)) {
        if(dim(x)[2] == 1) {
            x <- x[, 1]
        }
    }

    x <- as.vector(x)
    sample.est <- mean(x)
    names(sample.est) <- "mean"

    null.value <- mu
    names(null.value) <- "mean"
    
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

    if(ignoreNA == TRUE)
    {
        x <- x[!is.na(x)]
    }
    else if(any(is.na(x)) == TRUE)
    {
        stop("The data contains NA's!")
    }

    ## warnings
    if (min(x) < lower | max(x) > upper)
        stop("Some values are out of bounds (or NA)!")

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

    error <- 1
    rejMatrix <- vector(mode = "numeric", length = 0)

    if(alternative == "two.sided")
    {
        ##
        ## first the upper side, alternative = "greater"
        ##
        res <- try(optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                            c(0, 1), p = p, N = n,
                                            alpha = alpha / 2 - epsilon),
                   silent = TRUE)
        if(inherits(res, "try-error")) {
            ## pick up an error in the theta calculation
            return(structure(list(method = method,
                                  data.name = DNAME,
                                  alternative = alternative,
                                  null.hypothesis = null.hypothesis,
                                  alt.hypothesis = alt.hypothesis,
                                  estimate = sample.est,
                                  probrej = NULL,
                                  rejection = FALSE,
                                  alpha = NULL,
                                  theta = NULL,
                                  d.alternative = NULL,
                                  typeIIerror = NULL,
                                  iterations = NULL,
                                  pseudoalpha = NULL,
                                  bounds = NULL,
                                  null.value = null.value),
                             class = "nphtest"))
        }

        thetaUpper <- minTypeIIError(optimaltypeII[[1]],
                                     p = p, N = n, alpha = alpha / 2 - epsilon)
        pseudoalpha <- (alpha/2) * thetaUpper$theta

        while(error > epsilon & length(rejMatrix) <= max.iterations) {
        rejMatrix <- c(rejMatrix,
                       replicate(iterations,
                                 transBinomTest(x, p, xp, n,
                                                pseudoalpha)))
        rejUpper <- mean(rejMatrix)
        error <- exp(-2 * length(rejMatrix) * (rejUpper - thetaUpper$theta)^2)
        }
        
        rejectionUpper <- ifelse(rejUpper >= thetaUpper$theta, TRUE, FALSE)
        iterations.taken <- length(rejMatrix)


        ##
        ## secondly the lower side, alternative = "less"
        ##
        x <- 1 - x
        p <- 1 - p
        xp <- x - p
        
        res <- try(optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                            c(0, 1), p = p, N = n,
                                            alpha = alpha / 2 - epsilon),
                   silent = TRUE)
        if(inherits(res, "try-error")) {
            ## pick up an error in the theta calculation
            return(structure(list(method = method,
                                  data.name = DNAME,
                                  alternative = alternative,
                                  null.hypothesis = null.hypothesis,
                                  alt.hypothesis = alt.hypothesis,
                                  estimate = sample.est,
                                  probrej = NULL,
                                  rejection = FALSE,
                                  alpha = NULL,
                                  theta = NULL,
                                  d.alternative = NULL,
                                  typeIIerror = NULL,
                                  iterations = NULL,
                                  pseudoalpha = NULL,
                                  bounds = NULL,
                                  null.value = null.value),
                             class = "nphtest"))
        }
        
        thetaLess <- minTypeIIError(optimaltypeII[[1]],
                                    p = p, N = n,
                                    alpha = alpha / 2 - epsilon)
        pseudoalpha <- (alpha/2) * thetaLess$theta

        error <- 1
        rejMatrix <- vector(mode = "numeric", length = 0)
        while(error > epsilon & length(rejMatrix) <= max.iterations) {
        rejMatrix <- c(rejMatrix,
                       replicate(iterations,
                                 transBinomTest(x, p, xp, n,
                                                pseudoalpha)))
        rejLess <- mean(rejMatrix)
        error <- exp(-2 * length(rejMatrix) * (rejLess - thetaLess$theta)^2)
        }
        
        rejectionLess <- ifelse(rejLess >= thetaLess$theta, TRUE, FALSE)
        iterations.taken <- max(length(rejMatrix), iterations.taken)

        rej <- rejUpper + rejLess
        rejection <- ifelse(rejectionUpper + rejectionLess >= 1, TRUE, FALSE)
    }
    else ## alternative == "greater" => default
    {
        if(alternative == "less")
        {
            x <- 1 - x
            p <- 1 - p
            xp <- x - p
        }

        res <- try(optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                            c(0, 1), p = p, N = n,
                                            alpha = alpha - epsilon),
                   silent = TRUE)
        if(inherits(res, "try-error")) {
            ## pick up an error in the theta calculation
            return(structure(list(method = method,
                                  data.name = DNAME,
                                  alternative = alternative,
                                  null.hypothesis = null.hypothesis,
                                  alt.hypothesis = alt.hypothesis,
                                  estimate = sample.est,
                                  probrej = NULL,
                                  rejection = FALSE,
                                  alpha = NULL,
                                  theta = NULL,
                                  d.alternative = NULL,
                                  typeIIerror = NULL,
                                  iterations = NULL,
                                  pseudoalpha = NULL,
                                  bounds = NULL,
                                  null.value = null.value),
                             class = "nphtest"))
        }
        
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = p, N = n, alpha = alpha - epsilon)
        pseudoalpha <- alpha * theta$theta

        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     transBinomTest(x, p, xp, n,
                                                    pseudoalpha)))
            rej <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rej - theta$theta)^2)
        }
        rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
        iterations.taken <- length(rejMatrix)
    }

    if(!is.null(iterations) & iterations.taken < 1000)
        warning("Low number of iterations. Results may be inaccurate.")

    if(iterations.taken >= max.iterations)
        warning(paste("The maximum number of iterations (",
                      format(max.iterations, scientific = FALSE),
                      ") was reached. Rejection may be very sensible to the choice of the parameters.", sep = ""))

    ## if rejection in a two.sided setting, we inform the user of the
    ## side of rejection
    if(alternative == "two.sided")
    {
        if(rejection == TRUE)
        {
            alt.hypothesis <- paste("E(", DNAME, ")",
                                    ifelse(rejectionUpper == TRUE,
                                           " > ", " < "),
                                    mu, sep = "")
        }
        if(rejectionUpper == TRUE) {
            theta <- thetaUpper
        } else {
            theta <- thetaLess
        }
    }

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
                   d.alternative = optimaltypeII$root,
                   typeIIerror = theta$typeII,
                   iterations = iterations.taken,
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
        h2 <- dbinom(ones, zeros + ones, p)
        if (res.binomtest <= (pseudoalpha + h2))
        {
            res <- ((pseudoalpha - res.binomtest + h2)/h2)
        }
    }
    res
}
