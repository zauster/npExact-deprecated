##' A test for the mean difference between two bounded random variables given
##' matched pairs.
##' 
##' This test requires that the user knows bounds before gathering the data
##' such that the properties of the data generating process imply that all
##' observations will be within these bounds. The data input consists of pairs
##' of observations, each pair consisting of an observation of each random
##' variable, different pairs being independently generated. No further
##' distributional assumptions are made.
##' 
##' Under alternative = "greater", it is a test of the null hypothesis
##' \eqn{H_0: E(x_1) \le E(x_2)} against the alternative hypothesis \eqn{H_1:
##' E(x_1) > E(x_2)}.
##' 
##' 
##' This test uses the known bounds of the variables to transform the data into
##' [0, 1]. Then a random transformation is used to turn the data into
##' binary-valued variables. On this variables the exact McNemar Test with
##' level \code{pseudoalpha} is performed and the result recorded. The random
##' transformation and the test are then repeated \code{iterations} times. If
##' the average rejection probability \code{probrej} of the iterations is at
##' least \code{theta}, then the null hypothesis is rejected. If however
##' \code{probrej} is too close to the threshold \code{theta} then the number
##' of iterations is increased. The algorithm keeps increasing the number of
##' iterations until the bound on the mistake involved by running these
##' iterations is below \code{epsilon}. This error epsilon is incorporated into
##' the overall level \code{alpha} in order to maintain that the test is exact.
##' 
##' \code{theta} (and a value \code{mu} of the difference between the two means
##' in the set of the alternative hypothesis) is found in an optimization
##' procedure. \code{theta} and \code{mu} are chosen as to maximize the set of
##' data generating processes belonging to the alternative hypothesis that
##' yield type II error probability below 0.5. Please see the cited paper below
##' for further information.
##' 
##' @param x1,x2 the (non-empty) numerical data vectors which contain the
##' variables to be tested. The first values of the vectors are assumed to be
##' the first matched pair of observations, the second values the second
##' matched pair and so on.
##' @param lower,upper the theoretical lower and upper bounds on the data
##' outcomes known ex-ante before gathering the data.
##' @param alpha the type I error.
##' @param alternative a character string describing the alternative
##' hypothesis, can take values "greater", "less" or "two.sided".
##' @param iterations the number of iterations used, should not be changed if
##' the exact solution should be derived
##' @param epsilon the tolerance in terms of probability of the Monte Carlo
##' simulations.
##' @param max.iterations the maximum number of iterations that should be
##' carried out. This number could be increased to achieve greater accuracy in
##' cases where the difference between the threshold probability and theta is
##' small. Default: \code{10000}
##' @return A list with class "nphtest" containing the following components:
##' 
##' \item{method}{ a character string indicating the name and type of the test
##' that was performed.  } \item{data.name}{ a character string giving the
##' name(s) of the data.  } \item{alternative}{ a character string describing
##' the alternative hypothesis.  } \item{estimate}{ the sample means of the
##' given data.  } \item{probrej}{ numerical estimate of the rejection
##' probability of the randomized test, derived by taking an average of
##' \code{iterations} realizations of the rejection probability.  }
##' \item{bounds}{ the lower and upper bounds of the variables.  }
##' \item{null.value}{ the specified hypothesized value of the difference of
##' the variable means.  } \item{alpha}{ the type I error.  } \item{theta}{ the
##' parameter that minimizes the type II error.  } \item{pseudoalpha}{
##' \code{theta}*\code{alpha}, this is the level used when calculating the
##' average rejection probability during the iterations.  } \item{rejection}{
##' logical indicator for whether or not the null hypothesis can be rejected.
##' } \item{iterations}{ the number of iterations that were performed.  }
##' @author Karl Schlag, Christian Pechhacker and Oliver Reiter
##' @seealso
##' \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
##' @references Schlag, Karl H. 2008, A New Method for Constructing Exact Tests
##' without Making any Assumptions, Department of Economics and Business
##' Working Paper 1109, Universitat Pompeu Fabra. Available at
##' \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}.
##' @keywords pairwise mean test
##' @examples
##' 
##' ## test whether pain after the surgery is less than before the surgery
##' data(pain)
##' npMeanPaired(pain$before, pain$after, lower = 0, upper = 100)
##' 
##' ## when the computer was used in the surgery
##' before_pc <- pain[pain$pc == 1, "before"]
##' after_pc <- pain[pain$pc == 1, "after"]
##' npMeanPaired(before_pc, after_pc, lower = 0, upper = 100)
##' 
##' ## test whether uncertainty decreased from the first to the second round
##' data(uncertainty)
##' npMeanPaired(uncertainty$w1, uncertainty$w2, upper = 60) ## or
##' with(uncertainty, npMeanPaired(w1, w2, upper = 60))
##' 
##' @export npMeanPaired
npMeanPaired <- function(x1, x2, lower = 0, upper = 1, ## d = 0,
                         alpha = 0.05,
                         alternative = "two.sided",
                         epsilon = 1 * 10^(-6),
                         iterations = 5000,
                         max.iterations = 100000)
{
    method <- "Nonparametric Mean Test for Matched Pairs"
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

    ## if x1 (or x2) is a 1-column data.frame, convert it to a vector
    if(is.data.frame(x1)) {
        if(dim(x1)[2] == 1) {
            x1 <- x1[, 1]
        }
    }
    if(is.data.frame(x2)) {
        if(dim(x2)[2] == 1) {
            x2 <- x2[, 1]
        }
    }
    
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

    ## Sample estimates
    sample.est <- c(mean(x1), mean(x2))
    names(sample.est) <- c(paste("mean(", names.x1, ")", sep = ""),
                           paste("mean(", names.x2, ")", sep = ""))

    ## Null value, for prettier output later
    null.value <- 0
    names(null.value) <- "E[x2] - E[x1]"

    ## Normalize vectors to [0,1]
    x1 <- (x1 - lower)/(upper - lower)
    x2 <- (x2 - lower)/(upper - lower)

    if(alternative == "two.sided")
    {
        ##
        ## alternative "greater" at alpha / 2
        ##
        resultsGreater <- doTwoVariablesTest(alpha = alpha / 2,
                                             epsilon = epsilon,
                                             iterations = iterations,
                                             max.iterations = max.iterations,
                                             testFunction = McNemarTestRandom,
                                             x1 = x1, x2 = x2,
                                             p = 0.5, n = n)

        ##
        ## alternative "less" at alpha / 2
        ##
        x1 <- 1 - x1
        x2 <- 1 - x2

        resultsLess <- doTwoVariablesTest(alpha = alpha / 2,
                                          epsilon = epsilon,
                                          theta = resultsGreater[["theta"]],
                                          typeII = resultsGreater[["typeIIerror"]],
                                          d.alternative = resultsGreater[["d.alternative"]],
                                          iterations = iterations,
                                          max.iterations = max.iterations,
                                          testFunction = McNemarTestRandom,
                                          x1 = x1, x2 = x2,
                                          p = 0.5, n = n)


        ## "greater" rejects 
        if(resultsGreater[["rejection"]] == TRUE) {
            results <- resultsGreater
            theta <- resultsGreater[["theta"]]
        }
        ## "less" rejects
        else if(resultsLess[["rejection"]] == TRUE) {
            results <- resultsLess
            theta <- resultsLess[["theta"]]
        }
        ## none rejects:
        ## we take the one that is more likely to reject
        else {
            if((sample.est[1] - sample.est[2] > 0) & !is.null(resultsGreater[["theta"]])) {
                results <- resultsGreater
                theta <- resultsGreater[["theta"]]
            }
            else if((sample.est[1] - sample.est[2] < 0) & !is.null(resultsLess[["theta"]])) {
                results <- resultsLess
                theta <- resultsLess[["theta"]]
            } else {
                results <- resultsGreater                
                theta <- resultsGreater[["theta"]]
            }
        }

        results <- mergeTwoResultSets(results, resultsGreater, resultsLess)
         
        ## if rejection in a two.sided setting, we inform the user of the
        ## side of rejection
        if(results[["rejection"]] == TRUE)
        {
            alt.hypothesis <- paste("E(", names.x1, ")",
                                    ifelse(resultsGreater[["rejection"]] == TRUE, " < ", " > "),
                                    "E(", names.x2, ")", sep = "")
        }
    }
    else
    {
        if(alternative == "greater")
        {
            x1 <- 1 - x1
            x2 <- 1 - x2
        }

        results <- doTwoVariablesTest(alpha = alpha,
                                      epsilon = epsilon,
                                      iterations = iterations,
                                      max.iterations = max.iterations,
                                      testFunction = McNemarTestRandom,
                                      x1 = x1, x2 = x2,
                                      p = 0.5, n = n)

        theta <- results[["theta"]]
        if(alternative == "less" & !is.null(results[["d.alternative"]])) {
            results[["d.alternative"]] <- 1 - results[["d.alternative"]]
        }
    }

    if(!is.null(iterations) & results[["iterations.taken"]] < 1000)
        warning("Low number of iterations. Results may be inaccurate.")

    if(results[["iterations.taken"]] >= max.iterations)
        warning(paste("The maximum number of iterations (",
                      format(max.iterations, scientific = FALSE),
                      ") was reached. Rejection may be very sensible to the choice of the parameters.", sep = ""))

    bounds <- paste("[", lower, ", ", upper, "]", sep = "")

    structure(list(method = method,
                   data.name = DNAME,
                   alternative = alternative,
                   null.hypothesis = null.hypothesis,
                   alt.hypothesis = alt.hypothesis,
                   estimate = sample.est,
                   probrej = results[["probrej"]],
                   rejection = results[["rejection"]],
                   mc.error = results[["mc.error"]],
                   alpha = alpha,
                   theta = theta,
                   thetaValue = results[["theta"]],
                   d.alternative = (results[["d.alternative"]] - 0.5) * 2 * (upper - lower),
                   typeIIerror = results[["typeIIerror"]],
                   iterations = results[["iterations.taken"]],
                   pseudoalpha = results[["pseudoalpha"]],
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
} ## end of npMeanPaired


## performs the randomized McNemar test
## x1, x2 are bounden real-valued vectors 
McNemarTestRandom <- function(x1, x2, pseudoalpha, dots)
{
    ## b1, b2 are binary-valued vectors of equal length
    ## n10, n01 ... counts of (1,0) and (0,1) respectively.
    ##              (1, 1) and (0, 0) are dropped

    n <- dots[["n"]]
    
    ## returns either 1 (rejection), p (prob of rejection) or 0 (no
    ## rejection)

    b1 <- runif(n) < x1
    b2 <- runif(n) < x2

    n10 <- sum(b1 > b2)
    n01 <- sum(b1 < b2)

    k <- n01:(n10 + n01)
    prob <- sum(choose(n10 + n01, k)/(2^(n10 + n01)))

    res <- 0
    if(prob <= pseudoalpha)
    {
        res <- 1
    }
    else
    {
        h <- choose(n10 + n01, n01)/(2^(n10 + n01))
        if(prob <= (pseudoalpha + h))
        {
            res <- (pseudoalpha - prob + h)/h
        }
    }
    return(res)
}
