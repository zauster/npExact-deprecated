##' A test for comparing the means of two bounded random variables given two
##' independent samples
##' 
##' This test requires that the user knows upper and lower bounds before
##' gathering the data such that the properties of the data generating process
##' imply that all observations will be within these bounds. The data input
##' consists of a sequence of independent observations for each random
##' variable, the two sequences being generated independently. No further
##' distributional assumptions are made.
##' 
##' This is a test of the null hypothesis: \eqn{H_0: E(X_1) \le E(X_2)} against
##' \eqn{H_1: E(X_1) > E(X_2)}.
##' 
##' This test uses the known bounds of the variables to transform the data into
##' [0, 1]. Then a random transformation is used to turn the data into
##' binary-valued variables. On this variables the exact Fischer-Tocher Test
##' with level \code{pseudoalpha} is performed and the result recorded. The
##' random transformation and the test are then repeated \code{iterations}
##' times. If the average rejection probability \code{probrej} of the
##' iterations is at least \code{theta}, then the null hypothesis is rejected.
##' If however \code{probrej} is too close to the threshold \code{theta} then
##' the number of iterations is increased. The algorithm keeps increasing the
##' number of iterations until the bound on the mistake involved by running
##' these iterations is below \code{epsilon}. This error epsilon is
##' incorporated into the overall level \code{alpha} in order to maintain that
##' the test is exact.
##' 
##' \code{theta} is found in an optimization procedure. \code{theta} is chosen
##' as to bring the type II error to 0.5. Please see the cited paper below for
##' further information.
##' 
##' @param x1,x2 the (non-empty) numerical data vectors which contain the
##' variables to be tested.
##' @param lower,upper the theoretical lower and upper bounds on the data
##' outcomes known ex-ante before gathering the data.
##' @param iterations the number of iterations used, should not be changed if
##' the exact solution should be derived.
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
##' the alternative hypothesis.  } \item{estimate}{ the sample means of the two
##' variables.  } \item{probrej}{ numerical estimate of the rejection
##' probability of the randomized test, derived by taking an average of
##' \code{iterations} realizations of the rejection probability.  }
##' \item{bounds}{ the lower and upper bounds of the variables.  }
##' \item{null.value}{ the specified hypothesized value of the correlation
##' between the variables.  } \item{alpha}{ the type I error.  } \item{theta}{
##' the parameter that minimizes the type II error.  } \item{pseudoalpha}{
##' \code{theta}*\code{alpha}, this is the level used when calculating the
##' average rejection probability during the iterations } \item{rejection}{
##' logical indicator for whether or not the null hypothesis can be rejected }
##' \item{iterations}{ the number of iterations that were performed.  }
##' @author Karl Schlag, Christian Pechhacker, Peter Saffert and Oliver Reiter
##' @seealso
##' \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
##' @references Karl Schlag (2008), A New Method for Constructing Exact Tests
##' without Making any Assumptions. Available at
##' \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}.
##' @keywords unpaired data mean test
##' @examples
##' 
##' ## test whether countries with french origin score lower than
##' ## countries with no french origin
##' data(french)
##' npMeanUnpaired(french[,1], french[,2], alternative = "less", ignoreNA =
##' TRUE)
##' 
##' ## test whether American tend to be more generous than Isrealis
##' ## in a round of the Ultimatum game
##' data(bargaining)
##' npMeanUnpaired(bargaining$US, bargaining$IS, lower = 0, upper = 10, ignoreNA = TRUE)
##' 
##' @export npMeanUnpaired
npMeanUnpaired <- function(x1, x2,
                           lower = 0, upper = 1,
                           iterations = 5000,
                           alpha = 0.05,
                           alternative = "two.sided",
                           epsilon = 1 * 10^(-6),
                           ignoreNA = FALSE,
                           max.iterations = 100000)
{
    method <- "Nonparametric Mean Test for unpaired variables"

    null.value <- 0
    names(null.value) <- "E[x2] - E[x1]" ## "mean difference"

    names.x1 <- deparse(substitute(x1))
    names.x2 <- deparse(substitute(x2))

    DNAME <- paste(names.x1, "and", names.x2)

    null.hypothesis <- paste("E(", names.x1, ")",
                             ifelse(alternative == "less", " >= ",
                             ifelse(alternative == "greater", " <= ",
                                    " = ")),
                             "E(", names.x2, ")", sep = "")

    alt.hypothesis <- paste("E(", names.x1, ")",
                            ifelse(alternative == "less", " < ",
                            ifelse(alternative == "greater", " > ",
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

    if(ignoreNA == TRUE)
    {
        x1 <- x1[!is.na(x1)]
        x2 <- x2[!is.na(x2)]
    }
    else if(any(is.na(c(x1, x2))) == TRUE)
    {
        stop("The data contains NA's!")
    }

    if(min(x1, x2) < lower | max(x1, x2) > upper)
        stop("Some values are out of bounds!")

    if(alternative != "two.sided" & alternative != "greater" & alternative != "less")
        stop("Please specify which alternative hypothesis you want to test for: 'greater', 'less' or 'two.sided'")

    if(alpha >= 1 | alpha <= 0)
        stop("Please supply a sensible value for alpha.")


    sample.est <- c(mean(x1), mean(x2))
    names(sample.est) <- c(paste("mean(", names.x1, ")", sep = ""),
                           paste("mean(", names.x2, ")", sep = ""))

    ## standardize variables
    ## d <- d/(upper - lower)
    x1 <- (x1 - lower)/(upper - lower)
    x2 <- (x2 - lower)/(upper - lower)

    ## x1 <- as.matrix(x1)
    ## x2 <- as.matrix(x2)

    ## define local variables
    n1 <- length(x1)
    n2 <- length(x2)
    min.length <- min(n1, n2)

    error <- 1
    rejMatrix <- vector(mode = "numeric", length = 0)

    if(alternative == "two.sided")
    {
        ##
        ## alternative "greater" at alpha / 2
        ##
        resultsGreater <- doTwoVariablesTest(alpha = alpha / 2,
                                             epsilon = epsilon,
                                             iterations = iterations,
                                             max.iterations = max.iterations,
                                             testFunction = randomTest,
                                             x1 = x1, x2 = x2,
                                             n1 = n1, n2 = n2)

        ##
        ## alternative "less"
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
                                          testFunction = randomTest,
                                          x1 = x1, x2 = x2,
                                          n1 = n1, n2 = n2)


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
                                      testFunction = randomTest,
                                      x1 = x1, x2 = x2,
                                      n1 = n1, n2 = n2)

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
                   d.alternative = results[["d.alternative"]],
                   typeIIerror = results[["typeIIerror"]],
                   iterations = results[["iterations.taken"]],
                   pseudoalpha = results[["pseudoalpha"]],
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
} ## end of npMeanUnpaired


randomTest <- function(x1, x2, pseudoalpha, dots)
{
    n1 <- dots[["n1"]]
    n2 <- dots[["n2"]]
    
    s1 <- sum(x1 >= runif(n1))
    s2 <- sum(x2 >= runif(n2))
    s3 <- s2 + s1
    k <- max(0, s3 - n2):(s1 - 1)
    prob <- 0
    if (s1 >= (1 + k[1]))
    {
        prob <- sum(choose(n1,
                           k) * choose(n2,
                                       s3 - k)/choose(n1 + n2,
                                                      s3))
        ## h.alt <- phyper(s1 - 1, n1, n2, A)
    }

    res <- 0
    if (prob <= pseudoalpha)
    {
        ## h2 <- prob + choose(n1, s1) * choose(n2, s3 - s1)/choose(n1 +
        ## n2, s3)
        h2 <- prob + dhyper(s1, n1, n2, s3)
        if (h2 <= pseudoalpha)
        {
            res <- 1
        }
        else
        {
            res <- ((pseudoalpha - prob)/(h2 - prob))
        }
    }
    return(res)
}


########################################
## Theta function new
########################################


## calculates pvalue of Fisher's test
pvalueFisher <- function(n1, n2, s1, s2)
{
    ## if( s1 == -1 | s2 > n2)
    ##   return(0)
    ## else
    ##   {
    ##     phyper(s2 - 1, n2, n1, s1 + s2,
    ##                  lower.tail = FALSE)
    ##   }

    ## try to vectorize it -> seems to work
    ifelse(s1 == -1 | s2 > n2, 0, phyper(s2 - 1, n2, n1, s1 + s2,
                                         lower.tail = FALSE))
}

## ## calculates typeII error for given y1, y2 (and of course n1,n2, alpha)
## maxTypeII <- function(y1, d, n1, n2, y2 = y1 + d,
##                       alpha = alpha, theta = 0.2)
## {
##   pseudoalpha <- theta * alpha
##   ## exmat <- matrix(nrow = n1 + 1, ncol = n2 + 1)
##   res <- 0
##   for(s1 in 0:n1)
##     {
##   ##     for(s2 in 0:n2)
##   ##       {
##   ##         t1 <- pvalueFisher(n1, n2, s1, s2)
##   ##         t2 <- pvalueFisher(n1, n2, s1 - 1, s2 + 1)

##   ##         if( t2 >= pseudoalpha)
##   ##           {
##   ##             ## in this case, pr = zero
##   ##             ## so we can skip the calculation of the term
##   ##             exmat[s1 + 1, s2 + 1] <- 0
##   ##           }
##   ##         else
##   ##           {
##   ##             if( t1  > pseudoalpha & pseudoalpha > t2)
##   ##               {
##   ##                 ## pr <- (pseudoalpha - t2) / ((choose(n1, s1) *
##   ##                 ## choose(n2, s2))/ choose(n1 + n2, s1 + s2))
##   ##                 pr <- (pseudoalpha - t2) / dhyper(s1, n1, n2, s1 + s2)
##   ##               }
##   ##             else
##   ##               {
##   ##                 if(t1 <= pseudoalpha) pr <- 1
##   ##               }
##   ##             exmat[s1 + 1,
##   ##                   s2 + 1] <- dbinom(s1, n1, y1) * dbinom(s2, n2, y2) * pr
##   ##           }
##   ##       }
##   ## }
##       ## now instead of the second for clause -> vectorized if-clauses!
##       t1 <- pvalueFisher(n1, n2, s1, 0:n2)
##       t2 <- pvalueFisher(n1, n2, s1 - 1, 1:(n2 + 1))
##       res <- res + sum(ifelse(t2 >= pseudoalpha, 0,
##                     ifelse(t1 > pseudoalpha & pseudoalpha > t2,
##                            dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * (pseudoalpha - t2) / dhyper(s1, n1, n2, s1 + 0:n2),
##                            dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * 1)))
##     }
##   type2 <- (1 - res) / (1 - theta)
##   ## type2 <- (1 - sum(exmat)) / (1 - theta)
##   return(min(type2, 1))
## }

maxTypeII <- function(y1, d, n1, n2, y2 = y1 + d,
                      alpha = 0.05, theta = 0.2)
{
    pseudoalpha <- theta * alpha
    res <- 0
    f <- function(s1, n1, n2, y1, y2)
    {
        t1 <- pvalueFisher(n1, n2, s1, 0:n2)
        t2 <- pvalueFisher(n1, n2, s1 - 1, 1:(n2 + 1))
        res <- sum(ifelse(t2 >= pseudoalpha, 0,
                   ifelse(t1 > pseudoalpha & pseudoalpha > t2,
                          dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * (pseudoalpha - t2) / dhyper(s1, n1, n2, s1 + 0:n2),
                          dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * 1)))
        res
    }
    res <- sum(sapply(0:n1,
                      f, n1, n2, y1, y2))
    type2 <- (1 - res)/(1 - theta)
    return(min(type2, 1))
}

## same as function typeII error, only order of inputs changed,
## so that function "optimize" can be used
minTypeII <- function(theta, y1, y2, n1, n2, alpha)
{
    maxTypeII(y1, d = y2 - y1, n1, n2, y2,
              alpha = alpha, theta = theta)
}

### calculate theta
optimizeTheta <- function(n1, n2, diff, alpha = alpha)
{
    ## STEP 1)  maximize typeII error over y1, y2
    ## cat("\nmax ")
    maxexpect <- optimize(maxTypeII, c(0, 1 - diff),
                          ## tol = .Machine$double.eps^0.5,
                          d = diff, n1 = n1, n2 = n2,
                          alpha = alpha, maximum = T)
    e1opt <- maxexpect$maximum
    e2opt <- e1opt + diff
    ## cat(e1opt, " ")
    ## cat(e2opt, " ")

    ## STEP 2)  minimize typeII error over theta
    ## cat("min ")
    thetaval <- optimize(minTypeII, c(0,1),
                         ## tol = .Machine$double.eps^0.5,
                         n1 = n1, n2 = n2,
                         y1 = e1opt, y2 = e2opt, alpha = alpha)

    ## if(thetaval$objective == 1)
    ##   stop("TypeII error = 1. Increase difference d")

    ## cat(thetaval$minimum, " ", thetaval$objective)

    return(list(typeII = thetaval$objective,
                theta = thetaval$minimum))
}

npMeanUnpairedminTypeIIErrorWrapper <- function(d, n1, n2, alpha,
                                                typeIIgoal = 0.5)
{
    (optimizeTheta(n1, n2, d, alpha)$typeII - typeIIgoal)^2
}

## optimaltypeII <- optimize(npMeanUnpairedminTypeIIErrorWrapper,
##                           c(0, 1), n1 = 25, n2 = 29, alpha = .05)
## theta <- optimizeTheta(n1, n2, optimaltypeII$minimum, alpha)
## optimaltypeII <- uniroot(npMeanUnpairedminTypeIIErrorWrapper,
##                          c(0, 29), n1 = 25, n2 = 29,
##                          alpha = alpha)
## theta <- optimizeTheta(n1, n2,
##                        optimaltypeII[[1]], alpha = alpha)
## thetause <- theta(N1, N2 , diff = d,
##                   alpha = alpha)$theta
## pseudoalpha <- alpha * thetause
