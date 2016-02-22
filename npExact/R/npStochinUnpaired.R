##' A test of a stochastic inequality given two independent samples
##' 
##' The data input consists of a sequence of independent realizations
##' observations of each random variable, observations of the different
##' sequences also being independent.
##' 
##' Given \eqn{-1 < d < 1} it is a test of the null hypothesis \eqn{H_0 : P(X_2
##' > X_1) \le P(X_2 < X_1) + d} against the alternative hypothesis \eqn{H_1 :
##' P(X_2 > X_1) > P(X_2 < X_1) + d}.
##' 
##' 
##' The data is randomly matched into pairs and then treats them as matched
##' pairs. The number of pairs is equal to the number of observations in the
##' smaller sequence. The exact randomized test is then used to determine if
##' sufficiently many occurrences of \eqn{x_2 > x_1} occur when compared to how
##' often \eqn{x_2 < x_1} occurs, using level \code{theta}*\code{alpha}. The
##' matching into pairs is repeated \code{iterations} times. The test gives a
##' rejection of the average rejection probability in these iterations lies
##' above \code{theta}. If the average rejection probability lies too close to
##' theta then the number of iterations is increased.
##' 
##' \code{theta} is determined to maximize the set of differences
##' \eqn{P(X_2>X_1) - P(X_2<X_1)} belonging to the alternative hypothesis in
##' which the type II error probability lies below 0.5. For more details see
##' the paper.
##' 
##' @param x1,x2 the (non-empty) numerical data vectors which contain the
##' variables to be tested.
##' @param d the maximal difference in probabilities assumed \eqn{H_0 : P(X_2 >
##' X_1) - P(X_2 < X_1) <= d}. Default is 0.
##' @param alternative a character string describing the alternative
##' hypothesis. Default is "greater". If "less" is given, \code{x1} and
##' \code{x2} are switched for each other.
##' @param iterations the number of iterations used, should not be changed if
##' the exact solution should be derived.
##' @param alpha the type I error.
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
##' the alternative hypothesis.  } \item{estimate}{ an estimate of \eqn{P(x_2 >
##' x_1) - P(x_2 < x_1)}.  } \item{probrej}{ numerical estimate of the
##' rejection probability of the randomized test, derived by taking an average
##' of \code{iterations} realizations of the rejection probability.  }
##' \item{bounds}{ the lower and upper bounds of the variables.  }
##' \item{null.value}{ the specified hypothesized value of the correlation
##' between the variables.  } \item{alpha}{ the type I error.  } \item{theta}{
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
##' @keywords unpaired data stochastic inequality
##' @examples
##' 
##' data(french)
##' x <- french[, 1]
##' y <- french[, 2]
##' npStochinUnpaired(x, y, ignoreNA = TRUE)
##' 
##' @export npStochinUnpaired
npStochinUnpaired <- function(x1, x2, d = 0,
                              alternative = "two.sided",
                              iterations = 5000, alpha = 0.05,
                              epsilon = 1 * 10^(-6),
                              ignoreNA = FALSE,
                              max.iterations = 100000)
{
    method <- "Nonparametric Test for Stochastic Inequality"
    names.x1 <- deparse(substitute(x1))
    names.x2 <- deparse(substitute(x2))
    DNAME <- paste(names.x1, "and",
                   names.x2)

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

    if(alpha >= 1 | alpha <= 0)
        stop("Please supply a sensible value for alpha.")

    null.value <- d
    names(null.value) <- "relation P(x1 > x2) - P(x1 < x2)"
    
    ## swap variable if alternative is "less"
    if(alternative == "less")
    {
        names.x1.new <- names.x2
        names.x2 <- names.x1
        names.x1 <- names.x1.new
        x1.new <- x2
        x2 <- x1
        x1 <- x1.new

        d <- -d
    }

    ## define local variables
    N1 <- length(x1)
    N2 <- length(x2)
    min.length <- min(N1, N2)
    p <- (1 + d)/2

    ## compute the sample estimate
    if(alternative == "less") {
        count.x1 <- sum(tapply(x1, 1:N1, function(x.i) sum(x.i < x2)))
        count.x2 <- sum(tapply(x1, 1:N1, function(x.i) sum(x.i > x2)))
        sample.est <- (count.x1 - count.x2)/(N1 * N2)
        stochin.parameter <- paste("P(", names.x1, " < ", names.x2, ") - P(",
                                   names.x1, " > ", names.x2, ")",
                                   sep = "")
    } else {
        count.x1 <- sum(tapply(x1, 1:N1, function(x.i) sum(x.i > x2)))
        count.x2 <- sum(tapply(x1, 1:N1, function(x.i) sum(x.i < x2)))
        sample.est <- (count.x1 - count.x2)/(N1 * N2)
        stochin.parameter <- paste("P(", names.x1, " > ", names.x2, ") - P(",
                                   names.x1, " < ", names.x2, ")",
                                   sep = "")
    }

    ## set name of estimate
    names(sample.est) <- stochin.parameter 

    ## null and alternative hypothesis
    null.hypothesis <- paste("SI",
                             ifelse(alternative == "greater", " <= ",
                             ifelse(alternative == "less", " >= ",
                                    " = ")),
                             null.value, sep = "")
    alt.hypothesis <- paste("SI",
                            ifelse(alternative == "greater", " > ",
                            ifelse(alternative == "less", " < ", " != ")),
                            null.value, sep = "")
    
    if(alternative == "two.sided")
    {
        ##
        ## alternative = "greater" at alpha / 2
        ##
        resultsGreater <- doTwoVariablesTest(alpha = alpha / 2,
                                             epsilon = epsilon,
                                             iterations = iterations,
                                             max.iterations = max.iterations,
                                             testFunction = sampleBinomTest,
                                             x1 = x2, x2 = x1,
                                             p = p, n = min.length,
                                             diff = d)
        
        ##
        ## alternative = "less" at alpha / 2
        ##
        d <- -d
        p <- (1 + d)/2

        resultsLess <- doTwoVariablesTest(alpha = alpha / 2,
                                      epsilon = epsilon,
                                      theta = resultsGreater[["theta"]],
                                      typeII = resultsGreater[["typeIIerror"]],
                                      d.alternative = resultsGreater[["d.alternative"]],
                                      iterations = iterations,
                                      max.iterations = max.iterations,
                                      testFunction = sampleBinomTest,
                                      x1 = x1, x2 = x2,
                                      p = p, n = min.length,
                                      diff = d)


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
            if((sample.est < null.value) & !is.null(resultsGreater[["theta"]])) {
                results <- resultsGreater
                theta <- resultsGreater[["theta"]]
            }
            else if((sample.est > null.value) & !is.null(resultsLess[["theta"]])) {
                results <- resultsLess
                theta <- resultsLess[["theta"]]
            } else {
                results <- resultsGreater                
                theta <- resultsGreater[["theta"]]
            }
        }

        results <- mergeTwoResultSets(results, resultsGreater, resultsLess)
        
        if(results[["rejection"]] == TRUE) {
            alt.hypothesis <- paste("SI",
                                    ifelse(resultsGreater[["rejection"]] == TRUE, " > ", " < "),
                                    null.value, sep = "")
        }
        
    }
    else
    {
        results <- doTwoVariablesTest(alpha = alpha,
                                      epsilon = epsilon,
                                      iterations = iterations,
                                      max.iterations = max.iterations,
                                      testFunction = sampleBinomTest,
                                      x1 = x2, x2 = x1,
                                      p = p, n = min.length, diff = d)

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


    structure(list(method = method,
                   data.name = DNAME,
                   alternative = alternative,
                   stochin.parameter = stochin.parameter,
                   stochin.estimate = sample.est,
                   null.hypothesis = null.hypothesis,
                   alt.hypothesis = alt.hypothesis,
                   estimate = NULL,
                   probrej = results[["probrej"]],
                   rejection = results[["rejection"]],
                   mc.error = results[["mc.error"]],
                   alpha = alpha,
                   theta = theta,
                   thetaValue = results[["theta"]],
                   d.alternative = (results[["d.alternative"]] * 2 - 1),
                   typeIIerror = results[["typeIIerror"]],
                   iterations = results[["iterations.taken"]],
                   pseudoalpha = results[["pseudoalpha"]],
                   bounds = NULL,
                   null.value = null.value),
              class = "nphtest")
}


## sampleBinomTest
## executes a binomial test randomly sampled matched pairs

## x1, x2 ... data vectors
## n ... the minimum length of the data vectors

sampleBinomTest <- function(x1, x2, pseudoalpha, dots)
{
    n <- dots[["n"]]
    p <- dots[["p"]]
    d <- dots[["diff"]]
    
    c1 <- sample(x1, n)
    c2 <- sample(x2, n)

    s1 <- sum(c1 > c2) #counts how often c1 > c2
    s2 <- sum(c1 < c2) #counts how often c1 < c2

    if((s1 + s2) != n)
    {
        if(d > 0)
        {
            q <- runif(sum(c1 == c2)) #vector with draws from uniform
                                        #distribution of length=number of
                                        #times elements of c1=c2
            s1 <- s1 + sum(q < (d/(1 + d)))
        }
        else
        {
            if(d < 0)
            {
                q <- runif(sum(c1 == c2))  #vector with draws from
                                        #uniform distribution of
                                        #length=number of times
                                        #elements of c1=c2
                s1 <- s1 + sum(q < (-1) * (d/(1 - d))) ## right?
            }
        }
    }
    
    prob <- sum(dbinom(s2:(s1 + s2), (s1 + s2), p)) ## or
    ## prob <- 1 - pbinom(s2 - 1, s1 + s2, p) ## less exact than above?

    res <- 0
    if(prob <= pseudoalpha)
    {
        res <- 1
    }
    else
    {
        h2 <- (p^s2) * ((1 - p)^s1) * choose(s1 + s2, s2)
        if(prob <= (pseudoalpha + h2))
        {
            res <- ((pseudoalpha - prob + h2)/h2)
        }
    }
    return(res)
}
