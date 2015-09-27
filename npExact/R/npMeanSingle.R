## npMeanSingle
## a function to test the mean of a variable

## x ... the data vector
## mu ... the value to be tested, the supposed true value of the mean
## lower, upper ... logical bounds of x, thus x in [lower, upper]
## iterations ... number of iterations
## alpha ... the level of the test
## alternative ... a character string specifying the alternative
## hypothesis, must be one of 'two.sided', 'greater' or 'less'.

## example
## x <- runif(20)
## npMeanSingle(x, mu = 0.3)

## npMeanSingle(x, mu = .9, alternative = "two.sided")
## npMeanSingle(x, mu = .9, alternative = "greater")
## npMeanSingle(x, mu = .9, alternative = "less")

npMeanSingle <- function(x, mu,
                         lower = 0, upper = 1,
                         alternative = "two.sided",
                         iterations = 5000, alpha = 0.05,
                         epsilon = 1 * 10^(-6),
                         ignoreNA = FALSE,
                         max.iterations = 100000)
{
    DNAME <- deparse(substitute(x))

    if(is.data.frame(x))
    {
        if(ncol(x) > 1)
            stop("Please provide 'x' as a vector or a single column data.frame.")

        x <- x[,1]
    }

    x <- as.vector(x)
    sample.est <- mean(x)

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
        ## first the upper side
        ##
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = p, N = n,
                                 alpha = alpha / 2 - epsilon)
        thetaUpper <- minTypeIIError(optimaltypeII[[1]],
                                     p = p, N = n, alpha = alpha / 2 - epsilon)
        pseudoalpha <- (alpha/2) * thetaUpper$theta

        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     transBinomTest(x, p, xp, n,
                                                    pseudoalpha)))
            rejUpper <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rejUpper - thetaUpper$theta)^2)
        }
        rejectionUpper <- ifelse(rejUpper > thetaUpper$theta, TRUE, FALSE)


        ##
        ## secondly the lower side
        ##
        x <- 1 - x
        p <- 1 - p
        xp <- x - p
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = p, N = n,
                                 alpha = alpha / 2 - epsilon)
        thetaLess <- minTypeIIError(optimaltypeII[[1]],
                                     p = p, N = n,
                                     alpha = alpha / 2 - epsilon)
        pseudoalpha <- (alpha/2) * thetaLess$theta

        error <- 1
        rejMatrix <- vector(mode = "numeric", length = 0)
        while(error > epsilon & length(rejMatrix) <= max.iterations)
        {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     transBinomTest(x, p, xp, n,
                                                    pseudoalpha)))
            rejLess <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rejLess - thetaLess$theta)^2)
        }
        rejectionLess <- ifelse(rejLess > thetaLess$theta, TRUE, FALSE)


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

        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = p, N = n,
                                 alpha = alpha - epsilon)
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
    }

    if(!is.null(iterations) & length(rejMatrix) < 1000)
        warning("Low number of iterations. Results may be inaccurate.")

    if(length(rejMatrix) >= max.iterations)
        warning(paste("The maximum number of iterations (",
                      format(max.iterations, scientific = FALSE),
                      ") was reached. Rejection may be very sensible to the choice of the parameters.", sep = ""))

    method <- "Nonparametric Single Mean Test"
    names(sample.est) <- "mean"
    null.value <- mu
    names(null.value) <- "mean"

    ## if rejection in a two.sided setting, we inform the user of the
    ## side of rejection
    if(alternative == "two.sided")
    {
        if(rejection == TRUE)
        {
alt.hypothesis <- paste("E(", DNAME, ")",
                                ifelse(rejectionUpper == TRUE, " > ", " < "),
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
                   iterations = length(rejMatrix),
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
