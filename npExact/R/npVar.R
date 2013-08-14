###################################################
###         Testing Variance of Singe Sample     ##
###         Upper and Lower Confidence Bounds    ##
###                   R Code                     ##
###################################################

##-------------------------------------------------
## Reference: K. H. Schlag (2008), A New Method for Constructing Exact
## Tests without Making any Assumptions, Universitat Pompeu Fabra working
## paper Working 1109
##-------------------------------------------------

## This is a program that implements a nonrandomized one-sided test for
## the variance of a single sample of independent observations.

## Null hypothesis: H0: VarY >= w, if type="greater",
##                      finetuning test for H1: VarY <= w1
## Null hypothesis: H0: VarY <= w, if type="lower",
##                      finetuning test for H1: VarY >= w1

## This test can be used to derive (1 - alpha) UPPER CONFIDENCE BOUNDS,
## type="greater": Find largest w such that not rejects the null

## This test can be used to derive (1 - alpha) LOWER CONFIDENCE BOUNDS,
## type="lower": Find smallest w such that not rejects the null.

## The test transforms the data into a new sample that has mean equal to
## variance of orginal process and then applies the A test.

## example
## npVar(runif(50), low = 0, up = 1, var = 0.1, alt.var = 0.15)

npVarianceSingle <- function(x, lower = 0, upper = 1, variance,
                  alternative = "greater",
                  alpha = 0.05, iterations = 5000,
                  epsilon = 1 * 10^(-6),
                  ignoreNA = FALSE)
{
  DNAME <- deparse(substitute(x))
  sample.est <- var(x)
  x <- as.vector(x)

  null.hypothesis <- paste("Var(", DNAME, ") ",
                           ifelse(alternative == "greater", "<= ",
                                  ifelse(alternative == "less", ">= ",
                                         "= ")),
                           variance, sep = "")
  alt.hypothesis <- paste("Var(", DNAME, ") ",
                          ifelse(alternative == "greater", "> ",
                                 ifelse(alternative == "less", "< ",
                                        "!= ")),
                          variance, sep = "")

  if(ignoreNA == TRUE)
    {
      x <- x[!is.na(x)]
    }

  if (min(x) < lower | max(x) > upper)
    stop("Some values are out of bounds (or NA)!")

  if (variance > 0.25*(upper - lower)^2)
    stop("Hypothesized 'variance' is too large.")

  ## if (alt.variance >= variance & alternative == "less")
  ##   stop("A 'alt.variance' < 'variance' is needed when alternative is 'less'.")

  ## if (alt.variance <= variance & alternative == "greater")
  ##   stop("A 'alt.variance' > 'variance' is needed when alternative is 'greater'.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  ## Computation of sample mean and variance for output
  m <- floor(length(x) / 2)
  x <- (x - lower)/(upper - lower)  ## Normalization so that x in [0,1]
  p <- 2 * variance / (upper - lower)^2  ## normalized threshold

  error <- i <- 1
  rejMatrix <- NULL

  ## deterministic test
  if(alternative == "two.sided")
    {
        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = p, N = m,
                                 alpha = alpha)
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = p, N = m, alpha = alpha)
        pseudoalpha <- alpha * theta$theta / 2
        while(error > epsilon & i <= 20)
          {
              rejMatrix <- cbind(rejMatrix,
                                 replicate(iterations,
                                           sampleBinomTestnpVar(x, m, p,
                                                                alternative = "greater",
                                                                alpha = pseudoalpha)))
              rejUpper <- mean(rejMatrix)
              error <- exp(-2 * (iterations * i) * (rejUpper - theta$theta)^2)
              i <- i + 1
          }

        optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                 c(0, 1), p = 1 - p, N = m,
                                 alpha = alpha)
        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = 1 - p, N = m, alpha = alpha)
        pseudoalpha <- alpha * theta$theta / 2

        error <- i <- 1
        rejMatrix <- NULL

        while(error > epsilon & i <= 20)
          {
              rejMatrix <- cbind(rejMatrix,
                                 replicate(iterations,
                                           sampleBinomTestnpVar(x, m, p,
                                                                alternative = "less",
                                                                alpha = pseudoalpha)))
              rejLess <- mean(rejMatrix)
              error <- exp(-2 * (iterations * i) * (rejLess - theta$theta)^2)
              i <- i + 1
          }
        rej <- rejUpper + rejLess
    }
      else
        {
            optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                     c(0, 1),
                                     p = ifelse(alternative == "greater",
                                         p, 1 - p),
                                     N = m,
                                     alpha = alpha)
            theta <- minTypeIIError(optimaltypeII[[1]],
                                    p = ifelse(alternative == "greater",
                                        p, 1 - p),
                                    N = m, alpha = alpha)
            pseudoalpha <- alpha * theta$theta

            while(error > epsilon & i <= 20)
              {
                  rejMatrix <- cbind(rejMatrix,
                                     replicate(iterations,
                                               sampleBinomTestnpVar(x, m, p,
                                                                    alternative = alternative,
                                                                    alpha = pseudoalpha)))
                  rej <- mean(rejMatrix)
                  error <- exp(-2 * (iterations * i) * (rej - theta$theta)^2)
                  i <- i + 1
              }
        }

  if(i == 21)
    warning("The maximum number of iterations (100,000) was reached. Rejection may be very sensible to the choice of the parameters.")

  names(sample.est) <- "variance"
  null.value <- variance
  names(null.value) <- "variance"
  bounds <- paste("[", lower, ", ", upper, "]", sep = "")
  rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
  ## alternative <- ifelse(alternative == "greater", "less",
  ## "greater")

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
                 iterations = iterations * (i - 1),
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

    res <- 0
    if(prob <= alpha) ## reject with probability 1
      {
        res <- 1
      }
    else
      {
        ## h <- (p^zeros)*((1 - p)^ones) * choose(ones + zeros, zeros)
        h <- dbinom(zeros, zeros + ones, p) ## more efficient
        if (prob <= alpha + h) ##(reject with positive probability)
          {
            res <- ((alpha - prob + h) / h)
          }
      }
    return(res)
  }
