
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
## npMeanSingle(runif(20), mu = 0.3)

npMeanSingle <- function(x, mu,
                         lower = 0, upper = 1,
                         iterations = 5000, alpha = 0.05,
                         alternative = "greater",
                         epsilon = 1 * 10^(-6),
                         ignoreNA = FALSE)
{
  DNAME <- deparse(substitute(x))

  if(is.data.frame(x))
      {
          if(ncol(x) > 1)
              stop("Please provide 'x' as a vector.")

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
  i <- 1
  rejMatrix <- NULL

  if(alternative == "two.sided")
    {
      ## first the upper side
      optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                               c(0, 1), p = p, N = n,
                               alpha = alpha - epsilon)
      theta <- minTypeIIError(optimaltypeII[[1]],
                              p = p, N = n, alpha = alpha - epsilon)
      pseudoalpha <- (alpha/2) * theta$theta

      while(error > epsilon & i <= 20)
        {
          rejMatrix <- c(rejMatrix,
                              replicate(iterations,
                                        transBinomTest(x, p, xp, n,
                                                       pseudoalpha)))
          rejGreater <- mean(rejMatrix)
          error <- exp(-2 * (iterations * i) * (rejGreater - theta$theta)^2)
          i <- i + 1
        }

      ## secondly the lower side
      x <- 1 - x
      p <- 1 - p
      xp <- x - p
      optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                               c(0, 1), p = p, N = n,
                               alpha = alpha - epsilon)
      theta <- minTypeIIError(optimaltypeII[[1]],
                              p = p, N = n, alpha = alpha - epsilon)
      pseudoalpha <- (alpha/2) * theta$theta

      error <- 1
      i <- 1
      rejMatrix <- NULL
      while(error > epsilon & i <= 20)
        {
          rejMatrix <- c(rejMatrix,
                              replicate(iterations,
                                        transBinomTest(x, p, xp, n,
                                                       pseudoalpha)))
          rejLess <- mean(rejMatrix)
          error <- exp(-2 * (iterations * i) * (rejLess - theta$theta)^2)
          i <- i + 1
        }

      rej <- rejGreater + rejLess
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

      while(error > epsilon & (iterations * i <= 100000))
        {
          rejMatrix <- c(rejMatrix,
                              replicate(iterations,
                                        transBinomTest(x, p, xp, n,
                                                       pseudoalpha)))
          rej <- mean(rejMatrix)
          error <- exp(-2 * (iterations * i) * (rej - theta$theta)^2)
          i <- i + 1
        }
    }

  if(!is.null(iterations) & iterations * (i - 1) < 1000)
    warning("Low number of iterations. Results may be inaccurate.")

  if(iterations * i >= 100000)
    warning("The maximum number of iterations (100,000) was reached. Rejection may be very sensible to the choice of the parameters.")

  method <- "Nonparametric Single Mean Test"
  names(sample.est) <- "mean"
  null.value <- mu
  names(null.value) <- "mean"
  rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
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
                 iterations = iterations * (i - 1),
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
        ## h2 <- (p^ones) * ((1 - p)^zeros) * choose(zeros + ones,
          ## ones)
          h2 <- dbinom(ones, zeros + ones, p)
        if (res.binomtest <= (pseudoalpha + h2))
          {
            res <- ((pseudoalpha - res.binomtest + h2)/h2)
          }
      }
    res
  }
