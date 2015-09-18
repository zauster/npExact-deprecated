##########################################################################
# Program Name  stochin_exact                                            #
# Purpose       nonparametric test for stochastic inequality             #
#               based on two independent samples (Schlag, 2008)          #
#               to test of H0: P(Y_2>Y_1)-P(Y2<Y1)<=d    for -1 < d < 1  #
#                                                                        #
#                                                                        #
# Author            Karl Schlag                                          #
# Date              04.03.2012                                           #
#                                                                        #
# based on version of                                                    #
#                                                                        #
# Author            Peter Saffert                                        #
# Date              22.11.2009                                           #
#                                                                        #
#For documentation see                                                   #
#           (Schlag, Karl H. 2008, A New Method for Constructing Exact   #
#           Tests without Making any Assumptions, Department of          #
#           Economics and Business Working Paper 1109, Universitat       #
#           Pompeu Fabra)                                                #
#                                                                        #
##########################################################################

## npStochin
## a function to test for stochastic inequality, i.e.
## H0: P(x2 > x1) <= P(x2 < x1), or: P(x2 > x1) - P(x2 < x1) <= d
## H1: P(x2 > x1) > P(x2 < x1), or: P(x2 > x1) - P(x2 < x1) > d

## x1, x2 ... two data vectors
## alpha ... level of the test
## iterations ... number of iterations
## alternative = "greater" -> x1 >> x2

npStochinUnpaired <- function(x1, x2, d = 0,
                              alternative = "two.sided",
                              iterations = 5000, alpha = 0.05,
                              epsilon = 1 * 10^(-6),
                              ignoreNA = FALSE,
                              max.iterations = 100000)
{
  names.x1 <- deparse(substitute(x1))
  names.x2 <- deparse(substitute(x2))
  DNAME <- paste(names.x1, "and",
                 names.x2)

  x1 <- as.vector(x1)
  x2 <- as.vector(x2)

  if(alternative == "less")
    {
      names.x1.new <- names.x2
      names.x2 <- names.x1
      names.x1 <- names.x1.new
      x1.new <- x2
      x2 <- x1
      x1 <- x1.new
    }

  null.hypothesis <- paste("P(", names.x1, " > ", names.x2, ") - P(",
                           names.x1, " < ", names.x2, ")",
                           ifelse(alternative == "greater", " <= ",
                                  ifelse(alternative == "less", " >= ",
                                         " = ")),
                           d, sep = "")
  alt.hypothesis <- paste("P(", names.x1, " > ", names.x2, ") - P(",
                           names.x1, " < ", names.x2, ")",
                           ifelse(alternative == "greater", " > ",
                                  ifelse(alternative == "less", " < ",
                                         " != ")),
                           d, sep = "")

  if(ignoreNA == TRUE)
    {
      x1 <- x1[!is.na(x1)]
      x2 <- x2[!is.na(x2)]
    }
  else if(any(is.na(c(x1, x2))) == TRUE)
      {
          stop("The data contains NA's!")
      }

  ## if(alternative == "two.sided")
  ##   stop("Not yet implemented. Please test for greater and less at alpha/2.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  ## define local variables
  N1 <- length(x1)
  N2 <- length(x2)
  min.length <- min(N1, N2)
  p <- (1 + d)/2

  ## sample estimate
  count.x1 <- 0
  count.x2 <- 0
  for(i in 1:N1)
    {
      count.x1 <- count.x1 + sum(x1[i] > x2)
      count.x2 <- count.x2 + sum(x1[i] < x2)
    }
  sample.est <- (count.x1 - count.x2)/(N1 * N2)
  names(sample.est) <- paste("stochastic inequality: P(",
                             names.x1,
                             ifelse(alternative == "greater", " > ",
                                    ifelse(alternative == "less", " < ",
                                           " > ")),
                             names.x2, ") - P(",
                             names.x1,
                             ifelse(alternative == "greater", " < ",
                                    ifelse(alternative == "less", " > ",
                                           " < ")),
                             names.x2, ")",
                             sep = "")

  optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                           c(0, 1), p = p, N = min.length,
                           alpha = alpha - epsilon)
  theta <- minTypeIIError(optimaltypeII[[1]],
                          p = p, N = min.length, alpha = alpha - epsilon)
  pseudoalpha <- alpha * theta$theta

  error <- 1
  rejMatrix <- vector(mode = "numeric", length = 0)

  if(alternative == "two.sided")
  {
      while(error > epsilon & length(rejMatrix) <= max.iterations)
      {
          rejMatrix <- c(rejMatrix,
                         replicate(iterations,
                                   sampleBinomTest(x2, x1, min.length,
                                                   p, d, pseudoalpha)))
          rejUpper <- mean(rejMatrix)
          error <- exp(-2 * length(rejMatrix) * (rejUpper - theta$theta)^2)
      }
      error <- 1
      rejMatrix <- vector(mode = "numeric", length = 0)
      ## names.x1.new <- names.x2
      ## names.x2 <- names.x1
      ## names.x1 <- names.x1.new
      x1.new <- x2
      x2 <- x1
      x1 <- x1.new
      while(error > epsilon & length(rejMatrix) <= max.iterations)
      {
          rejMatrix <- c(rejMatrix,
                         replicate(iterations,
                                   sampleBinomTest(x2, x1, min.length,
                                                   p, d, pseudoalpha)))
          rejLess <- mean(rejMatrix)
          error <- exp(-2 * length(rejMatrix) * (rejLess - theta$theta)^2)
      }
      rej <- rejUpper + rejLess

  }
  else
  {
      while(error > epsilon & length(rejMatrix) <= max.iterations)
      {
          rejMatrix <- c(rejMatrix,
                         replicate(iterations,
                                   sampleBinomTest(x2, x1, min.length,
                                                   p, d, pseudoalpha)))
          rej <- mean(rejMatrix)
          error <- exp(-2 * length(rejMatrix) * (rej - theta$theta)^2)
      }
  }

  if(!is.null(iterations) & length(rejMatrix) < 1000)
    warning("Low number of iterations. Results may be inaccurate.")

  if(length(rejMatrix) >= max.iterations)
            warning(paste("The maximum number of iterations (",
                          format(max.iterations, scientific = FALSE),
                          ") was reached. Rejection may be very sensible to the choice of the parameters.", sep = ""))


## H0: P(x1 > x2) <= P(x1 < x2), or: P(x1 > x2) - P(x1 < x2) <= d
## H1: P(x1 > x2) > P(x1 < x2), or: P(x1 > x2) - P(x1 < x2) > d
  ## names(d) <- "difference in 'greater-than'-Probabilities"
  names(d) <- "relation P(x1 > x2) - P(x1 < x2)"
  rejection <- ifelse(rej > theta$theta, TRUE, FALSE)

    ## if rejection in a two.sided setting, we inform the user of the
    ## side of rejection
    if(rejection == TRUE & alternative == "two.sided")
    {
        alt.hypothesis <- paste("P(", names.x1, " > ", names.x2, ") - P(",
                                names.x1, " < ", names.x2, ")",
                                ifelse(rejUpper >= theta$theta, " > ", " < "),
                                d, sep = "")
    }

  structure(list(method = "Nonparametric Test for Stochastic Inequality",
                 data.name = DNAME,
                 alternative = alternative,
                 null.hypothesis = null.hypothesis,
                 alt.hypothesis = alt.hypothesis,
                 estimate = sample.est,
                 probrej = rej,
                 rejection = rejection,
                 alpha = alpha,
                 theta = theta$theta,
                 d.alternative = (optimaltypeII$root*2 - 1),
                 typeIIerror = theta$typeII,
                 iterations = length(rejMatrix),
                 pseudoalpha = pseudoalpha,
                 bounds = NULL,
                 null.value = d),
            class = "nphtest")
}


## sampleBinomTest
## executes a binomial test randomly sampled matched pairs

## x1, x2 ... data vectors
## n ... the minimum length of the data vectors

sampleBinomTest <- function(x1, x2, n, p, d, pseudoalpha)
  {
    c1 <- sample(x1, n)
    c2 <- sample(x2, n)

    s1 <- sum(c1 > c2) #counts how often c1 > c2
    s2 <- sum(c1 < c2) #counts how often c1 < c2

    ## print(s1)
    ## print(s2)
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
                s1 <- s1 + sum(q < (-1)*(d/(1 - d))) ## right?
              }
          }
      }
    prob <- sum(dbinom(s2:(s1+s2), (s1+s2), p)) ## or
    ## prob <- 1 - pbinom(s2 - 1, s1 + s2, p) ## less exact than above?

    res <- 0
    if(prob <= pseudoalpha)
      {
        res <- 1
      }
    else
      {
        h2 <- (p^s2)*((1-p)^s1)*choose(s1+s2,s2)
        if(prob <= (pseudoalpha + h2))
          {
            res <- ((pseudoalpha - prob + h2)/h2)
          }
      }
    return(res)
  }
