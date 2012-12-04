##########################################################################
# Program Name  stochin_exact                                            #
# Purpose       nonparametric test for stochastic inequality             #
#               based on two independent samples (Schlag, 2008)          #
#               to test of H0: P(Y_2>Y_1)-P(Y2<Y1)<=d    for -1 < d < 1  #
#                                                                        #
# Input variables                                                        #
#        required:       data_vector1                                    #
#                        data_vector2                                    #
#        optional        d      (0) threshold value in null              #
#                        d1     (0.3) test designed to minimize typeII   #
#                               error when H1: P(Y_2>Y_1)-P(Y2<Y1)>=d1   #
#                               so d<d1<1                                #
#                        iterations  (30000) number of Monte Carlo iterations #
#                        alpha  (0.05) significance level of test        #
# Return Values                                                          #
#                        obs1   number of observations in sample1        #
#                        obs2   number of observations in sample2        #
#                        theta  cutoff level theta                       #
#                        alpha  significance level alpha                 #
#                        rj     probability of rejection                 #
#                        sg1    P(Y2>Y1)-P(Y2<Y1)                        #
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

npStochin <- function(x1, x2, d = 0, d1 = 0.3,
                      alternative = "greater",
                      iterations = 5000, alpha = 0.05)
{
  DNAME <- paste(deparse(substitute(x1)), "and",
                 deparse(substitute(x2)))

  x1 <- as.vector(x1)
  x2 <- as.vector(x2)

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  if(alternative == "two.sided")
    stop("Not sensible in this case.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  if (min(d + 1, d1 - d,1 - d1) <= 0)
    stop("we need that -1 < d < d1 < 1")

  if(alternative == "less")
    {
      x1.new <- x2
      x2 <- x1
      x1 <- x1.new
    }

  ## define local variables
  N1 <- length(x1)
  N2 <- length(x2)
  min.length <- min(N1, N2)
  p <- (1 + d)/2

  ## sample estimate
  count.x1 <- 0
  for(i in 1:N1)
    {
      count.x1 <- count.x1 + sum(x1[i] > x2)
    }
  count.x2 <- N1 * N2 - count.x1
  sample.est <- abs(count.x2 - count.x1)/(N1 * N2)
  names(sample.est) <- paste("stochastic inequality:",
                             ifelse(alternative == "greater",
                                    "P(x2 > x1)",
                                    "P(x2 < x1)"))

  ## it <- as.numeric(min_value(n=mi, p=p, p1=(1+d1)/2, alpha=alpha))
  ## if (it[2]>=0.99) stop("increase d1 so that typeII is below 1")
  ## theta <- it[1]
  theta <- .4
  pseudoalpha <- alpha * theta

  rej <- mean(replicate(iterations,
                        sampleBinomTest(x1, x2, min.length,
                                        p, d, pseudoalpha)))

## H0: P(x2 > x1) <= P(x2 < x1), or: P(x2 > x1) - P(x2 < x1) <= d
## H1: P(x2 > x1) > P(x2 < x1), or: P(x2 > x1) - P(x2 < x1) > d
  ## names(d) <- "difference in 'greater-than'-Probabilities"
  names(d) <- "relation P(x2 > x1) - P(x2 < x1)"
  rejection <- ifelse(rej > theta, TRUE, FALSE)

  structure(list(method = "Nonparametric Test for Stochastic Inequality",
                 data.name = DNAME,
                 alternative = alternative,
                 estimate = sample.est,
                 probrej = rej,
                 rejection = rejection,
                 alpha = alpha,
                 theta = theta,
                 iterations = iterations,
                 pseudoalpha = pseudoalpha,
                 ## bounds = bounds,
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
