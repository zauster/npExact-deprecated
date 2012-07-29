######################################################################
##  Program Name      nplowerconf_exact                              #
##  Purpose           exact test for mean of single sample           #
##                    based on Schlag 2008                           #
##                    test of EY<=mu based on single sample          #
##                                                                   #
##  Input variables                                                  #
##         required:  x                                              #
##                    mu     test value for EY<=mu                   #
##                    mu1    test optimized to minimize typeII error #
##                           when EY>=mu1 (so mu1>mu)                #
##                    low    exogenously known lower                 #
##                           bound for any outcome that can          #
##                           be generated in either sample           #
##                    up     exogenously known upper                 #
##                           bound for any outcome that can          #
##                           be generated in either sample           #
##  optional (default value):
##                    times (20000) number of Monte Carlo iterations #
##                    alpha (.05) significance level of test         #
##
##  Return Values                                                    #
##                    N      number of observations                  #
##                    low    lower bound for data values             #
##                    up     upper bound for data values             #
##                    theta  cutoff level theta                      #
##                    alpha  significance level alpha                #
##                    mu     mean to be tested                       #
##                    avg    average of values in x                  #
##                    rj1    probability of rejection                #
##                                                                   #
##  Author            Peter Saffert                                  #
##  Date              19.11.2009                                     #
##  Version           1.0                                            #
##                                                                   #
## For documentation see                                             #
##          (Schlag, Karl H. 2008, A New Method for Constructing Exact
##          Tests without Making any Assumptions, Department of      #
##          Economics and Business Working Paper 1109, Universitat   #
##          Pompeu Fabra)                                            #
##                                                                   #
######################################################################

########## Exact test for H0: E(x) <= mu based on Schlag (2008)

npMeanSingle <- function(x, mu, mu1,
                         low = 0, up = 1,
                         iterations = 20000, alpha = 0.05,
                         alternative = "greater")
{
  require(stats)
  DNAME <- deparse(substitute(x))
  sample.est <- mean(x)

  x <- as.vector(x)

  if (min(x) < low | max(x) > up)
    stop("Some values are out of bounds!")

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  if(alternative == "two.sided")
    stop("Currently not supported. Please test for greater and less at alpha/2.")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")

  N <- length(x)

  ## standardize variables
  y <- (x - low)/(up - low)
  p <- (mu - low)/(up - low)
  p1 <- (mu1 - low)/(up - low)
  y <- as.matrix(y)

  it <- as.numeric(min_value(n = N, p = p, p1 = p1, alpha = alpha))
  theta <- it[1]
  pseudoalpha <- alpha * theta

  Pm <- matrix(p, N, 1)  ## N x 1 matrix with all elements = p
  Am <- y - Pm  ##
  rj <- 0  ## rejection probability under size pseudoalpha

  for (t in 1:iterations)
    {
      ## iteration over t
      q <- runif(N)  ## vector of length N1 with draws from uniform distribution [0,1]
      s1 <- 0
      s2 <- 0
      help2 <- Am > (q * (1 - p))
      s2 <- sum(help2)  ## counts how often values of Am > q*(1-p)
      help1 <- y < (q * p)
      s1 <- sum(help1)  ## counts how often values of y < q*p

      ## binomial test {alternative: R command binom.test}
      h1 <- 0
      ## for(k in s2:(s1+s2)){ h1 <- h1 + (p^k) *
      ## ((1-p)^(s1+s2-k)) * choose(s1+s2,k) ##alternative use
      ## dbinom or pbinom command which are computational more
      ## efficient than a for loop } ##end for k in s2:s1+s2
      h1 <- sum(dbinom(s2:(s1 + s2), (s1 + s2), p))

      if (h1 <= pseudoalpha)
        {
          rj <- rj + (1/iterations)
        }
      else
        {
          h2 <- (p^s2) * ((1 - p)^s1) * choose(s1 + s2, s2)
          if (h1 <= (pseudoalpha + h2))
            {
              rj <- rj + (((pseudoalpha - h1 + h2)/h2)/iterations)
            }
        }

    }  ##end for t in 1:iterations

  ## rouding might not be neccessary in R -> just limit the
  ## output


  method <- "Nonparametric Single Mean Test"
  names(sample.est) <- "mean"
  null.value <- mu
  names(null.value) <- "mean"
  rejection <- ifelse(rj >= theta, TRUE, FALSE)
  bounds <- paste("[", low, ", ", up, "]", sep = "")

  alternative <- "greater"

  structure(list(method = method,
                   data.name = DNAME,
                   alternative = alternative,
                   estimate = sample.est,
                   probrej = rj,
                   rejection = rejection,
                   alpha = alpha,
                   theta = theta,
                   iterations = iterations,
                   pseudoalpha = pseudoalpha,
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
}  ##end function npMeanSingle

## example
## npMeanSingle(runif(20), mu = 0.5, mu1 = 0.6)
