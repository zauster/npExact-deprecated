##########################################################################
##   Program Name      npMeanUnpaired                                    #
##   Purpose           Exact nonparametric test of Schlag (2008)         #
##                     for comparing expected values given               #
##                     two independent samples                           #
##                     for H0: E(Y1) >= E(Y2)  if Y1,Y2 in
##                     [lower,upper]
##                                                                       #
##   Input variables                                                     #
##        required:      y1       independent sample of Y1               #
##                       y2       independent sample of Y2               #
##                       lower      exogenously known lower              #
##                                bound for any outcome that can         #
##                                be generated in either sample          #
##                       upper       exogenously known upper             #
##                                bound for any outcome that can         #
##                                be generated in either sample          #
##                       d        difference (d= EY2 - EY1)              #
##                                s.t. type II is minimized              #
##                                when H1: E(Y1)+d <= E(Y2)              #
##                                                                       #
##        optional (default value):
##                  iterations    (30000) number of Monte Carlo          #
##                                iterations                             #
##                  alpha    (0.05) significance level of test           #
##                  plotif   (F) if true then plots typeII error         #
##                           when E(Y2) = E(Y1) + d and minimal          #
##                           typeII error across different values        #
##                           of theta                                    #
##                  helppoints (100)  number of points the               #
##                           graphic uses when plotif=T                  #
##                                                                       #
##   Return Values                                                       #
##                       obs1   number of observations in sample1        #
##                       obs2   number of observations in sample2        #
##                       lower    lower bound for data values            #
##                       upper     upper bound for data values           #
##                       theta  cutoff level theta                       #
##                       alpha  significance level alpha                 #
##                       avg1   average of values in data1               #
##                       avg2   average of values in data2               #
##                       P(rej)  probability of rejection in the under   #
##                              lying randomdomized test                 #
##                                                                       #
##   Author            Christian Pechhacker                              #
##   Date              18.02.2012                                        #
##                                                                       #
##   using program written by                                            #
##                                                                       #
##   Author            Peter Saffert                                     #
##   Date              19.11.2009                                        #
##                                                                       #
## For documentation see                                                 #
##           (Schlag, Karl H. 2008, A New Method for Constructing Exact  #
##           Tests without Making any Assumptions, Department of         #
##           Economics and Business Working Paper 1109, Universitat      #
##           Pompeu Fabra)                                               #
##                                                                       #
##########################################################################

## Exact test for H0: E(Y1) >= E(Y2) based on Schlag (2008)

## Test for H0: E(x1) >= E(x2), or: E(x1) - E(x2) >= 0 against
## H1: E(x1) < E(x2), or E(x1) - E(x2) < 0

npMeanUnpaired <- function(x1, x2,
                           lower = 0, upper = 1,
                           iterations = 5000,
                           alpha = 0.05,
                           alternative = "greater",
                           epsilon = 1 * 10^(-6),
                           ignoreNA = FALSE)
{
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

  x1 <- as.vector(x1)
  x2 <- as.vector(x2)

  if(ignoreNA == TRUE)
    {
      x1 <- x1[!is.na(x1)]
      x2 <- x2[!is.na(x2)]
    }

  if (min(x1, x2) < lower | max(x1, x2) > upper)
    stop("Some values are out of bounds (or NA)!")

  if(iterations < 500)
    warning("Low number of iterations. Results may be inaccurate.")

  if(alternative != "two.sided" & alternative != "greater" & alternative != "less")
    stop("Please specify which alternative hypothesis you want to test for: 'greater', 'less' or 'two.sided'")

  if(alpha >= 1 | alpha <= 0)
    stop("Please supply a sensible value for alpha.")


  sample.est <- c(mean(x1), mean(x2))

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

  optimaltypeII <- optimize(npMeanUnpairedminTypeIIErrorWrapper,
                            c(0, 1), n1 = n1, n2 = n2, alpha = alpha,
                            tol = .Machine$double.eps^0.25)
  theta <- optimizeTheta(n1, n2, optimaltypeII$minimum, alpha)
  pseudoalpha <- alpha * theta$theta

  error <- i <- 1
  rejMatrix <- NULL

  if(alternative == "two.sided")
    {
        pseudoalpha <- pseudoalpha/2
        while(error > epsilon & i <= 20)
          {
              rejMatrix <- cbind(rejMatrix,
                                 replicate(iterations,
                                           randomTest(x1, x2, n1, n2, pseudoalpha)))
              rejUpper <- mean(rejMatrix)
              error <- exp(-2 * (iterations * i) * (rejUpper - theta$theta)^2)
              i <- i + 1
          }
        x1 <- 1 - x1
        x2 <- 1 - x2

        error <- i <- 1
        rejMatrix <- NULL
        while(error > epsilon & i <= 20)
          {
              rejMatrix <- cbind(rejMatrix,
                                 replicate(iterations,
                                           randomTest(x1, x2, n1, n2, pseudoalpha)))
              rejLess <- mean(rejMatrix)
              error <- exp(-2 * (iterations * i) * (rejLess - theta$theta)^2)
              i <- i + 1
          }

        rej <- rejUpper + rejLess
    }
      else
        {
            if(alternative == "greater")
              {
                  x1 <- 1 - x1
                  x2 <- 1 - x2
              }
            while(error > epsilon & i <= 20)
              {
                  rejMatrix <- cbind(rejMatrix,
                                     replicate(iterations,
                                               randomTest(x1, x2,
                                                          n1, n2,
                                                          pseudoalpha)))
                  rej <- mean(rejMatrix)
                  error <- exp(-2 * (iterations * i) * (rej - theta$theta)^2)
                  i <- i + 1
              }
        }

  if(i == 21)
    warning("The maximum number of iterations (100,000) was reached. Rejection may be very sensible to the choice of the parameters.")

  method <- "Nonparametric Mean Test for unpaired variables"
  names(sample.est) <- c(paste("mean(", names.x1, ")", sep = ""),
                         paste("mean(", names.x2, ")", sep = ""))
  null.value <- 0
  names(null.value) <- "E[x2] - E[x1]" ## "mean difference"
  rejection <- ifelse(rej >= theta$theta, TRUE, FALSE)
  bounds <- paste("[", lower, ", ", upper, "]", sep = "")

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
                 d.alternative = optimaltypeII$minimum,
                 typeIIerror = theta$typeII,
                 iterations = iterations * (i - 1),
                 pseudoalpha = pseudoalpha,
                 bounds = bounds,
                 null.value = null.value),
            class = "nphtest")
} ## end of npMeanUnpaired

## example
## npMeanUnpaired(runif(20), runif(20))

randomTest <- function(x1, x2, n1, n2, alpha)
  {
    s1 <- sum(x1 >= runif(n1))
    s2 <- sum(x2 >= runif(n2))
    s3 <- s2 + s1
    k <- max(0, s3 - n2):(s1 - 1)
    prob <- 0 ## was 1. Why?
    if (s1 >= (1 + k[1]))
      {
        prob <- sum(choose(n1,
                           k) * choose(n2,
                                       s3 - k)/choose(n1 + n2,
                                                      s3))
        ## h.alt <- phyper(s1 - 1, n1, n2, A)
      }

    res <- 0
    if (prob <= alpha)
      {
        ## h2 <- prob + choose(n1, s1) * choose(n2, s3 - s1)/choose(n1 +
        ## n2, s3)
        h2 <- prob + dhyper(s1, n1, n2, s3)
        if (h2 <= alpha)
          {
            res <- 1
          }
        else
          {
            res <- ((alpha - prob)/((h2 - prob)))
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

## calculates typeII error for given y1, y2 (and of course n1,n2, alpha)
maxTypeII <- function(y1, d, n1, n2, y2 = y1 + d,
                      alpha = alpha, theta = 0.2)
{
  pseudoalpha <- theta * alpha
  exmat <- matrix(nrow = n1 + 1, ncol = n2 + 1)
  ## needed to calculate typeII error
  res <- 0
  for(s1 in 0:n1)
    {
  ##     for(s2 in 0:n2)
  ##       {
  ##         t1 <- pvalueFisher(n1, n2, s1, s2)
  ##         t2 <- pvalueFisher(n1, n2, s1 - 1, s2 + 1)

  ##         if( t2 >= pseudoalpha)
  ##           {
  ##             ## in this case, pr = zero
  ##             ## so we can skip the calculation of the term
  ##             exmat[s1 + 1, s2 + 1] <- 0
  ##           }
  ##         else
  ##           {
  ##             if( t1  > pseudoalpha & pseudoalpha > t2)
  ##               {
  ##                 ## pr <- (pseudoalpha - t2) / ((choose(n1, s1) *
  ##                 ## choose(n2, s2))/ choose(n1 + n2, s1 + s2))
  ##                 pr <- (pseudoalpha - t2) / dhyper(s1, n1, n2, s1 + s2)
  ##               }
  ##             else
  ##               {
  ##                 if(t1 <= pseudoalpha) pr <- 1
  ##               }
  ##             exmat[s1 + 1,
  ##                   s2 + 1] <- dbinom(s1, n1, y1) * dbinom(s2, n2, y2) * pr
  ##           }
  ##       }
  ## }
      ## now instead of the second for clause -> vectorized if-clauses!
      t1 <- pvalueFisher(n1, n2, s1, 0:n2)
      t2 <- pvalueFisher(n1, n2, s1 - 1, 1:(n2 + 1))
      res <- res + sum(ifelse(t2 >= pseudoalpha, 0,
                    ifelse(t1 > pseudoalpha & pseudoalpha > t2,
                           dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * (pseudoalpha - t2) / dhyper(s1, n1, n2, s1 + 0:n2),
                           dbinom(s1, n1, y1) * dbinom(0:n2, n2, y2) * 1)))
    }
  type2 <- (1 - res) / (1 - theta)
  ## type2 <- (1 - sum(exmat)) / (1 - theta)
  return(min(type2, 1))
}

## mean(system.time(maxTypeII))

## same as function typeII error,  only order of inputs changed,
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
  maxexpect <- optimize(maxTypeII, c(0, 1 - diff),
                        ## tol = .Machine$double.eps^0.5,
                        d = diff, n1 = n1, n2 = n2,
                        alpha = alpha, maximum = T)
  e1opt <- maxexpect$maximum
  e2opt <- e1opt + diff

  ## STEP 2)  minimize typeII error over theta
  thetaval <- optimize(minTypeII, c(0,1),
                       ## tol = .Machine$double.eps^0.5,
                       n1 = n1, n2 = n2,
                       y1 = e1opt, y2 = e2opt, alpha = alpha)

  ## if(thetaval$objective == 1)
  ##   stop("TypeII error = 1. Increase difference d")

  return(list(typeII = thetaval$objective,
              theta = thetaval$minimum ))
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
