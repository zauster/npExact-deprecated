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
##                       lower    lower bound for data values              #
##                       upper     upper bound for data values              #
##                       theta  cutoff level theta                       #
##                       alpha  significance level alpha                 #
##                       avg1   average of values in data1               #
##                       avg2   average of values in data2               #
##                       P(rej)  probability of rejection in the under    #
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
                           lower = 0, upper = 1, ## d = 0.3,
                           iterations = 5000,
                           alpha = 0.05,
                           ## plotif = F,
                           ## helppoints = 100,
                           alternative = "greater")
{
    names.x1 <- deparse(substitute(x1))
    names.x2 <- deparse(substitute(x2))

    DNAME <- paste(names.x1, "and", names.x2)

    null.hypothesis <- paste("E(", names.x2, ")",
                            ifelse(alternative == "greater", " <= ",
                                   ifelse(alternative == "less", " >= ",
                                          " = ")),
                            "E(", names.x1, ")", sep = "")

    alt.hypothesis <- paste("E(", names.x2, ")",
                            ifelse(alternative == "greater", " > ",
                                   ifelse(alternative == "less", " < ",
                                          " != ")),
                            "E(", names.x1, ")", sep = "")

  x1 <- as.vector(x1)
  x2 <- as.vector(x2)

  if (min(x1, x2) < lower | max(x1, x2) > upper)
    stop("Some values are out of bounds!")

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

  optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                           c(0, 1), p = 0.5, N = min.length,
                           alpha = alpha)
  theta <- minTypeIIError(optimaltypeII[[1]],
                          p = 0.5, N = min.length, alpha = alpha)
  pseudoalpha <- alpha * theta$theta

  ## thetause <- theta(n1, n2, diff = d,
  ##                   alpha = alpha, plotif = plotif)$theta
  ## thetause <- .4


  if(alternative == "two.sided")
    {
      pseudoalpha <- pseudoalpha/2
      rej.upper <- mean(replicate(iterations,
                            randomTest(x1, x2, n1, n2, pseudoalpha)))
      x1 <- 1 - x1
      x2 <- 1 - x2
      rej.lower <- mean(replicate(iterations,
                            randomTest(x1, x2, n1, n2, pseudoalpha)))
      rej <- rej.upper + rej.lower
    }
  else
    {
      if(alternative == "less")
        {
          x1 <- 1 - x1
          x2 <- 1 - x2
        }
      rej <- mean(replicate(iterations,
                            randomTest(x1, x2, n1, n2, pseudoalpha)))
    }

  method <- "Nonparametric Mean Test for unpaired variables"
  names(sample.est) <- c("mean", "mean")
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
                 iterations = iterations,
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

## ######## define help functions
## #### calculates pvalue of Fisher's test
## pvalue.fisher <- function(n1, n2, s1, s2)
##   {
##     if (s1 == -1 | s2 > n2)
##       return(0)
##     else {
##       term <- function(j)
##         {
##           choose(n1, s1 + s2 - j) * choose(n2, j)/choose(n1 + n2, s1 + s2)
##         }
##       j <- as.matrix(s2:min(n2, s1 + s2))
##       return(sum(apply(j, 1, term)))
##     }
##   }

## ## calculates typeII error for given y1, y2
## ## (and of course n1, n2, alpha)

## typeII <- function(y1, d, n1, n2, y2 = y1 + d,
##                    alpha = alpha,
##                    theta = 0.2)
##   {
##     pseudoalpha <- theta * alpha
##     exmat <- matrix(nrow = n1 + 1,
##                     ncol = n2 + 1)
## #### needed to calculate typeII error

##     for (s1 in 0:n1)
##       {
##         for (s2 in 0:n2)
##           {
##             t1 <- pvalue.fisher(n1, n2, s1, s2)
##             t2 <- pvalue.fisher(n1, n2, s1 - 1, s2 + 1)
##             if (t1 > pseudoalpha & pseudoalpha > t2)
##               {
##                 pr <- (pseudoalpha - t2)/((choose(n1, s1) * choose(n2, s2))/choose(n1 + n2, s1 + s2))
##               }
##             else
##               {
##                 if (t1 <= pseudoalpha)
##                   pr <- 1
##                 if (t2 >= pseudoalpha)
##                   pr <- 0
##               }
##             exmat[s1 + 1, s2 + 1] <- choose(n1, s1) * y1^s1 *
##               (1 - y1)^(n1 - s1) * choose(n2, s2) * y2^s2 *
##                 (1 - y2)^(n2 - s2) * pr
##           }
##       }
##     type2 <- (1 - sum(exmat))/(1 - theta)
##     return(min(type2, 1))
##   }

## ## same as function typeII error, only order of inputs changed, so that
## ## function "optimize" can be used

## calctheta <- function(theta, y1, y2,
##                       n1, n2, alpha)
##   {

##     pseudoalpha <- theta * alpha
##     exmat <- matrix(nrow = n1 + 1, ncol = n2 + 1)

##     for (s1 in 0:n1)
##       {
##         for (s2 in 0:n2)
##           {
##             t1 <- pvalue.fisher(n1, n2, s1, s2)
##             t2 <- pvalue.fisher(n1, n2, s1 - 1, s2 + 1)
##             if (t1 > pseudoalpha & pseudoalpha > t2)
##               {
##                 pr <- (pseudoalpha - t2)/((choose(n1, s1) *
##                                            choose(n2, s2))/choose(n1 + n2, s1 + s2))
##               }
##             else
##               {
##                 if (t1 <= pseudoalpha)
##                   pr <- 1
##                 if (t2 >= pseudoalpha)
##                   pr <- 0
##               }
##             exmat[s1 + 1, s2 + 1] <- choose(n1, s1) * y1^s1 *
##               (1 - y1)^(n1 - s1) * choose(n2, s2) * y2^s2 *
##                 (1 - y2)^(n2 - s2) * pr
##           }
##       }
##     type2 <- (1 - sum(exmat))/(1 - theta)
##     return(min(type2, 1))
##   }

## ###### calculate theta
## theta <- function(n1, n2, diff, alpha = alpha,
##                   plotif = plotif, helppoints = 100)
##   {

##     ## STEP 1) maximize typeII error over y1, y2
##     maxexpect <- optimize(typeII, c(0, 1 - diff), tol = 0.001,
##                           d = diff, n1 = n1, n2 = n2,
##                           alpha = alpha, maximum = T)
##     e1opt <- maxexpect$maximum
##     e2opt <- e1opt + diff
##     ## print(paste('Given theta=0.2 the maximal typeII error
##     ## (=', round(maxexpect$objective,3), ') found for EY1
##     ## =',round(e1opt,3)))
##     ## if (plotif == T)
##     ##   {
##     ##     par(mfrow = c(2, 1))
##     ##     y1 <- seq(0, 1 - diff, length = helppoints)
##     ##     typeIIerror <- c()
##     ##     for (i in 1:helppoints)
##     ##       {
##     ##         typeIIerror[i] <- typeII(y1[i],
##     ##                                  diff, n1, n2,
##     ##                                  y2 = y1[i] + diff,
##     ##                                  alpha = alpha,
##     ##                                  theta = 0.2)
##     ##       }
##     ##     plot(typeIIerror ~ y1, type = "l",
##     ##          xlab = "(EY1-lower)/(upper-lower)",
##     ##          ylab = "TypeII Error",
##     ##          main = paste("TypeII Error Plot"))
##     ##     abline(v = e1opt, h = maxexpect$objective, lty = 2)
##     ## }

##     ## STEP 2) minimize typeII error over theta
##     thetaval <- optimize(calctheta, c(0, 1), tol = 0.001,
##                          n1 = n1, n2 = n2,
##                          y1 = e1opt, y2 = e2opt,
##                          alpha = alpha)

##     if (thetaval$objective == 1)
##       stop("TypeII error = 1. Increase difference d.")

##     ## print(paste("TypeII error equals",
##     ##             round(thetaval$objective, 3),
##     ##             " given choice of d and optimal theta."))

##     ## if (plotif == T)
##     ##   {
##     ##     theta1 <- seq(0, 1, length = helppoints)
##     ##     typeIIerror1 <- c()
##     ##     for (i in 1:helppoints)
##     ##       {
##     ##         typeIIerror1[i] <- calctheta(theta1[i],
##     ##         y1 = e1opt, y2 = e2opt, n1 = n1, n2 = n2, alpha =
##     ##     alpha)
##     ##       }

##     ##     plot(typeIIerror1 ~ theta1, type = "l", xlab = "theta",
##     ##         ylab = "TypeII Error", main = paste("TypeII Error Plot"),
##     ##         xlim = c(0, 1), ylim = c(0, 0.9))
##     ##     abline(v = thetaval$minimum, h = thetaval$objective,
##     ##         lty = 2)
##     ## }

##     return(list(typeII = thetaval$objective,
##                 theta = thetaval$minimum))
##   }
