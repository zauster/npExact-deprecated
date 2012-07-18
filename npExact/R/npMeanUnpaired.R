######################################################################################
##   Program Name      npMeanUnpaired                                          #
##   Purpose           Exact nonparametric test of Schlag (2008)                      #
##                     for comparing expected values given two independent samples    #
##                     for H0: E(Y1) >= E(Y2)  if Y1,Y2 in [low,up]                   #
##                                                                                    #
##   Input variables                                                                  #
##        required:                   y1       independent sample of Y1               #
##                                    y2       independent sample of Y2               #
##                                    low      exogenously known lower                #
##                                             bound for any outcome that can         #
##                                             be generated in either sample          #
##                                    up       exogenously known upper                #
##                                             bound for any outcome that can         #
##                                             be generated in either sample          #
##                                    d        difference (d= EY2 - EY1)              #
##                                             s.t. type II is minimized              #
##                                             when H1: E(Y1)+d <= E(Y2)              #
##                                                                                    #
##        optional (default value):   times    (30000) number of Monte Carlo          #
##                                             iterations                             #
##                                    alpha    (0.05) significance level of test      #
##                                    plotif   (F) if true then plots typeII error    #
##                                             when E(Y2) = E(Y1) + d and minimal     #
##                                             typeII error across different values   #
##                                             of theta                               #
##                                  helppoints (100)  number of points the            #
##                                             graphic uses when plotif=T             #
##                                                                                    #
##   Return Values                                                                    #
##                                    obs1   number of observations in sample1        #
##                                    obs2   number of observations in sample2        #
##                                    low    lower bound for data values              #
##                                    up     upper bound for data values              #
##                                    theta  cutoff level theta                       #
##                                    alpha  significance level alpha                 #
##                                    avg1   average of values in data1               #
##                                    avg2   average of values in data2               #
##                                    P(rj)  probability of rejection in the under    #
##                                           lying randomdomized test                 #
##                                                                                    #
##   Author            Christian Pechhacker                                           #
##   Date              18.02.2012                                                     #
##                                                                                    #
##   using program written by                                                         #
##                                                                                    #
##   Author            Peter Saffert                                                  #
##   Date              19.11.2009                                                     #
##                                                                                    #
## For documentation see                                                              #
##           (Schlag, Karl H. 2008, A New Method for Constructing Exact               #
##           Tests without Making any Assumptions, Department of                      #
##           Economics and Business Working Paper 1109, Universitat                   #
##           Pompeu Fabra)                                                            #
##                                                                                    #
######################################################################################

## data input

## set working directory, eg
## setwd("C:/Users/Schlag/Documents/karl/statistics/R")

## put data as *.dat file (copy from excel into notepad),
## then input
## y1 <- as.matrix(read.table("civil law.dat"))
## y2 <- as.matrix(read.table("common.dat"))

##if your file has comma as decimal separator then do
##
##y1 <- as.matrix(read.table("civil law.dat",dec = ","))
##y2 <- as.matrix(read.table("common.dat",dec = ","))

## Then copy this programm file into R

## Then call function with y1 and y2 for vector1 and vector2.

## Exact test for H0: E(Y1) >= E(Y2) based on Schlag (2008)

npMeanUnpaired <- function(vector1, vector2,
                                  low, up, d, times = 30000,
                                  alpha = 0.05, plotif = F,
                                  helppoints = 100) {
    require(stats)  ##neccessary package for calculation of binomial coefficient

    var1 <- as.vector(vector1)  ##eventuell weitere mögliche Falscheingaben abfangen
    var2 <- as.vector(vector2)
    if (min(var1, var2) < low | max(var1, var2) > up)
        stop("Your variables include values outside the specified bound")

    ## standardize variables
    d <- d/(up - low)
    x1 <- (var1 - low)/(up - low)
    x2 <- (var2 - low)/(up - low)

    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)

###### calculate theta
    theta <- function(n1, n2, diff = d, alpha = alpha,
                      plotif = plotif)
      {

######## define help functions

#### calculates pvalue of Fisher's test
        pvalue.fisher <- function(n1, n2, s1, s2) {
          if (s1 == -1 | s2 > n2)
            return(0)
          else {
            term <- function(j)
              {
                choose(n1, s1 + s2 - j) * choose(n2, j)/choose(n1 +
                                                               n2, s1 + s2)
              }
            j <- as.matrix(s2:min(n2, s1 + s2))
            return(sum(apply(j, 1, term)))
          }
        }

   ## calculates typeII error for given y1, y2
   ## (and of course n1,n2, alpha)

        typeII <- function(y1, d, n1, n2, y2 = y1 + d,
                           alpha = alpha,
                           theta = 0.2)
          {
            pseudoalpha <- theta * alpha
            exmat <- matrix(nrow = n1 + 1,
                            ncol = n2 + 1)
#### needed to calculate typeII error

            for (s1 in 0:n1)
              {
                for (s2 in 0:n2)
                  {
                    t1 <- pvalue.fisher(n1, n2, s1, s2)
                    t2 <- pvalue.fisher(n1, n2, s1 - 1, s2 + 1)
                    if (t1 > pseudoalpha & pseudoalpha > t2)
                      {
                        pr <- (pseudoalpha - t2)/((choose(n1, s1) *
                                                   choose(n2, s2))/choose(n1 + n2, s1 + s2))
                      }
                    else
                      {
                        if (t1 <= pseudoalpha)
                          pr <- 1
                        if (t2 >= pseudoalpha)
                          pr <- 0
                      }
                    exmat[s1 + 1, s2 + 1] <- choose(n1, s1) * y1^s1 *
                      (1 - y1)^(n1 - s1) * choose(n2, s2) * y2^s2 *
                        (1 - y2)^(n2 - s2) * pr
                  }
              }
            type2 <- (1 - sum(exmat))/(1 - theta)
            return(min(type2, 1))
          }

     #### same as function typeII error,  only order of inputs changed, so that function "optimize" can be used
        calctheta <- function(theta, y1 = e1opt, y2 = e2opt,
                              n1, n2, alpha)
          {

            pseudoalpha <- theta * alpha
            exmat <- matrix(nrow = n1 + 1, ncol = n2 + 1)

            for (s1 in 0:n1)
              {
                for (s2 in 0:n2)
                  {
                    t1 <- pvalue.fisher(n1, n2, s1, s2)
                    t2 <- pvalue.fisher(n1, n2, s1 - 1, s2 + 1)
                    if (t1 > pseudoalpha & pseudoalpha > t2)
                      {
                        pr <- (pseudoalpha - t2)/((choose(n1, s1) *
                                                   choose(n2, s2))/choose(n1 + n2, s1 + s2))
                      }
                    else
                      {
                        if (t1 <= pseudoalpha)
                          pr <- 1
                        if (t2 >= pseudoalpha)
                          pr <- 0
                      }
                    exmat[s1 + 1, s2 + 1] <- choose(n1, s1) * y1^s1 *
                      (1 - y1)^(n1 - s1) * choose(n2, s2) * y2^s2 *
                        (1 - y2)^(n2 - s2) * pr
                  }
              }
            type2 <- (1 - sum(exmat))/(1 - theta)
            return(min(type2, 1))
          }


        ## STEP 1) maximize typeII error over y1, y2
        maxexpect <- optimize(typeII, c(0, 1 - diff), tol = 0.001,
                              d = diff, n1 = n1, n2 = n2,
                              alpha = alpha, maximum = T)
        e1opt <- maxexpect$maximum
        e2opt <- e1opt + diff
        ## print(paste('Given theta=0.2 the maximal typeII error
        ## (=', round(maxexpect$objective,3), ') found for EY1
        ## =',round(e1opt,3)))
        if (plotif == T)
          {
            par(mfrow = c(2, 1))
            y1 <- seq(0, 1 - diff, length = helppoints)
            typeIIerror <- c()
            for (i in 1:helppoints)
              {
                typeIIerror[i] <- typeII(y1[i],
                                         diff, n1, n2,
                                         y2 = y1[i] + diff,
                                         alpha = alpha,
                                         theta = 0.2)
              }
            plot(typeIIerror ~ y1, type = "l",
                 xlab = "(EY1-low)/(up-low)",
                 ylab = "TypeII Error",
                 main = paste("TypeII Error Plot"))
            abline(v = e1opt, h = maxexpect$objective, lty = 2)
        }

        ## STEP 2) minimize typeII error over theta
        thetaval <- optimize(calctheta, c(0, 1), tol = 0.001,
                             n1 = n1, n2 = n2,
                             y1 = e1opt, y2 = e2opt,
                             alpha = alpha)

        if (thetaval$objective == 1)
            stop("TypeII error = 1. Increase difference d")

        print(paste("TypeII error equals", round(thetaval$objective,
            3), " given choice of d and optimal theta"))

        if (plotif == T) {
            theta1 <- seq(0, 1, length = helppoints)
            typeIIerror1 <- c()
            for (i in 1:helppoints)
              {
                typeIIerror1[i] <- calctheta(theta1[i],
                y1 = e1opt, y2 = e2opt, n1 = n1, n2 = n2, alpha =
            alpha)
              }

            plot(typeIIerror1 ~ theta1, type = "l", xlab = "theta",
                ylab = "TypeII Error", main = paste("TypeII Error Plot"),
                xlim = c(0, 1), ylim = c(0, 0.9))
            abline(v = thetaval$minimum, h = thetaval$objective,
                lty = 2)
        }

        return(list(typeII = thetaval$objective, theta = thetaval$minimum))
    }

    ## define local variables

    N1 <- length(var1)  ##the stata version uses m instead of N
    N2 <- length(var2)
    mean1 <- mean(var1)
    mean2 <- mean(var2)

    thetause <- theta(N1, N2, diff = d, alpha = alpha, plotif = plotif)$theta
    pseudoalpha <- alpha * thetause

    rj <- 0  ## rejection probability under size pseudoalpha

    for (t in 1:times) {
        ## iteration over t
        q1 <- runif(N1)  ## vector of length N1 with draws from uniform distribution [0,1]
        q2 <- runif(N2)  ## vector of length N2 with draws from uniform distribution [0,1]
        s1 <- 0
        s2 <- 0
        help1 <- x1 >= q1
        s1 <- sum(help1)  ## counts how often values of x1 >= q1
        help2 <- x2 >= q2
        s2 <- sum(help2)  ## counts how often values of x2 >= q2

        ## binomial test {alternative: R command binom.test}
        h1 <- 0
        s3 <- s2 + s1
        s4 <- s3 - N2
        A <- max(0, s4)  ##in the Stata file A=(0  s4) and afterwards only max(A) is used
        if (s1 >= (1 + A)) {
            for (i in A:(s1 - 1)) {
                h1 <- h1 + choose(N1, i) * choose(N2, s3 - i)/choose(N1 +
                  N2, s3)
            }
        }
        if (h1 <= pseudoalpha) {
            h2 <- h1 + choose(N1, s1) * choose(N2, s3 - s1)/choose(N1 +
                N2, s3)
            if (h2 <= pseudoalpha) {
                rj <- rj + (1/times)
            } else {
                rj <- rj + ((pseudoalpha - h1)/(times * (h2 -
                  h1)))
            }
        }
    }  ##end for t in 1:times

    ## rounding is not neccessary -> just limit the output

    DNAME1 <- deparse(substitute(vector1))
    DNAME2 <- deparse(substitute(vector2))

    if (rj >= thetause) {
        cat("H0: E(", DNAME1, ")>= E(", DNAME2, ") rejected (because P(rj)>=theta)\n")
    } else {
        cat("H0: E(", DNAME1, ") >= E(", DNAME2, ") not rejected (because P(rj)<theta)\n")
    }

    ## Output list(obs1=N1, obs2=N2, low=low, up=up,
    ## theta=theta, alpha=alpha, avg1=mean1, avg2=mean2,
    ## Prj=rj) ##nicer output if all variables arranged in a
    ## dataset OR
    names(rj) <- "P(rj)"
    names(N1) <- "obs1"
    names(N2) <- "obs2"
    names(low) <- "low"
    names(up) <- "up"
    names(alpha) <- "alpha"
    names(thetause) <- "theta"
    names(pseudoalpha) <- "pseudoalpha (=alpha*theta)"
    names(mean1) <- "avg1"
    names(mean2) <- "avg2"
    return(c(N1, N2, low, up, format(thetause, digits = 3), format(alpha,
        digits = 3), format(mean1, digits = 3), format(mean2,
        digits = 3), format(rj, digits = 3)))  ##nicer output if all variables arranged in a dataset

}  ##end function npMeanUnpaired

## example 1
## y1 <- rnorm(51,3, 1)
## i <- 1
## while (i <= 51) {
##   y1[i] <- min(max(y1[i],0),10)
##   i <- i+1 }
## y2 <- rnorm(21,6,2)
## i <- 1
## while (i <= 21) {
##   y2[i] <- min(max(y2[i],0),10)
##   i <- i+1 }
## npMeanUnpaired(y1, y2, 0, 10, d=3.3, alpha=0.05)

## npMeanUnpaired(y1, y2, 0, 10, d=3.3, alpha=0.05, plotif=T)
