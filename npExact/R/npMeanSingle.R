#######################################################################################
##   Program Name      nplowerconf_exact                                               #
##   Purpose           exact test for mean of single sample                            #
##                     based on Schlag 2008                                            #
##                     test of EY<=mu based on single sample                           #
##                                                                                     #
##   Input variables                                                                   #
##          required:                   data_vector                                    #
##                                      mu     test value for EY<=mu                   #
##                                      mu1    test optimized to minimize typeII error #
##                                             when EY>=mu1 (so mu1>mu)                #
##                                      low    exogenously known lower                 #
##                                             bound for any outcome that can          #
##                                             be generated in either sample           #
##                                      up     exogenously known upper                 #
##                                             bound for any outcome that can          #
##                                             be generated in either sample           #
##   optional (default value):          times (20000) number of Monte Carlo iterations #
##                                      alpha (.05) significance level of test         #
##   Return Values                                                                     #
##                                      N      number of observations                  #
##                                      low    lower bound for data values             #
##                                      up     upper bound for data values             #
##                                      theta  cutoff level theta                      #
##                                      alpha  significance level alpha                #
##                                      mu     mean to be tested                       #
##                                      avg    average of values in data_vector        #
##                                      rj1    probability of rejection                #
##                                                                                     #
##   Author            Peter Saffert                                                   #
##   Date              19.11.2009                                                      #
##   Version           1.0                                                             #
##                                                                                     #
## For documentation see                                                               #
##           (Schlag, Karl H. 2008, A New Method for Constructing Exact                #
##           Tests without Making any Assumptions, Department of                       #
##           Economics and Business Working Paper 1109, Universitat                    #
##           Pompeu Fabra)                                                             #
##                                                                                     #
#######################################################################################

## Data input

## set working directory, eg
## setwd("C:/Users/Schlag/Documents/karl/statistics/R")

## put data as *.dat file (copy from excel into notepad) into the working directory, then input
## z <- as.matrix(read.table("civil law.dat"))

## if your file has comma as decimal separator then do
## z <- as.matrix(read.table("civil law.dat",dec = ","))

## Then copy this programm file into R

## Then put z below where data_vector is written.


########## Exact test for H0: E(data_vector) <= mu based on Schlag (2008)

npMeanSingle <- function(data_vector, mu, mu1,
    low, up, times = 20000, alpha = 0.05) {
    require(stats)  ##neccessary package for calculation of binomial coefficient

    varname <- as.vector(data_vector)  ##eventuell weitere mögliche Falscheingaben abfangen


    ## define local variables
    N <- length(varname)
    mean1 <- mean(varname)

    if (min(varname) < low)
        stop("Your variable includes values below the specified lower bound (low). This contradicts the theoretical requirements of this test.\n")
    if (max(varname) > up)
        stop("Your variable includes values above the specified upper bound (up). This contradicts the theoretical requirements of this test.\n")

    ## standardize variables
    y <- (varname - low)/(up - low)
    p <- (mu - low)/(up - low)
    p1 <- (mu1 - low)/(up - low)
    y <- as.matrix(y)

    it <- as.numeric(min_value(n = N, p = p, p1 = p1, alpha = alpha))
    theta <- it[1]
    pseudoalpha <- alpha * theta

    Pm <- matrix(p, N, 1)  ## N x 1 matrix with all elements = p
    Am <- y - Pm  ##
    rj <- 0  ## rejection probability under size pseudoalpha

    for (t in 1:times) {
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

        if (h1 <= pseudoalpha) {
            rj <- rj + (1/times)
        } else {
            h2 <- (p^s2) * ((1 - p)^s1) * choose(s1 + s2, s2)
            if (h1 <= (pseudoalpha + h2)) {
                rj <- rj + (((pseudoalpha - h1 + h2)/h2)/times)
            }
        }

    }  ##end for t in 1:times

    ## rouding might not be neccessary in R -> just limit the
    ## output

    DNAME <- deparse(substitute(data_vector))

    if (rj >= theta) {
        cat("H0: E(", DNAME, ")<=", mu, " rejected (because P(rj)>=theta), typeII=",
            it[2], "\n")
    } else {
        cat("H0: E(", DNAME, ")<=", mu, " not rejected (because P(rj)<theta), typeII=",
            it[2], "\n")
    }

    ## list(obs=N, low=low, up=up, theta=theta, alpha=alpha,
    ## mu=mu, avg=mean2, Prj=rj1) OR
    names(rj) <- "P(rj)"
    names(N) <- "obs"
    names(low) <- "low"
    names(up) <- "up"
    names(alpha) <- "alpha"
    names(theta) <- "theta"
    names(pseudoalpha) <- "pseudoalpha (=alpha*theta)"
    names(mean1) <- "avg"
    names(mu) <- "mu"

    return(c(N, low, up, format(theta, digits = 3), format(alpha,
        digits = 3), mu, format(mean1, digits = 3), format(rj,
        digits = 3)))  ##nicer output if all variables arranged in a dataset
}  ##end function npMeanSingle

## example
## z <- rnorm(21,3, 1)
## i <- 1
## while (i <= 21) {
## z[i] <- min(max(z[i],0),10)
## i <- i+1 }
## npMeanSingle(z,mu=1,mu1=3,0,10,alpha=0.05)
