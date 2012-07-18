\name{npMeanUnpaired}
\alias{npMeanUnpaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the mean of two variables, where the observations do not have
to be paired.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npMeanUnpaired(vector1, vector2, low, up, d, times = 30000, alpha = 0.05, plotif = F, helppoints = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vector1}{
%%     ~~Describe \code{vector1} here~~
}
  \item{vector2}{
%%     ~~Describe \code{vector2} here~~
}
  \item{low}{
%%     ~~Describe \code{low} here~~
}
  \item{up}{
%%     ~~Describe \code{up} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{plotif}{
%%     ~~Describe \code{plotif} here~~
}
  \item{helppoints}{
%%     ~~Describe \code{helppoints} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vector1, vector2, low, up, d, times = 30000, alpha = 0.05, 
    plotif = F, helppoints = 100) 
{
    require(stats)
    var1 <- as.vector(vector1)
    var2 <- as.vector(vector2)
    if (min(var1, var2) < low | max(var1, var2) > up) 
        stop("Your variables include values outside the specified bound")
    d <- d/(up - low)
    x1 <- (var1 - low)/(up - low)
    x2 <- (var2 - low)/(up - low)
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    theta <- function(n1, n2, diff = d, alpha = alpha, plotif = plotif) {
        pvalue.fisher <- function(n1, n2, s1, s2) {
            if (s1 == -1 | s2 > n2) 
                return(0)
            else {
                term <- function(j) {
                  choose(n1, s1 + s2 - j) * choose(n2, j)/choose(n1 + 
                    n2, s1 + s2)
                }
                j <- as.matrix(s2:min(n2, s1 + s2))
                return(sum(apply(j, 1, term)))
            }
        }
        typeII <- function(y1, d, n1, n2, y2 = y1 + d, alpha = alpha, 
            theta = 0.2) {
            pseudoalpha <- theta * alpha
            exmat <- matrix(nrow = n1 + 1, ncol = n2 + 1)
            for (s1 in 0:n1) {
                for (s2 in 0:n2) {
                  t1 <- pvalue.fisher(n1, n2, s1, s2)
                  t2 <- pvalue.fisher(n1, n2, s1 - 1, s2 + 1)
                  if (t1 > pseudoalpha & pseudoalpha > t2) {
                    pr <- (pseudoalpha - t2)/((choose(n1, s1) * 
                      choose(n2, s2))/choose(n1 + n2, s1 + s2))
                  }
                  else {
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
        calctheta <- function(theta, y1 = e1opt, y2 = e2opt, 
            n1, n2, alpha) {
            pseudoalpha <- theta * alpha
            exmat <- matrix(nrow = n1 + 1, ncol = n2 + 1)
            for (s1 in 0:n1) {
                for (s2 in 0:n2) {
                  t1 <- pvalue.fisher(n1, n2, s1, s2)
                  t2 <- pvalue.fisher(n1, n2, s1 - 1, s2 + 1)
                  if (t1 > pseudoalpha & pseudoalpha > t2) {
                    pr <- (pseudoalpha - t2)/((choose(n1, s1) * 
                      choose(n2, s2))/choose(n1 + n2, s1 + s2))
                  }
                  else {
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
        maxexpect <- optimize(typeII, c(0, 1 - diff), tol = 0.001, 
            d = diff, n1 = n1, n2 = n2, alpha = alpha, maximum = T)
        e1opt <- maxexpect$maximum
        e2opt <- e1opt + diff
        if (plotif == T) {
            par(mfrow = c(2, 1))
            y1 <- seq(0, 1 - diff, length = helppoints)
            typeIIerror <- c()
            for (i in 1:helppoints) {
                typeIIerror[i] <- typeII(y1[i], diff, n1, n2, 
                  y2 = y1[i] + diff, alpha = alpha, theta = 0.2)
            }
            plot(typeIIerror ~ y1, type = "l", xlab = "(EY1-low)/(up-low)", 
                ylab = "TypeII Error", main = paste("TypeII Error Plot"))
            abline(v = e1opt, h = maxexpect$objective, lty = 2)
        }
        thetaval <- optimize(calctheta, c(0, 1), tol = 0.001, 
            n1 = n1, n2 = n2, y1 = e1opt, y2 = e2opt, alpha = alpha)
        if (thetaval$objective == 1) 
            stop("TypeII error = 1. Increase difference d")
        print(paste("TypeII error equals", round(thetaval$objective, 
            3), " given choice of d and optimal theta"))
        if (plotif == T) {
            theta1 <- seq(0, 1, length = helppoints)
            typeIIerror1 <- c()
            for (i in 1:helppoints) {
                typeIIerror1[i] <- calctheta(theta1[i], y1 = e1opt, 
                  y2 = e2opt, n1 = n1, n2 = n2, alpha = alpha)
            }
            plot(typeIIerror1 ~ theta1, type = "l", xlab = "theta", 
                ylab = "TypeII Error", main = paste("TypeII Error Plot"), 
                xlim = c(0, 1), ylim = c(0, 0.9))
            abline(v = thetaval$minimum, h = thetaval$objective, 
                lty = 2)
        }
        return(list(typeII = thetaval$objective, theta = thetaval$minimum))
    }
    N1 <- length(var1)
    N2 <- length(var2)
    mean1 <- mean(var1)
    mean2 <- mean(var2)
    thetause <- theta(N1, N2, diff = d, alpha = alpha, plotif = plotif)$theta
    pseudoalpha <- alpha * thetause
    rj <- 0
    for (t in 1:times) {
        q1 <- runif(N1)
        q2 <- runif(N2)
        s1 <- 0
        s2 <- 0
        help1 <- x1 >= q1
        s1 <- sum(help1)
        help2 <- x2 >= q2
        s2 <- sum(help2)
        h1 <- 0
        s3 <- s2 + s1
        s4 <- s3 - N2
        A <- max(0, s4)
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
            }
            else {
                rj <- rj + ((pseudoalpha - h1)/(times * (h2 - 
                  h1)))
            }
        }
    }
    DNAME1 <- deparse(substitute(vector1))
    DNAME2 <- deparse(substitute(vector2))
    if (rj >= thetause) {
        cat("H0: E(", DNAME1, ")>= E(", DNAME2, ") rejected (because P(rj)>=theta)\n")
    }
    else {
        cat("H0: E(", DNAME1, ") >= E(", DNAME2, ") not rejected (because P(rj)<theta)\n")
    }
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
        digits = 3), format(rj, digits = 3)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
