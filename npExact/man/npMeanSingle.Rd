\name{npMeanSingle}
\alias{npMeanSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the location of a single variable.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npMeanSingle(data_vector, mu, mu1, low, up, times = 20000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_vector}{
The data vector to be tested.
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{mu1}{
%%     ~~Describe \code{mu1} here~~
}
  \item{low}{
%%     ~~Describe \code{low} here~~
}
  \item{up}{
%%     ~~Describe \code{up} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (data_vector, mu, mu1, low, up, times = 20000, alpha = 0.05) 
{
    require(stats)
    varname <- as.vector(data_vector)
    N <- length(varname)
    mean1 <- mean(varname)
    if (min(varname) < low) 
        stop("Your variable includes values below the specified lower bound (low). This contradicts the theoretical requirements of this test.\n")
    if (max(varname) > up) 
        stop("Your variable includes values above the specified upper bound (up). This contradicts the theoretical requirements of this test.\n")
    y <- (varname - low)/(up - low)
    p <- (mu - low)/(up - low)
    p1 <- (mu1 - low)/(up - low)
    y <- as.matrix(y)
    it <- as.numeric(min_value(n = N, p = p, p1 = p1, alpha = alpha))
    theta <- it[1]
    pseudoalpha <- alpha * theta
    Pm <- matrix(p, N, 1)
    Am <- y - Pm
    rj <- 0
    for (t in 1:times) {
        q <- runif(N)
        s1 <- 0
        s2 <- 0
        help2 <- Am > (q * (1 - p))
        s2 <- sum(help2)
        help1 <- y < (q * p)
        s1 <- sum(help1)
        h1 <- 0
        h1 <- sum(dbinom(s2:(s1 + s2), (s1 + s2), p))
        if (h1 <= pseudoalpha) {
            rj <- rj + (1/times)
        }
        else {
            h2 <- (p^s2) * ((1 - p)^s1) * choose(s1 + s2, s2)
            if (h1 <= (pseudoalpha + h2)) {
                rj <- rj + (((pseudoalpha - h1 + h2)/h2)/times)
            }
        }
    }
    DNAME <- deparse(substitute(data_vector))
    if (rj >= theta) {
        cat("H0: E(", DNAME, ")<=", mu, " rejected (because P(rj)>=theta), typeII=", 
            it[2], "\n")
    }
    else {
        cat("H0: E(", DNAME, ")<=", mu, " not rejected (because P(rj)<theta), typeII=", 
            it[2], "\n")
    }
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
        digits = 3)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
