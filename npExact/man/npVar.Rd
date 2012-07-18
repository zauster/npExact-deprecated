\name{npVar}
\alias{npVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the variance of a variable.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npVar(test1, low, up, w, w1, type = "upper", alpha = 0.025, T = 50000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test1}{
%%     ~~Describe \code{test1} here~~
}
  \item{low}{
%%     ~~Describe \code{low} here~~
}
  \item{up}{
%%     ~~Describe \code{up} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{w1}{
%%     ~~Describe \code{w1} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
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
function (test1, low, up, w, w1, type = "upper", alpha = 0.025, 
    T = 50000) 
{
    if (min(test1 - low, up - test1) < 0) 
        stop("values have to lie in [low,up]")
    if (w > (up - low)/4) 
        stop("w too large")
    if (w1 >= w & type == "upper") 
        stop("need w1 < w when type is upper")
    if (w1 <= w & type == "lower") 
        stop("need w1 > w when type is lower")
    m1 <- length(test1)
    m <- floor(m1/2)
    x <- (test1 - low)/(up - low)
    w0 <- w/(up - low)^2
    p <- 2 * w0
    p1 <- 2 * w1/(up - low)^2
    mean1 <- sum(test1/m1)
    var1 <- sum((test1 - mean1)^2/(m1 - 1))
    if (type == "upper") {
        it <- as.numeric(min_value(n = m, p = 1 - p, p1 = 1 - 
            p1, alpha = alpha))
        if (it[2] >= 0.99) 
            stop("decrease w1 so that typeII is below 1")
    }
    else {
        it <- as.numeric(min_value(n = m, p = p, p1 = p1, alpha = alpha))
        if (it[2] >= 0.99) 
            stop("increase w1 so that typeII is below 1")
    }
    theta <- it[1]
    pseudoalpha <- alpha * theta
    cat(paste("Sample N = ", m1, ",   avg = ", round(mean1, 3), 
        ",   var = ", round(var1, 3), ",  theta =  ", round(theta, 
            3), ",  alpha =  ", alpha, sep = ""))
    cat("\n")
    rj <- 0
    for (t in 1:T) {
        c1 <- sample(x)
        y <- (c1[c(1:m) * 2] - c1[(c(1:m) * 2 - 1)])^2
        s1 <- 0
        s2 <- 0
        q <- runif(m)
        s2 <- length(which(y - p > (q * (1 - p))))
        s1 <- length(which(y < q * p))
        h1 <- 0
        k <- switch(type, upper = 0:s2, lower = s2:(s1 + s2))
        h1 <- sum((p^k) * ((1 - p)^(s1 + s2 - k)) * choose((s1 + 
            s2), k))
        if (h1 <= pseudoalpha) {
            rj <- rj + 1/T
        }
        else {
            h <- (p^s2) * ((1 - p)^s1) * choose(s1 + s2, s2)
            if (h1 <= pseudoalpha + h) {
                rj = rj + ((pseudoalpha - h1 + h)/h)/T
            }
        }
    }
    if (rj >= theta) {
        cat(paste("H0: VarY ", switch(type, upper = ">=", lower = "<="), 
            w, " REJECTED  (rj=  ", round(rj, 4), "), typeII =", 
            round(it[2], 4), sep = ""))
    }
    else {
        cat(paste("H0: VarY ", switch(type, upper = ">=", lower = "<="), 
            w, " NOT rejected  (rj=  ", round(rj, 4), "), typeII =", 
            round(it[2], 4), sep = ""))
    }
    cat("\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
