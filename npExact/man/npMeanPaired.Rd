\name{npMeanPaired}
\alias{npMeanPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the mean difference of two variables, if the observations are paired.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npMeanPaired(var1, var2, low = 0, up = 1, d, alpha = 0.05, alternative = "greater", iter = 30000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var1}{
%%     ~~Describe \code{var1} here~~
}
  \item{var2}{
%%     ~~Describe \code{var2} here~~
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
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
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
function (var1, var2, low = 0, up = 1, d, alpha = 0.05, alternative = "greater", 
    iter = 30000) 
{
    DNAME <- paste(deparse(substitute(var1)), "and", deparse(substitute(var2)))
    var1 <- as.vector(var1)
    var2 <- as.vector(var2)
    if (min(var1) < low) 
        stop("Your variable includes values below the specified lower bound (low). This contradicts the theoretical requirements of this test.\n")
    if (max(var2) > up) 
        stop("Your variable includes values above the specified upper bound (up). This contradicts the theoretical requirements of this test.\n")
    n <- length(var1)
    mean1 <- mean(var1)
    mean2 <- mean(var2)
    var1 <- (var1 - low)/(up - low)
    var2 <- (var2 - low)/(up - low)
    d <- d/(up - low)
    theta <- 0.4
    pseudoalpha <- alpha * theta
    res <- replicate(iter, doMatching(randTransformation(var1, 
        n), randTransformation(var2, n), pseudoalpha))
    rej <- mean(res)
    method <- "Nonparametric Paired Mean Test"
    sample.est <- c(mean(var1), mean(var2))
    names(sample.est) <- c("mean", "mean")
    null.value <- 0
    names(null.value) <- "mean difference"
    rejection <- ifelse(rej > theta, TRUE, FALSE)
    bounds <- paste("[", low, ", ", up, "]", sep = "")
    structure(list(method = method, data.name = DNAME, alternative = alternative, 
        estimate = sample.est, probrej = rej, rejection = rejection, 
        alpha = alpha, theta = theta, iterations = iter, pseudoalpha = pseudoalpha, 
        bounds = bounds, null.value = null.value), class = "nphtest")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
