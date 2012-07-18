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
  \item{var1, var2}{
    the numerical vectors which contain the variables to be tested.
}
  \item{low, up}{
    the logical bounds of the variables
}
  \item{d}{
    the hypothesised mean difference, which will be used to minimize the
    type II error.
}
  \item{alpha}{
    the type I error
}
  \item{alternative}{
    a character string describing the alternative hypothesis.
}
  \item{iter}{
Number of iterations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with class '"nphtest"' containing the following components:

  \item{method}{
    a character string indicating the name and type of the test that was
    performed.
  }
  \item{data.name}{
    a character string giving the name(s) of the data.
  }
  \item{alternative}{
    a character string describing the alternative hypothesis.
  }
  \item{estimate}{
    the estimated mean or difference in means depending on
    whether it was a one-sample test or a two-sample test.
  }
  \item{probrej}{
    the probability with which the null hypothesis can be rejected.
  }
  \item{bounds}{
    the lower and upper bounds of the variables.
  }
  \item{null.value}{
    the specified hypothesized value of the correlation between the variables.
  }
  \item{alpha}{
    the type I error
  }
  \item{theta}{
    the parameter that minimizes the type II error.
  }
  \item{pseudoalpha}{
    theta * alpha, this is the alpha that is used during the iterations.
  }
  \item{rejection}{
    Logical. Given alpha, the type I error, and theta, can a rejection be achieved?
  }
  \item{iterations}{
    the number of iterations that were performed.
  }
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
