\name{npVarianceSingle}
\alias{npVarianceSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test for the variance of a bounded random variable based on a single sample of iid observations.
}
\description{
This test requires that the user knows upper and lower bounds before gathering the data
such that the properties of the data generating process imply that all observations will be within these bounds.
The data input consists of a sequence of observations, each being an independent realization of the random variable.
No further distributional assumptions are made.

This is a test of the null hypothesis \eqn{H_0: Var(X) <= v} against
\eqn{H_1 : Var(X) > v}. 
}

\usage{
npVarianceSingle(x, lower = 0, upper = 1, v,
      alternative = "greater", alpha = 0.05, iterations = 5000,
      epsilon = 1 * 10^(-6), ignoreNA = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a (non-empty) numeric vector of data values.
}
  \item{lower, upper}{
the theoretical lower and upper bounds on the data outcomes known ex-ante before gathering the data.
}
  \item{v}{
the value of the variance to be tested as \eqn{H_0: Var(x) >= v}.
}
  \item{alternative}{
a character string describing the alternative hypothesis, can take values "greater", "less" or "two.sided"
}
  \item{alpha}{
the type I error.
}
  \item{iterations}{
the number of iterations used, should not be changed if the exact solution should be derived.
}
  \item{epsilon}{
    the tolerance in terms of probability of the Monte Carlo simulations.
}
\item{ignoreNA}{
  if \code{TRUE}, NA values will be omitted. Default: \code{FALSE}
}
}
\details{
  This test randomly matches the data into pairs, then computes for each
  pair the square of the difference and continues with the resulting sequence with half as many observations as npMeanSingle. See the cited paper
for more information. 
}
\value{
A list with class "nphtest" containing the following components:

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
numerical estimate of the rejection probability of the randomized
    test, derived by taking an average of \code{iterations} realizations of the rejection probability. 
  }
  \item{bounds}{
    the lower and upper bounds of the variables.
  }
  \item{null.value}{
    the specified hypothesized value of the correlation between the variables.
  }
  \item{alpha}{
    the type I error.
  }
  \item{theta}{
    the parameter that minimizes the type II error.
  }
  \item{pseudoalpha}{
    \code{theta}*\code{alpha}, this is the level used when calculating
    the average rejection probability during the iterations.
  }
  \item{rejection}{
    logical indicator for whether or not the null hypothesis can be
    rejected.
  }
  \item{iterations}{
    the number of iterations that were performed.
  }
}
\references{
Karl Schlag (2008).
  Exact tests for correlation and for the slope in simple linear
  regressions without making assumptions. Available at \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1097}.  
}
\author{
Karl Schlag and Oliver Reiter
}


\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
\examples{
## see if the minority share holder shores have a variance greater
## than 0.05
data(mshscores)

scores <- as.vector(as.matrix(mshscores))
npVarianceSingle(scores, lower = 0, upper = 1, v = 0.05, ignoreNA = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variance test}
\keyword{single sample}
