\name{npVarianceSingle}
\alias{npVarianceSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the variance of a variable.
}
\description{
This is a program that implements a nonrandomized one-sided test for
the variance of a single sample of independent observations. The null
hypothesis that is being tested is \eqn{H_0: Var(x) >= w} against
\eqn{H_1: Var(x) < w}
}

\usage{
npVarianceSingle(x, lower = 0, upper = 1, variance,
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
  \item{variance}{
the value of the variance to be tested as \eqn{H_0: Var(x) >= w}.
}
  \item{alternative}{
a character string describing the alternative hypothesis, can take values "greater", "less" or "two.sided"
}
  \item{alpha}{
the type I error.
}
  \item{iterations}{
the number of iterations.
}
  \item{epsilon}{
    the tolerance in terms of probability of the monte carlo simulations.
}
\item{ignoreNA}{
  if \code{TRUE}, NA values will be omitted. Default: \code{FALSE}
}
}
\details{
  This function does the following steps in each iteration: First, it
  permutes the data vector. Then it calculates pairwise
  differences and squares them: \eqn{(x_{1:(n/2)} - x_{(n/2 +
      1):n})^2}. These values all lie in [0, 1] and are then randomly
  transformed into {0, 1}. A binomial test with \eqn{p = 2 *
    \textrm{variance} / (\textrm{upper} - \textrm{lower})^2} is then
  performed on the data. If the fraction of recorded rejections lies
  above the value \code{theta}, the null hypothesis can be rejected.

  \code{theta} (and a \code{w} in the set of the alternative hypothesis)
is found in an optimization procedure. \code{theta} and the
alternative \code{w} are chosen as to bring the type II error to
0.5. Please see the cited paper below for further information.
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
    the probability with which the null hypothesis can be rejected.
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
    \code{theta} * \code{alpha}, this is the alpha that is used during the iterations.
  }
  \item{rejection}{
    Logical. Given \code{alpha}, the type I error, and \code{theta}, can a rejection be achieved?
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
npVarianceSingle(scores, lower = 0, upper = 1, var = 0.05, ignoreNA = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variance test}
\keyword{single sample}
