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
the data vector to be tested.
}
  \item{lower, upper}{
the logical bounds of the variable.
}
  \item{variance}{
the value of the variance to be tested as \eqn{H_0: Var(x) >= w}.
}
  \item{alternative}{
the alternative hypothesis, default = "greater"
}
  \item{alpha}{
the type I error.
}
  \item{iterations}{
the number of iterations.
}
  \item{epsilon}{
    the algorithms keeps repeating until \eqn{exp(-2 * \textrm{iterations} *
    (\textrm{rejectionProbability} - \textrm{theta})^2)} is smaller than \code{epsilon}
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
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
\examples{
data(mshscores)
npVarianceSingle(mshscores$civil, lower = 0, upper = 1, var = 0.01)
npVarianceSingle(mshscores$civil, lower = 0, upper = 1, var = 0.1, alternative = "less")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variance test}
\keyword{single sample}
