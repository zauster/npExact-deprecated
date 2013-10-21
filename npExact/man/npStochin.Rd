\name{npStochin}
\alias{npStochin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A test for stochastic inequality.
}
\description{
  This test performs the exact McNemar test on a bootstrapped sample of
  the data. The corresponding null hypothesis for this test is \eqn{H_0:
    P(x_2 > x_1) <= P(x_2 < x_1)}, or equivalently: \eqn{P(x_2 > x_1) -
    P(x_2 < x_1) <= d}. 
}
\usage{
npStochin(x1, x2, d = 0, alternative = "greater", iterations = 5000,
          alpha = 0.05, epsilon = 1 * 10^(-6), ignoreNA = FALSE)
}
\arguments{
  \item{x1, x2}{
    the (non-empty) numerical data vectors which contain the variables to be tested.
  }
  \item{d}{
    the difference in probabilities of the null hypothesis: \eqn{P(x_2 >
      x_1) - P(x_2 < x_1) <= d}. Default is 0. 
  }
  \item{alternative}{
    a character string describing the alternative hypothesis. Default is
    "greater". If "less" is given, \code{x1} and \code{x2} are switched for each other.
  }
  \item{iterations}{
    Number of iterations.
  }
  \item{alpha}{
    the type I error.
  }
  \item{epsilon}{
    the tolerance in terms of probability of the monte carlo simulations.
}
\item{ignoreNA}{
  if \code{TRUE}, NA values will be omitted. Default: \code{FALSE}
}  
}
\details{
  The function draws two random samples (with the length of the shorter
  vector) of the data vectors and sees them as matched pairs. These
  matched pairs are then compared and the occurrences of \eqn{x_1 > x_2} and
  \eqn{x_1 < x_2} are counted. The exact binomial test then tests whether
  the success probability of \eqn{x_1 > x_2} is greater than 0.5 (in the
  default case). These steps are then repeated \code{iterations} times
  and rejection at level \code{alpha} times \code{theta} is recorded (for
  \code{theta}, see below). If the rejection probability (the fraction
  of rejections in all iterations) is above the threshold \code{theta},
  the null hypothesis is rejected.

  \code{theta} (and a \code{d} in the set of the alternative hypothesis)
  is found in an optimization procedure. \code{theta} and the
  alternative \code{d} are chosen so as to bring the type II error to
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
    an estimate of \eqn{P(x_2 > x_1) - P(x_2 < x_1)}.
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
  Schlag, Karl H. 2008, A New Method for Constructing Exact Tests without
  Making any Assumptions, Department of Economics and Business Working
  Paper 1109, Universitat Pompeu Fabra. Available at
  \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}. 
  
}
\author{
  Karl Schlag and Peter Saffert
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
data(french)
x <- french[, 1]
y <- french[, 2]
npStochin(x, y, ignoreNA = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{unpaired data}
\keyword{stochastic inequality}% __ONLY ONE__ keyword per line
