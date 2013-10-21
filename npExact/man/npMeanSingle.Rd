\name{npMeanSingle}
\alias{npMeanSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the mean of a single variable.
}
\description{
\code{npMeanSingle} provides an exact test for the mean of a single
sample based on Schlag (2008).
}
\usage{
npMeanSingle(x, mu, lower = 0, upper = 1, iterations = 5000,
             alpha = 0.05, alternative = "greater",
             epsilon = 1 * 10^(-6), ignoreNA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a (non-empty) numeric vector of data values.
}
  \item{mu}{
test value for \eqn{E(x) <= \mu}.
}
  \item{lower, upper}{
the theoretical lower and upper bounds on the data outcomes known ex-ante before gathering the data.
}
  \item{iterations}{
number of iterations.
}
  \item{alpha}{
    the type I error.
}
\item{alternative}{
    a character string describing the alternative hypothesis, can take values "greater", "less" or "two.sided".
}
  \item{epsilon}{
    the tolerance in terms of probability of the monte carlo simulations.
}
\item{ignoreNA}{
  if \code{TRUE}, NA values will be omitted. Default: \code{FALSE}
}
}
\details{
This test performs a test on the mean of a variable, of which only
exogenous bounds must be known, but not the distribution of the
variable.

With the known bounds, the data is transformed to lie in [0, 1]. Then,
0 and 1 are randomly assigned to each value, according to its
probability (which is its value after the transformation into [0,
1]). The exact binomial test then calculates the rejection probability
of this new data. These steps are repeated very often. If the fraction
of rejection is greater than \code{theta}, one can reject the null
hypothesis.

\code{theta} (and a \code{mu} in the set of the alternative hypothesis)
is found in an optimization procedure. \code{theta} and the
alternative \code{mu} are chosen as to bring the type II error to
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
    the type I error
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
  Schlag, Karl H. 2008, A New Method for Constructing Exact Tests
without Making any Assumptions, Department of Economics and Business
Working Paper 1109, Universitat Pompeu Fabra. Available at
  \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}. 
}

\author{
Karl Schlag, Peter Saffert and Oliver Reiter
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
\examples{
## test whether Americans gave more than 5 dollars in a round of
## the Ultimatum game
data(bargaining)
us_offers <- bargaining$US
npMeanSingle(us_offers, mu = 5, lower = 0, upper = 10, alternative =
"greater", ignoreNA = TRUE) ## no rejection

## test if the decrease in pain before and after the surgery is smaller
## than 50
data(pain)
pain$decrease <- with(pain, before - after)
without_pc <- pain[pain$pc == 0, "decrease"]
npMeanSingle(without_pc, mu = 50, lower = 0, upper = 100,
alternative = "less")
}
\keyword{single sample}
\keyword{mean test}
