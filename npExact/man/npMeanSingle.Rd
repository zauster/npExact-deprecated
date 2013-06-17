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
npMeanSingle(x, mu, lower = 0, upper = 1,
             iterations = 5000, alpha = 0.05,
             alternative = "greater",
             epsilon = 1 * 10^(-6),
             ignoreNA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The data vector to be tested.
}
  \item{mu}{
test value for EY <= mu
}
  \item{lower, upper}{
exogenously known bounds for any outcome that can be generated in
either sample
}
  \item{iterations}{
number of iterations
}
  \item{alpha}{
type I error, level of the test
}
\item{alternative}{
the alternative hypothesis (currently unsupported)
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
data(pain)
pain$decrease <- with(pain, before - after)
without_pc <- pain[pain$pc == 0, "decrease"]
npMeanSingle(without_pc, mu = 50, lower = 0, upper = 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
