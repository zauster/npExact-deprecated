\name{npCorr}
\alias{npCorr}
\title{
A test on the correlation between two variables.
}
\description{
This test uses the exact Fisher-Tocher test to test whether the
correlation between two variables is greater or smaller than 0.
}
\usage{
npCorr(x1, x2, lower.x1 = 0, upper.x1 = 1, lower.x2 = 0, upper.x2 = 1,
       iterations = 5000, alpha = 0.05, alternative = "greater",
       conditional = FALSE)
}
\arguments{
  \item{x1, x2}{
Two numerical vectors which contain the variables to be tested.
}
  \item{lower.x1, upper.x1, lower.x2, upper.x2}{
The logical bounds of the two variables, whereas \code{lower.x1} and \code{lower.x2}
refer to the lower bounds and \code{upper.x2} and \code{upper.x2} to the upper bounds.
}
  \item{iterations}{
Number of iterations.
}
%   \item{theta}{
%     the parameter that minimizes the type II error.
% }
  \item{alpha}{
    the type I error
}
%   \item{pseudoalpha}{
%     theta * alpha, this is the alpha that is used during the iterations.
% }
  \item{alternative}{
    a character string describing the alternative hypothesis.
}
  \item{conditional}{
Logical, should the conditional or the unconditional version of the test
be used?
}
}
\details{
  This function carries out a randomized version of the Fisher-Tocher
  test on the covariance of two variables with binary outcomes. Using
  the lower and upper bounds, the data in [lower, upper] is transformed
  to [0, 1]. Then, using a random transformation, the data is
  transformed again to \{0, 1\}. On this data, the exact Fisher-Tocher
  test is performed with a level of 'pseudoalpha' and the rejection
  TRUE/FALSE is recorded. This last step is repeated 'iterations'
  times. If the proportion of rejected iterations is at least 'theta',
  then the null hypothesis is rejected.
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
  Schlag, Karl H. 2008, Exact tests for correlation and for the slope in
    simple linear regressions without making assumptions. Available at \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1097}.
}
\author{
Karl Schlag and Oliver Reiter
}
\note{
This function is still being tested.
}

\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}

\examples{
x1 <- runif(20)
x2 <- runif(20)
npCorr(x1, x2) ## should NOT reject
npCorr(sort(x1), sort(x2)) ## should reject
}
% Add one or more standard keywords, see file 'KExWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
