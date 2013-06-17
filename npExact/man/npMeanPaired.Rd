\name{npMeanPaired}
\alias{npMeanPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A test on the mean difference of two variables, if the observations are paired.
}
\description{
  This test performs the exact McNemar Test on a random transformation of
  the data. The null hypothesis to be tested is H_0: E(x_1) >= E(x_2),
  with the corresponding alternative hypothesis H_1: E(x_2) < E(x_2) or
  E(x_2) - E(x_1) > 0
}
\usage{
npMeanPaired(x1, x2, lower = 0, upper = 1, alpha = 0.05,
             alternative = "greater", epsilon = 1 * 10^(-6),
             iterations = 5000, completeCases = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{
    the numerical vectors which contain the variables to be tested.
  }
  \item{lower, upper}{
    the logical bounds of the variables
  }
  % \item{d}{
  %   the hypothesised mean difference, which will be used to minimize the
  %   type II error.
  % }
  \item{alpha}{
    the type I error
  }
  \item{alternative}{
    a character string describing the alternative hypothesis.
  }
\item{iterations}{
    Number of iterations.
  }
  \item{epsilon}{
    the algorithms keeps repeating until \eqn{exp(-2 * \textrm{iterations} *
    (\textrm{rejectionProbability} - \textrm{theta})^2)} is smaller than \code{epsilon}
}
\item{completeCases}{
  if \code{TRUE}, only Pairs of the two data vectors, which are
  complete, are being used for the computation. Default \code{FALSE}.
}

}
\details{
  This test uses the known bounds of the variables to transform the data
  into [0, 1]. Then a random transformation is used to
  turn the data into binary-valued variables. On this variables the
  exact McNemar Test with level \code{pseudoalpha} is performed and the
  result recorded. The last steps are repeated very often and if the
  proportion of rejected iterations is at least \code{theta}, then the null
  hypothesis is rejected.
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
  Schlag, Karl H. 2008, A New Method for Constructing Exact Tests without
  Making any Assumptions, Department of Economics and Business Working
  Paper 1109, Universitat Pompeu Fabra. Available at
  \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}. 
}
\author{
  Karl Schlag and Christian Pechhacker
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
\examples{
data(uncertainty)
first <- uncertainty$w1
second <- uncertainty$w2
npMeanPaired(first, second, lower = 0, upper = 60, alternative = "less")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
