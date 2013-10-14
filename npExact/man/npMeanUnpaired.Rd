\name{npMeanUnpaired}
\alias{npMeanUnpaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the means of two unpaired variables.
}
\description{
A test on the mean of two variables, where the observations do not have
to be paired.
}
\usage{
npMeanUnpaired(x1, x2, lower = 0, upper = 1, iterations = 5000, alpha = 0.05,
               alternative = "greater", epsilon = 1 * 10^(-6), ignoreNA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{
    the (non-empty) numerical data vectors which contain the variables to be tested.
}
  \item{lower, upper}{
the logical bounds of the variables.
}
  \item{iterations}{
the number of iterations.
}
  \item{alpha}{
    the type I error.
}
  \item{alternative}{
    a character string describing the alternative hypothesis.    
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
This function performs a test on the means of two variables, \eqn{H_0:
  E(x_1) >= E(x_2)}, with the corresponding alternative hypothesis
\eqn{H_1:  E(x_1) < E(x_2)} or \eqn{E(x_2) - E(x_1) > 0}. 

  First, the data is transformed to lie in [0, 1]. This is done with the
  ex ante known bounds of the variables. Afterwards, the iteration phase
  takes place: In each iteration, each observation is randomly assigned
  a value {0, 1} (according to its probability, the value in [0, 1]). A
  randomized Fischer-Tocher-Test with level \code{pseudoalpha} is then
  applied to the data and its result recorded. If the fraction of the
  rejections is at least \code{theta}, then one can reject the null
  hypothesis.

  \code{theta} (and a \code{mu} in the set of the alternative
  hypothesis) is found in an optimization procedure. \code{theta} and the
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
    the sample means of the two variables.
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
Karl Schlag (2008), A New Method for Constructing Exact Tests without
Making any Assumptions. Available at \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1109}.
}
\author{
Karl Schlag, Christian Pechhacker, Peter Saffert and Oliver Reiter
}
\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}

\examples{
data(bargaining)
npMeanUnpaired(bargaining$US, bargaining$IS, lower = 0, upper = 10,
               alternative = "less", ignoreNA = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{unpaired data}
\keyword{mean test}% __ONLY ONE__ keyword per line
