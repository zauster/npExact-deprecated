\name{npMeanUnpaired}
\alias{npMeanUnpaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on two unpaired variables.
}
\description{
A test on the mean of two variables, where the observations do not have
to be paired.
}
\usage{
npMeanUnpaired(x1, x2,
               lower = 0, upper = 1, 
               iterations = 5000,
               alpha = 0.05,
               alternative = "greater")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{
the two variables to be tested.
}
  \item{lower, upper}{
the logical bounds of the variables.
}
%   \item{d}{
% difference (d = EY2 - EY1) such that type II is minimized when H1:
% E(Y1)+d <= E(Y2)
% }
  \item{iterations}{
the number of iterations
}
  \item{alpha}{
the type I error
}
%   \item{plotif}{
% logical, should a plot be generated?
% }
%   \item{helppoints}{
% the number of helping points for the plot.
% }
\item{alternative}{
the alternative hypothesis (currently unsupported)
}
}
\details{
Details to be added...
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
Karl Schlag (2008).
  Exact tests for correlation and for the slope in simple linear
  regressions without making assumptions. Available at \url{http://www.econ.upf.edu/en/research/onepaper.php?id=1097}.
}
\author{
Karl Schlag, Christian Pechhacker and Peter Saffert
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}

\examples{
npMeanUnpaired(runif(20), runif(20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
