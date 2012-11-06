\name{npMeanSingle}
\alias{npMeanSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the location of a single variable.
}
\description{
\code{npMeanSingle} provides an exact test for the mean of a single
sample based on Schlag (2008).
}
\usage{
npMeanSingle(x, mu, lower = 0, upper = 1,
             iterations = 5000, alpha = 0.05,
             alternative = "greater")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The data vector to be tested.
}
  \item{mu}{
Test value for EY <= mu
}
  \item{lower}{
exogenously known lower bound for any outcome that can be generated in
either sample
}
  \item{upper}{
exogenously known upper bound for any outcome that can be generated in
either sample
}
  \item{iterations}{
number of iterations
}
  \item{alpha}{
type I error
}
\item{alternative}{
the alternative hypothesis (currently unsupported)
}
}
\details{
Details missing currently ...
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
Working Paper 1109, Universitat Pompeu Fabra
}

\author{
Karl Schlag and Peter Saffert
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
x <- runif(30)
npMeanSingle(x, mu = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
