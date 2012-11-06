\name{npVar}
\alias{npVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the variance of a variable.
}
\description{
This is a program that implements a nonrandomized one-sided test for
the variance of a single sample of independent observations.
}
\usage{
npVar(x, low = 0, up = 1, w, w1, alternative = "greater", alpha = 0.025,
      iterations = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the data vector to be tested.
}
  \item{low, up}{
the logical bounds of the variable
}
  \item{w}{
the value of the variance to be tested as H0: VarY >= w
}
  \item{w1}{
finetuning the type II error for H1: VarY <= w1
}
  \item{alternative}{
the alternative hypothesis, default = "greater"
}
  \item{alpha}{
the type I error
}
  \item{iterations}{
the number of iterations
}
}
\details{
More details to come...
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
  regressions without making assumptions
  
  
}
\author{
Karl Schlag and Oliver Reiter
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  http://www.econ.upf.edu/en/research/onepaper.php?id=1097
}
\examples{
npVar(runif(20), w = 0.1, w1 = 0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
