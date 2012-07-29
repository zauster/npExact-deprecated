\name{npMeanPaired}
\alias{npMeanPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A test on the mean difference of two variables, if the observations are paired.
}
\description{
  This test performs the exact McNemar Test on a random transformation of
  the data. The null hypothesis to be tested is H_0: E(x_1) >= E(x_2).
}
\usage{
npMeanPaired(var1, var2, low = 0, up = 1, d, alpha = 0.05,
             alternative = "greater", iter = 30000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var1, var2}{
    the numerical vectors which contain the variables to be tested.
  }
  \item{low, up}{
    the logical bounds of the variables
  }
  \item{d}{
    the hypothesised mean difference, which will be used to minimize the
    type II error.
  }
  \item{alpha}{
    the type I error
  }
  \item{alternative}{
    a character string describing the alternative hypothesis.
  }
  \item{iter}{
    Number of iterations.
  }
}
\details{
  This test uses the known bounds of the variables to transform the data
  into [0, 1]. Then a random transformation is used to
  turn the data into a binary-valued variables. On this variables the
  exact McNemar Test with level 'pseudoalpha' is performed and the
  result recorded. The last steps are repeated very often and if the
  proportion of rejected iterations is at least 'theta', then the null
  hypothesis is rejected.
}
\value{
  A list with class '"nphtest"' containing the following components:

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
  Paper 1109, Universitat Pompeu Fabra
}
\author{
  Karl Schlag and Christian Pechhacker
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
npMeanPaired(runif(20), runif(20), low = 0, up = 1, d = 0.4,
             alternative = "greater", iterations = 2000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
