\name{npExact-package}
\alias{npExact-package}
\alias{npExact}
\docType{package}
\title{Nonparametric hypothesis tests}
\description{
  \code{npExact} provides distribution-free hypothesis tests.
}
\details{
  \tabular{ll}{
    Package: \tab npExact\cr
    Type: \tab Package\cr
    Version: \tab 0.2 \cr
    Date: \tab 2013-09-26 \cr
    License: \tab GPL-2\cr
  }
  This package contains several new hypothesis tests, which
  do not require that the user makes assumptions on the underlying distributions.

  However, all tests except \code{npStochin} can only be applied if
  there are exogenously given bounds known to the user before gathering
  the data such that it is known by definition of the underlying process
  that all observations lie within these bounds. 

  So for instance, if the data involves percentages then the lower bound
  is 0 and the upper bound is 100, by definition of the data and not
  something (like normality) that cannot be deduced from the properties
  of the data.  
}
\author{
  Karl Schlag, Oliver Reiter, Peter Saffert, Christian Pechhacker, Simona
  Jokubauskaite, Tautvilas Janusauskas
}
\references{
  Karl Schlag, A New Method for Constructing Exact Tests without Making any Assumptions (August, 2008)
  Department of Economics and Business Working Paper 1109, Universitat Pompeu Fabra
}
%% \keyword{nonparametric, distribution-free}

\seealso{
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics/exacthypothesistesting8.pdf}
  
  \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
}
% \examples{
% }
