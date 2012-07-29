\name{npMeanSingle}
\alias{npMeanSingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A test on the location of a single variable.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
npMeanSingle(data_vector, mu, mu1, low, up, times = 20000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_vector}{
The data vector to be tested.
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{mu1}{
%%     ~~Describe \code{mu1} here~~
}
  \item{low}{
%%     ~~Describe \code{low} here~~
}
  \item{up}{
%%     ~~Describe \code{up} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
npMeanSingle(runif(20), mu = 0.5, mu1 = 0.6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
