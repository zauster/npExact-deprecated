\name{uncertainty}
\docType{data}
\alias{uncertainty}
\title{Uncertainty in a game theoretical experiment.}
\description{
  In an experiment, subjects played a similar game twice. Choices could
  be between 110 and 170. Each time, before they made their own choice,
  they had to indicate an interval [L, U] that they believed would
  contain the choice of their opponent. They paid some additional money
  if the choice of their opponent was in the interval they specified,
  and were paid more the smaller this interval was. So the width W_i of
  this interval in round i gives an indication of how uncertain they are
  in round i. The data contains the interval width in round 1 and 2
  which makes this a sample of matched pairs. 
}
\usage{uncertainty}
\format{A dataframe containing the 25 intervals in each round of the game.}
\references{
  Galbiati, R., Schlag, K., & van der Weele, J. Sanctions that
  Signal: an Experiment. Journal of Economic Behavior and Organization,
  Forthcoming 
  }
\keyword{datasets}
