

##' Amount sent in the Ultimatum Game
##' 
##' The Ultimatum game was played separately in four different countries. This
##' data contains the offers of 30 students in Israel and 27 in the United
##' States on a scale from 0 to 10. This dataset is taken from Roth et al.
##' (1991).
##' 
##' 
##' @name bargaining
##' @docType data
##' @format A data frame containing 30 observations for Israel and 27 for the
##' US.
##' @references Roth, A. E., Prasnikar, V., Okuno-Fujiwara, M., & Zamir, S.
##' (1991). Bargaining and market behavior in Jerusalem, Ljubljana, Pittsburgh,
##' and Tokyo: An experimental study. The American Economic Review, 1068-1095.
##' @keywords datasets
NULL





##' Indices of minority shareholder protection of countries with civil law with
##' and without french origin.
##' 
##' This data contains the indices of minority shareholder protection on a
##' scale from 0 to 1 in 51 countries with civil law, differentiating between
##' those with (32 observations) and those without (19 observations) french
##' origin. A higher value of the index means that country is more protected.
##' The data set is taken from Djankov et al. (2008).
##' 
##' 
##' @name french
##' @docType data
##' @format A dataframe containing 32 observations of countries with french
##' origin and 19 countries without french origin.
##' @references Djankov, S., La Porta, R., Lopez-de-Silanes, F., & Shleifer, A.
##' (2008). The law and economics of self-dealing. Journal of financial
##' economics, 88(3), 430-465.
##' @keywords datasets
NULL





##' Indices of minority shareholder protection of countries with common and
##' with civil law.
##' 
##' This data contains the indices of minority shareholder protection on a
##' scale from 0 to 1 in 51 countries with civil law and 21 countries with
##' common loaw. A higher value of the index means that country is more
##' protected. The data set is taken from Djankov et al. (2008).
##' 
##' 
##' @name mshscores
##' @docType data
##' @format A dataframe containing 51 observations for civil law and 21 for
##' common law.
##' @references Djankov, S., La Porta, R., Lopez-de-Silanes, F., & Shleifer, A.
##' (2008). The law and economics of self-dealing. Journal of financial
##' economics, 88(3), 430-465.
##' @keywords datasets
NULL





##' Nonparametric hypothesis tests
##' 
##' \code{npExact} provides distribution-free hypothesis tests.
##' 
##' This package contains several new hypothesis tests, which do not
##' require that the user makes assumptions on the underlying
##' distributions.
##' 
##' However, all tests except \code{npStochin} can only be applied if
##' there are exogenously given bounds known to the user before
##' gathering the data such that it is known by definition of the
##' underlying process that all observations lie within these bounds.
##' 
##' So for instance, if the data involves percentages then the lower
##' bound is 0 and the upper bound is 100, by definition of the data
##' and not something (like normality) that cannot be deduced from the
##' properties of the data.
##' 
##' @name npExact-package
##' @aliases npExact-package npExact
##' @docType package
##' @author Karl Schlag, Oliver Reiter, Peter Saffert, Christian Pechhacker,
##' Simona Jokubauskaite, Tautvilas Janusauskas
##' @seealso
##' \url{http://homepage.univie.ac.at/karl.schlag/research/statistics/exacthypothesistesting8.pdf}
##' 
##' \url{http://homepage.univie.ac.at/karl.schlag/research/statistics.html}
##' @references Karl Schlag, A New Method for Constructing Exact Tests without
##' Making any Assumptions (August, 2008) Department of Economics and Business
##' Working Paper 1109, Universitat Pompeu Fabra
##' @examples
##' 
##' 
##' ## npMeanPaired
##' ## test whether pain after the surgery is less than before the surgery
##' data(pain)
##' npMeanPaired(pain$before, pain$after, lower = 0, upper = 100)
##' 
##' ## npMeanSingle
##' ## test whether Americans gave more than 5 dollars in a round of
##' ## the Ultimatum game
##' data(bargaining)
##' us_offers <- bargaining$US
##' npMeanSingle(us_offers, mu = 5, lower = 0, upper = 10, alternative =
##' "greater", ignoreNA = TRUE) ## no rejection
##' 
##' ## npMeanUnpaired
##' ## test whether countries with french origin score lower than
##' ## countries with no french origin
##' data(french)
##' npMeanUnpaired(french[,1], french[,2], alternative = "less", ignoreNA =
##' TRUE)
##' 
##' ## npStochin
##' data(french)
##' x <- french[, 1]
##' y <- french[, 2]
##' npStochinUnpaired(x, y, ignoreNA = TRUE)
##' 
##' ## npVarianceSingle
##' ## see if the minority share holder shores have a variance greater
##' ## than 0.05
##' data(mshscores)
##' scores <- as.vector(as.matrix(mshscores))
##' npVarianceSingle(scores, lower = 0, upper = 1, v = 0.05, ignoreNA = TRUE)
##' 
NULL





##' Pain experienced before and after a knie operation
##' 
##' There are two ways to determine where to start an operation on a knee,
##' either with a computer or manually. The data describes the pain experienced
##' by the patients before and after the surgery.
##' 
##' 
##' @name pain
##' @docType data
##' @format A dataframe containing 50 observations. Column "pc" indicates if a
##' computer was used (coded with "1") or not (coded with "0")
##' @references Sabeti-Aschraf, M., Dorotka, R., Goll, A., & Trieb, K. (2005).
##' Extracorporeal shock wave therapy in the treatment of calcific tendinitis
##' of the rotator cuff. The American journal of sports medicine, 33(9),
##' 1365-1368.
##' @keywords datasets
NULL





##' Uncertainty in a game theoretical experiment.
##' 
##' In an experiment, subjects played a similar game twice. Choices could be
##' between 110 and 170. Each time, before they made their own choice, they had
##' to indicate an interval [L, U] that they believed would contain the choice
##' of their opponent. They paid some additional money if the choice of their
##' opponent was in the interval they specified, and were paid more the smaller
##' this interval was. So the width W_i of this interval in round i gives an
##' indication of how uncertain they are in round i. The data contains the
##' interval width in round 1 and 2 which makes this a sample of matched pairs.
##' 
##' 
##' @name uncertainty
##' @docType data
##' @format A dataframe containing the 25 intervals in each round of the game.
##' @references Galbiati, R., Schlag, K., & van der Weele, J. Sanctions that
##' Signal: an Experiment. Journal of Economic Behavior and Organization,
##' Forthcoming
##' @keywords datasets
NULL



