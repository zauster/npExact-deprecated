
## source("npMeanPaired.R")
## source("npMeanSingle.R")
## source("npMeanUnpaired.R")
## source("npVarianceSingle.R")
## source("npStochin.R")
## source("theta_functions.R")
## source("nphtest.R")

## ## x <- runif(20)
## ## npConfInt(x, lower = 0, upper = 1)

## ## Computes confidence intervals for a given result of a test

## npConfInt <- function(x, y = NULL, lower = 0, upper = 1, conf.level = 0.95,
##                       bound = "lower", # or "upper", "two.sided"
##                       test = npMeanSingle, ...)
##     {
##         DNAME <- deparse(substitute(x))
##         if(!is.null(y))
##             DNAME <- paste(DNAME, "and", deparse(substitute(y)))

##         method <- paste("Confidence Interval with", deparse(substitute(test)))
##         bounds <- paste("[", round(lower, digits = 3), ", ",
##                         round(upper, digits = 3), "]", sep = "")


##         ## Single variable CIs
##         if(is.null(y))
##             {
##                 upperint <- upper
##                 epsilon <- (upper - lower)/sqrt(length(x))

##                 if(deparse(substitute(test)) == "npVarianceSingle")
##                     {
##                         upperint <- 0.25*(upper - lower)^2
##                         epsilon <- (upperint - lower)/length(x)
##                     }

##                 cat("epsilon: ", epsilon, "\tupperint: ", upperint, "\n")
##                 ## for two-sided testing
##                 if(bound == "two.sided")
##                     {
##                         alpha <- (1 - conf.level)/2

##                         lowerbound <- uniroot(npSingleVarWrapper,
##                                        interval = c(lower + epsilon, upperint - epsilon),
##                                        x = x, lower, upper,
##                                        alpha = alpha, test = test,
##                                        alternative = "greater", ...)$root
##                         upperbound <- uniroot(npSingleVarWrapper,
##                                        interval = c(lower + epsilon, upperint - epsilon),
##                                        x = x, lower, upper,
##                                        alpha = alpha, test = test,
##                                        alternative = "less", ...)$root
##                         CI <- c(lowerbound, upperbound, conf.level)
##                     }
##                 else
##                     {
##                         alpha <- 1 - conf.level
##                         alternative <- ifelse(bound == "lower", "greater", "less")
##                         res <- uniroot(npSingleVarWrapper,
##                                        interval = c(lower + epsilon, upperint - epsilon),
##                                        x = x, lower, upper,
##                                        alpha = alpha, test = test,
##                                        alternative, ...)$root
##                         if(bound == "lower")
##                             {
##                                 CI <- c(res, upperint, conf.level)
##                             }
##                         else
##                             {
##                                 CI <- c(lower, res, conf.level)
##                             }
##                     }

##             }

##         structure(list(method = method,
##                        data.name = DNAME,
##                        n = length(x),
##                        bounds = bounds,
##                        CI = CI),
##                   class = "npConfInt")

##   ## structure(list(method = method,
##   ##                data.name = DNAME,
##   ##                alternative = alternative,
##   ##                null.hypothesis = null.hypothesis,
##   ##                alt.hypothesis = alt.hypothesis,
##   ##                estimate = sample.est,
##   ##                probrej = rej,
##   ##                rejection = rejection,
##   ##                alpha = alpha,
##   ##                theta = theta$theta,
##   ##                d.alternative = optimaltypeII$root,
##   ##                typeIIerror = theta$typeII,
##   ##                iterations = length(rejMatrix),
##   ##                pseudoalpha = pseudoalpha,
##   ##                bounds = bounds,
##   ##                null.value = null.value),
##   ##           class = "nphtest")

##     }

## npSingleVarWrapper <- function(nullvalue, x, lower, upper,
##                                alpha = 0.025, test = npMeanSingle,
##                                alternative, ...)
##     {
##         res <- test(x = x, nullvalue, lower, upper, alternative, alpha = alpha,
##                     max.iterations = 25000, ...)
##         res <- res$probrej - res$theta

##         print("ITERITERITERITERITERITERITERITERITER")
##         print(nullvalue)

##         return(res)
##     }











## ## npMultVarWrapper <- function(nullvalue, x1, x2, alpha = 0.025, test = npMeanPaired, ...)
## ##     {
## ##         res <- test(x1, x2, ...)
## ##     }

## ## npMeanPairedWrapper <- function(nullvalue, x1, x2, alpha = 0.025, ...)
## ##     {
## ##         res <- npMeanPaired(x1, x2 + nullvalue,
## ##                             lower = lower + nullvalue,
## ##                             upper = upper + nullvalue, ...)

## ##         res <- res$probrej - res$theta

## ##         return(res)
## ##     }
