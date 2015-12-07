
## print.nphtest
## a function to display the results of a nonparametric hypothesis test

## x ... object of class nphtest, a list containing:
##
## method ... a character string with the name of the test, character
## data.name ... name of the input variable(s), character
## null.value ... tested value, real
## probrej ... probability of rejection, real
## rejection ... if the null has been rejected, logical
## alpha ... type I error, real
## theta ... correction parameter, real
## pseudoalpha ... type I error in the randomized test, real
## bounds ... lower and upper bounds of the input variables, real
## estimate ... sample estimate of the statistic to be tested, real
## alternative ... the alternative hypothesis, character string with
##                 possible values (two.sided, greater, less)


print.nphtest <- function(x, digits = 4, prefix = "", ...)
{
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep="\n")
    cat("\n")
    cat("data: ", x$data.name, "\n")

    if(!is.null(x$alternative) && !is.null(x$null.hypothesis) && !is.null(x$alt.hypothesis))
    {
        ## stochin extra
        if(!is.null(x$stochin.parameter)) {
            cat("parameter: SI =", x$stochin.parameter, "\n")
            cat("\t   estimated value of SI:", x$stochin.estimate, "\n")
            cat("\n")
        }
        alt.char <- ifelse(x$rejection == TRUE,
                           "greater", "less")
        thus <- ifelse(x$rejection == TRUE,
                       "REJECTED", "NOT rejected")
        cat("null hypothesis H0: ", x$null.hypothesis, "\n", "\t",
            thus, " in favour of", "\n", sep = "")
        cat("alternative hypothesis H1: ", x$alt.hypothesis, "\n\n",
            sep = "")
        cat("as threshold probability:", x$probrej, "is", alt.char,
            "than theta:", x$theta, "\n")
      }

    if(!is.null(x$p.value))
      {
        fp <- format.pval(x$p.value, digits = digits)
        fp <- paste("p-value",
                         if(substr(fp, 1L, 1L) == "<") fp
                         else paste("=",fp))
        cat(fp, "\n")
    }

    cat("\ngiven parameters:\n")
    if(!is.null(x$bounds))
      {
        cat(paste("   ", x$data.name,
                  " in ", x$bounds,
                  "\n", sep = ""))
      }
    cat("   alpha:", x$alpha)
    cat("\n   theta:", x$theta)

    ## if(!is.null(x$pseudoalpha))
    ##   cat("\n   pseudoalpha:", x$pseudoalpha)
    if(!is.null(x$d.alternative))
      {
        cat("\n   d.alt:", x$d.alternative)
        cat("\n   typeII:", x$typeIIerror)
      }

    if(!is.null(x$iterations))
      cat("\n   iterations:", x$iterations)

    ## if(!is.null(x$conf.int)) {
    ##     cat(format(100 * attr(x$conf.int, "conf.level")),
    ##         "percent confidence interval:\n",
    ##         format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
    ## }
    if(!is.null(x$estimate)) {
        cat("\n\nsample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}


## print.npConfInt <- function(x, ..., verbose = NULL,
##                       digits = max(3L, getOption("digits") - 3L))
##     {
##         cat("\n")
##         cat(strwrap(x$method, prefix = "\t"), sep="\n")
##         cat("\n")
##         cat("data:", x$data.name, "in", x$bounds, "\n")
##         cat("n = ", x$n, "\n")

##         CI <- format(round(x$CI, digits = digits), digits = digits)

##         names(CI) <- c("lower", "upper", "Level")

##         cat("\n")
##         print(CI, quote = FALSE, right = TRUE)
##         ## cat("\nparameters:\n")
##         ## if(!is.null(x$parameter$bounds))
##         ##     {
##         ##         cat(paste("   ", x$yname, " in ", x$parameter$bounds, "\n", sep = ""))
##         ##     }
##         ## cat("   alpha:", x$parameter$alpha)
##         ## if(!is.null(x$parameter$iterations))
##         ##     {
##         ##         cat("\n   iterations:", x$parameter$iterations)
##         ##     }
##         cat("\n")

##     }

##     Confidence Interval with npMeanSingle

## data: x and y
## n = 50

##       lower    upper          test
##        0.65     0.98  npMeanSingle
