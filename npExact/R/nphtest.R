
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

    if(!is.null(x$alternative))
      {
        alt.char <- ifelse(x$rejection == TRUE,
                           "greater", "less")
        thus <- ifelse(x$rejection == TRUE,
                       "REJECTED", "NOT rejected")
        cat("Probability of rejection: ", x$probrej, "is", alt.char,
            "than theta:", x$theta, "\n\tthus", thus, "\n")
      }
    ## out <- character()
    ## if(!is.null(x$statistic))
    ##     out <- c(out, paste(names(x$statistic), "=",
    ##                      format(round(x$statistic, 4))))
    ## if(!is.null(x$parameter))
    ##     out <- c(out, paste(names(x$parameter), "=",
    ##                      format(round(x$parameter, 3))))
    if(!is.null(x$p.value))
      {
        fp <- format.pval(x$p.value, digits = digits)
        fp <- paste("p-value",
                         if(substr(fp, 1L, 1L) == "<") fp
                         else paste("=",fp))
        cat(fp, "\n")
    }

    if(!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if(!is.null(x$null.value)) {
            if(length(x$null.value) == 1L) {
                alt.char <-
                    switch(x$alternative,
                           two.sided = "not equal to",
                           less = "less than",
                           greater = "greater than")
                cat("true", names(x$null.value), "is", alt.char,
                    x$null.value, "\n")
            }
            else {
                cat(x$alternative, "\nnull values:\n")
                print(x$null.value, ...)
            }
        }
        else cat(x$alternative, "\n")
    }

    cat("\ngiven parameters:\n")
    if(!is.null(x$bounds))
      {
        cat(paste("   ", x$data.name,
                  " in ", x$bounds, sep = ""))
      }
    cat("\n   alpha:", x$alpha)
    cat("\n   theta:", x$theta)

    ## if(!is.null(x$pseudoalpha))
    ##   cat("\n   pseudoalpha:", x$pseudoalpha)

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

## Example:

##      Nonparametric Test on the Correlation

## data: Y1 and Y2
## Probability of rejection: ... is greater/lower than theta = ...
## (P-Value, where applicable)
## alternative hypothesis: true correlation is greater/lower than 0

## Given parameters:
##   Y1 in [lower, upper]
##   Y2 in [lower, upper]
##   alpha: ...
##   theta: ...
##   pseudoalpha: ...

## Sample estimate:
##   Covariance: ...
