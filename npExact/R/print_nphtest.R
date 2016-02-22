print.nphtest <- function(x, digits = 4, prefix = "", ...)
{
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep="\n")
    cat("\n")
    cat("data:", x$data.name, "\n")

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

        ## continue only if a valid theta was found
        if(!is.null(x$theta)) {
            cat("as threshold probability:", x$probrej, "is", alt.char,
                "than theta:", x$theta, "\n")
        } else {
            cat("given parameters:\n")
            if(!is.null(x$bounds))
            {
                cat(paste("   ", x$data.name,
                          " in ", x$bounds,
                          "\n", sep = ""))
            }
            cat("   alpha:", x$alpha, "\n")

            cat("\nNote:\n")
            cat("The sample is so small that the null hypothesis will not be rejected for any values in the sample (and theta cannot be calculated).\n")
            
            if(!is.null(x$estimate)) {
                cat("\nsample estimates:\n")
                print(x$estimate, ...)
            }
        }
      }

    if(!is.null(x$theta)) {
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
        cat("\n   theta:", x$thetaValue)

        if(!is.null(x$d.alternative))
        {
            cat("\n   d.alt:", x$d.alternative)
            cat("\n   typeII:", x$typeIIerror)
        }

        if(!is.null(x$iterations))
            cat("\n   iterations:", x$iterations)

        if(!is.null(x$mc.error))
            cat("\n   max error prob of MC:", round(x$mc.error, digits = digits))
        
        if(!is.null(x$estimate)) {
            cat("\n\nsample estimates:\n")
            print(x$estimate, ...)
        }
    }
    cat("\n")
    invisible(x)
}
