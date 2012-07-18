\name{print.nphtest}
\alias{print.nphtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A method to display the results of the hypothesis tests beautifully.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{nphtest}(x, digits = 4, prefix = "", ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{prefix}{
%%     ~~Describe \code{prefix} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, digits = 4, quote = TRUE, prefix = "", ...) 
{
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data: ", x$data.name, "\n")
    if (!is.null(x$alternative)) {
        alt.char <- ifelse(x$rejection == TRUE, "greater", "less")
        thus <- ifelse(x$rejection == TRUE, "REJECTED", "NOT rejected")
        cat("Probability of rejection: ", x$probrej, "is", alt.char, 
            "than theta:", x$theta, "\n\tthus", thus, "\n")
    }
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        fp <- paste("p-value", if (substr(fp, 1L, 1L) == "<") 
            fp
        else paste("=", fp))
        cat(fp, "\n")
    }
    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative, two.sided = "not equal to", 
                  less = "less than", greater = "greater than")
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
    if (!is.null(x$bounds)) {
        cat(paste("   ", x$data.name, " in ", x$bounds, sep = ""))
    }
    cat("\n   alpha:", x$alpha)
    cat("\n   theta:", x$theta)
    if (!is.null(x$pseudoalpha)) 
        cat("\n   pseudoalpha:", x$pseudoalpha)
    if (!is.null(x$iterations)) 
        cat("\n   iterations:", x$iterations)
    if (!is.null(x$estimate)) {
        cat("\n\nsample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
