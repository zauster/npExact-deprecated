
## does a one-sided test for a single variable
## i.e. for the npMeanSingle or the npVarianceSingle test

doOneVariableTest <- function(alpha, epsilon,
                              iterations, max.iterations,
                              testFunction, p, n, ...)
{
    testFunction <- match.fun(testFunction)
    dots <- list(...)
   
    error <- 1
    rejMatrix <- vector(mode = "numeric", length = 0)

    ## cat("\np: ", p)
    ## cat("\nn: ", n)

    tryRes <- try(optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                           c(0, 1), p = p, N = n,
                                           alpha = alpha - epsilon),
                  silent = TRUE)
    if(inherits(tryRes, "try-error")) {
        
        ## pick up an error in the theta calculation
        ## and return a non-rejection

        results <- list(probrej = 0,
                        rejection = FALSE,
                        alpha = alpha,
                        theta = NULL,
                        d.alternative = NULL,
                        typeIIerror = NULL,
                        iterations.taken = 1000,
                        pseudoalpha = NULL)
        
    } else {

        theta <- minTypeIIError(optimaltypeII[[1]],
                                p = p, N = n,
                                alpha = alpha - epsilon)
        pseudoalpha <- alpha * theta$theta

        while(error > epsilon & length(rejMatrix) <= max.iterations) {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     testFunction(p = p, n = n,
                                                  pseudoalpha = pseudoalpha, dots)))
            rej <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rej - theta$theta)^2)
        }
        
        results <- list(probrej = rej,
                        rejection = ifelse(rej >= theta$theta, TRUE, FALSE),
                        alpha = alpha,
                        theta = theta$theta,
                        d.alternative = optimaltypeII$root,
                        typeIIerror = theta$typeII,
                        mc.error = error,
                        iterations.taken = length(rejMatrix), 
                        pseudoalpha = pseudoalpha)

    }

    return(results)
}
