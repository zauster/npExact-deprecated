
doTwoVariablesTest <- function(alpha, epsilon,
                               theta = NULL, typeII = NULL,
                               d.alternative = NULL,
                               iterations, max.iterations,
                               testFunction, x1, x2, 
                               ...)
{
    dots <- list(...)
    
    error <- 1
    rejMatrix <- vector(mode = "numeric", length = 0)

    results <- NULL
    if(is.null(theta)) {
        if(deparse(substitute(testFunction)) != "randomTest") {
            
            ## npStochinUnpaired and npMeanPaired
            res <- try({
                optimaltypeII <- uniroot(minTypeIIErrorWrapper,
                                         c(0, 1), p = dots[["p"]],
                                         N = dots[["n"]],
                                         alpha = alpha - epsilon)
                
                theta <- minTypeIIError(optimaltypeII[[1]],
                                        p = dots[["p"]], N = dots[["n"]],
                                        alpha = alpha - epsilon)
                d.alternative <- optimaltypeII$root
            }, silent = TRUE)
        }
        else {

            ## print("here")

            ## npMeanUnpaired
            res <- try({
                optimaltypeII <- optimize(npMeanUnpairedminTypeIIErrorWrapper,
                                          c(0, 1),
                                          n1 = dots[["n1"]],
                                          n2 = dots[["n2"]],
                                          alpha = alpha - epsilon,
                                          tol = .Machine$double.eps^0.15)
                theta <- optimizeTheta(dots[["n1"]],
                                       dots[["n2"]],
                                       optimaltypeII$minimum,
                                       alpha - epsilon)
                
                d.alternative <- optimaltypeII$minimum
                ## print("here")
                ## print(theta)
                if(round(theta$theta, digits = 4) >= 1L | round(theta$typeII >= 1L)) {
                    stop("No valid theta found")
                }
            }, silent = TRUE)
        }
        
        ## pick up an error in the theta calculation
        if(inherits(res, "try-error")) {
            results <- list(probrej = 0,
                            rejection = FALSE,
                            alpha = alpha,
                            theta = NULL,
                            d.alternative = NULL,
                            typeIIerror = NULL,
                            iterations.taken = 1000,
                            pseudoalpha = NULL)
            
        } else { ## if everything worked out, can compute theta
            typeII <- theta$typeII
            theta <- theta$theta
        }
    }

    ## if there was no error, we can proceed with the test
    if(is.null(results)) {
        testFunction <- match.fun(testFunction)
        pseudoalpha <- alpha * theta

        ## calculate the probability of rejection
        while(error > epsilon & length(rejMatrix) <= max.iterations) {
            rejMatrix <- c(rejMatrix,
                           replicate(iterations,
                                     testFunction(x1, x2, pseudoalpha,
                                                  dots)))
            rej <- mean(rejMatrix)
            error <- exp(-2 * length(rejMatrix) * (rej - theta)^2)
        }
        
        results <- list(probrej = rej,
                        rejection = ifelse(rej >= theta, TRUE, FALSE),
                        alpha = alpha,
                        theta = theta,
                        d.alternative = d.alternative,
                        typeIIerror = typeII,
                        mc.error = error,
                        iterations.taken = length(rejMatrix), 
                        pseudoalpha = pseudoalpha)
    }

    return(results)
    
}
