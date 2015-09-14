
###--------------------------------------------------
### Functions to derive the optimal theta
###--------------------------------------------------


## w
## helper function, to ease the reading of the code

w <- function(x)
    {
        as.numeric(x >= 0)
    }

## g1
## helper function
## g1 <- function(k, N, z, alpha = 0.05)
##   {
##     summationterm1 <- alpha - pbinom(k - 1, N, z, lower.tail = FALSE)
##     summationterm2 <- alpha - pbinom(k, N, z, lower.tail = FALSE)
##     term3 <- summationterm2/dbinom(k, N, z)
##     res <- w(summationterm1) + (1 - w(summationterm1)) *
##   (w(summationterm2)) * term3
##     res
##   }

## second version is faster

g1 <- function(k, N, z, alpha = 0.05)
    {
        ## if alpha < term2 -> 0
        ## if alpha > term1 -> 1
        ## else calculate the fraction term
        term2 <- pbinom(k, N, z, lower.tail = FALSE)
        ifelse(alpha < term2, 0,
               ifelse(alpha > pbinom(k - 1, N, z, lower.tail = FALSE),
                      1,
                      (alpha - term2)/dbinom(k, N, z)))
    }

## g2
## helper function
g2 <- function(mu, N, z, alpha)
    {
        k <- 0:(N - 1)

        ## term2
        ## if alpha >= z^N -> 1
        ## if alpha < z^N -> alpha/z^N
        term1 <- ifelse(alpha >= z^N, 1, alpha/z^N)

        ## g1 is 0 until alpha >= term2 from above
        ## how to find k^* fast? -> could then sum from there!
        res <- sum(dbinom(k, N, mu) * g1(k, N, z, alpha)) + mu^N*term1
        res
    }

## g2 <- function(mu, N, z, alpha)
##     {
##         k <- 0:(N - 1)
##         term2 <- w(alpha - z^N) + (1 - w(alpha - z^N))*(alpha/z^N)
##         res <- sum(dbinom(k, N, mu) * g1(k, N, z, alpha)) + mu^N*term2
##         res
##     }

## possibleTheta
## Calculates possible values of theta, which are in interval (0,1)
possibleTheta <- function(N, p, alpha)
{
    k <- 0:N
    ## j <- lapply(k, function(x)x:N)
    ## theta <- sapply(j,
    ##                 function(x)
    ##                 (1/alpha)*sum(choose(N,x)*p^x*(1-p)^(N-x)))
    ## r <- sapply(theta, function(x)if(x<1 && x>0) x else NA)
    ## res <- rbind(k, r)
    ## res1 <- res[, !is.na(res[2,])]

    theta <- (1/alpha) * pbinom(k - 1, N, p,
                                lower.tail = FALSE)
    sensibletheta <- theta < 1 & theta > 0
    res <- theta[sensibletheta]
    res <- rbind(k[sensibletheta], res)

    row.names(res) <- c("k","theta")
    res
    ## print(identical(round(res[2], 5), round(res1[2], 5)))
    ## list(res, res1)
}


minTypeIIError <- function(p.alt, p, N, alpha)
    {
        ## Calculates minimum value, for given difference d
        ## uses possibleTheta, g2
        theta <- possibleTheta(N, p, alpha)
        ## f <- function(x)
        ##   {
        ##     g2(p.alt, n, p, alpha*x[2])
        ##   }
        ## typeIIerrors <- apply(theta, 2,
        ##              function(x)(1 - f(x))/(1 - x[2]))

        cat(p.alt, "\t", p, "\t", N, "\t", alpha, "\n")

        f <- function(x)
            {
                (1 - g2(p.alt, N, p, alpha * x))/(1 - x)
            }
        typeIIerrors <- sapply(theta[2,], f)

        print(typeIIerrors)
        print(str(typeIIerrors))

        if(!is.numeric(typeIIerrors))
            {
                stop("The null will not be rejected as the value of the mean under the null is too large or too small.")
            }

        mintypeII <- min(typeIIerrors, 1)

        righttheta <- theta[2, which(typeIIerrors == mintypeII)]
        righttheta <- ifelse(length(righttheta) == 0, NA, righttheta)

        list(theta = righttheta,
             typeII = mintypeII)
    }

minTypeIIErrorWrapper <- function(p.alt, p, N, alpha,
                                  typeIIgoal = .5)
    {
        minTypeIIError(p.alt, p, N, alpha)$typeII - typeIIgoal
    }

## uniroot(minTypeIIErrorWrapper, c(0, 1), p = .1, N = 50, alpha = .05)
