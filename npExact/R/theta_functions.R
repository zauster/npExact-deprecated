
###--------------------------------------------------
### Functions to derive the optimal theta
###--------------------------------------------------


## w
## helper function, to ease the reading of the code

w <- function(x)
  {
    as.numeric(x >= 0)
  }

## g1fun
## helper function, to calculate the type II error in
## npMeanSingleTypeIIerror

g1fun <- function(k, N, z, alpha)
  {
    summationterm1 <- alpha - (1 - pbinom(k - 1, N, z))
    summationterm2 <- alpha - (1 - pbinom(k, N, z))
    term3 <- summationterm2/dbinom(k, N, z)
    res <- w(summationterm1) + (1 - w(summationterm1)) *
  (w(summationterm2)) * term3
    res
  }

## g2fun
## helper function, to calculate the type II error in
## npMeanSingleTypeIIerror

g2fun <- function(mu, N, z, alpha)
  {
    k <- 0:(N - 1)
    term2 <- w(alpha - z^N) + (1 - w(alpha - z^N))*(alpha/z^N)
    res <- sum(dbinom(k, N, mu) * g1fun(k, N, z, alpha)) + mu^N*term2
    res
  }

g1 <- function(k,n,z,alpha)
#Calculates g1 term from formula
#k
#z
{
    j <- k:n
    sum1 <- choose(n,j)*(z^j)*(1-z)^(n-j)
    t1 <- as.numeric(alpha-sum(sum1)>=0)

    t2 <- 1-as.numeric(alpha-sum(sum1)>=0)
    j <- (k+1):n
    sum2 <-  choose(n,j)*(z^j)*(1-z)^(n-j)
    t3 <- as.numeric(alpha-sum(sum2)>=0)
    t4 <- (alpha-sum(sum2))/(choose(n,k)*(z^k)*(1-z)^(n-k))
    res <- t1+t2*t3*t4
    res
}

g2 <- function(miu,n,z,alpha)
#Calculates g2 term from formula
#k
#z
#uses g1
{
    k <- 0:(n-1)
    t1 <- sum(sapply(k, function(x)choose(n,x)*(miu^x)*(1-miu)^(n-x)*g1(x,n,z,alpha)),na.rm =TRUE)
    t2 <- miu^n*(as.numeric(alpha-z^n>=0)+(1-as.numeric(alpha-z^n>=0))*alpha/(z^n))
    res <- t1+t2
    res
}

## mean(replicate(100, system.time(g2fun(.4, 50, .6, .05))["elapsed"]))
## mean(replicate(100, system.time(g2(.4, 50, .6, .05))["elapsed"]))

## mean(replicate(100, system.time(g1fun(27, 50, .45, .05))["elapsed"]))
## mean(replicate(100, system.time(g1(27, 50, .45, .05))["elapsed"]))

## for(i in 1:100)
##   {
##     print(identical(round(g2(.4, i, .6, .05), 5),
##                     round(g2fun(.4, i, .6, .05), 5)))
##   }
## bei i in z: bei 1 nicht gleiches ergebnis


## for(i in 1:49/50)
##   {
##     print(identical(round(g1(24, 50, 0.5, i), 5),
##                     round(g1fun(24, 50, 0.5, i), 5)))
##   }

## fazit
## gleiches ergebnis, g2fun (und g1fun) schneller

possibleTheta <- function(n,p,alpha)
#Calculates possible values of theta, which are in interval (0,1)
{
    k <- 0:n
    ## j <- lapply(k, function(x)x:n)
    ## theta <- sapply(j,
    ##                 function(x)
    ##                 (1/alpha)*sum(choose(n,x)*p^x*(1-p)^(n-x)))
    ## r <- sapply(theta, function(x)if(x<1 && x>0) x else NA)
    ## res <- rbind(k, r)
    ## res1 <- res[, !is.na(res[2,])]

    theta <- (1/alpha) * pbinom(k - 1, n, p, lower = FALSE)
    sensibletheta <- theta < 1 & theta > 0
    res <- theta[sensibletheta]
    res <- rbind(k[sensibletheta], res)

    row.names(res) <- c("k","theta")
    res
    ## print(identical(round(res[2], 5), round(res1[2], 5)))
    ## list(res, res1)
}


minValue <- function(p1, p, n, alpha)
#Calculates minimum value, for given difference d
#uses possibleTheta, g2

{
    theta <- possibleTheta(n, p, alpha)
    f <- function(x)
      {
        g2(p1, n, p, alpha*x[2])
      }
    res <- apply(theta, 2,
                 function(x)(1-f(x))/(1-x[2]))
    res <- min(res, 1)

    righttheta <- theta[2, which(res == res1)]
    righttheta <- if(length(righttheta) == 0) NA
                  else r_theta
    ## list(theta = r_theta,
    ##      type2 = res)
    res - .5
}

## uniroot(minValue, c(0, 1), p = .1, n = 50, alpha = .05)
