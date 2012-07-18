#######################################################################
##   Program Name   npMeanPaired                                      #
##   Purpose        Exact nonparametric test of Schlag (2008)         #
##                  for comparing expected values given matched pairs #
##                  for H0: E(Y1) >= E(Y2)  if Y1,Y2 in [low,up]      #
##                                                                    #
##   Input variables                                                  #
##        required:   y1       independent sample of Y1               #
##                    y2       independent sample of Y2               #
##                    low      exogenously known lower                #
##                             bound for any outcome that can         #
##                             be generated in either sample          #
##                    up       exogenously known upper                #
##                             bound for any outcome that can         #
##                             be generated in either sample          #
##                    d        difference (d= EY2 - EY1)              #
##                             s.t. type II is minimized              #
##                             when H1: E(Y1)+d <= E(Y2)              #
##                             so 0< d <= up-low                      #
##                                                                    #
##        optional (default): times  (30000) number of Monte Carlo    #
##                                   iterations                       #
##                            alpha  (0.05) significance level of test#
##                                                                    #
##   Return Values                                                    #
##                            obs    number of observations in sample #
##                            low    lower bound for data values      #
##                            up     upper bound for data values      #
##                            theta  cutoff level theta               #
##                            alpha  significance level alpha         #
##                            avg1   average of values in y1          #
##                            avg2   average of values in y2          #
##                            P(rj)  probability of rejection in the  #
##                                   under lying randomdomized test   #
##                                                                    #
##   Authors           Christian Pechhacker and Karl Schlag           #
##   Date              04.03.2012                                     #
##                                                                    #
##                                                                    #
## For documentation see                                              #
##   (Schlag, Karl H. 2008, A New Method for Constructing Exact       #
##   Tests without Making any Assumptions, Department of              #
##   Economics and Business Working Paper 1109, Universitat           #
##   Pompeu Fabra)                                                    #
##                                                                    #
#######################################################################


## Examples:

## y1 <- sample(c(1,2), size = 100, replace = TRUE)
## y2 <- sample(c(1,2,3,4), size = 100, replace = TRUE)
## npMeanPaired(y1, y2, low = 0, up = 5, d=1.27, alpha = 0.05)
## npMeanPaired(runif(20), runif(20), low = 0, up = 1, d = 0.1,
##              alternative = "greater", iter = 2000)


npMeanPaired <- function(var1, var2, low = 0, up = 1,
                         d, alpha = 0.05,
                         alternative = "greater", iter = 30000)
  {
    DNAME <- paste(deparse(substitute(var1)), "and",
                   deparse(substitute(var2)))

    var1 <- as.vector(var1)
    var2 <- as.vector(var2)

    if(min(var1) < low)
      stop("Your variable includes values below the specified lower bound (low). This contradicts the theoretical requirements of this test.\n")
    if(max(var2) > up)
      stop("Your variable includes values above the specified upper bound (up). This contradicts the theoretical requirements of this test.\n")

    n <- length(var1)
    mean1 <- mean(var1)
    mean2 <- mean(var2)


#### Normalize vectors to [0,1]
    var1 <- (var1 - low)/(up - low)
    var2 <- (var2 - low)/(up - low)
    d <- d/(up - low)

    ## it <- as.numeric(min_value(n=n, p=1/2, p1=1/2+d/2, alpha=alpha))
    ## print(it)
    ## theta <- it[1]
    theta <- 0.4
    pseudoalpha <- alpha * theta

    res <- replicate(iter, doMatching(randTransformation(var1, n),
                                      randTransformation(var2, n),
                                      pseudoalpha))
    rej <- mean(res)

    method <- "Nonparametric Paired Mean Test"
    sample.est <- c(mean(var1), mean(var2))
    names(sample.est) <- c("mean", "mean")
    null.value <- 0
    names(null.value) <- "mean difference"
    rejection <- ifelse(rej > theta, TRUE, FALSE)
    bounds <- paste("[", low, ", ", up, "]", sep = "")

    structure(list(method = method,
                   data.name = DNAME,
                   alternative = alternative,
                   estimate = sample.est,
                   probrej = rej,
                   rejection = rejection,
                   alpha = alpha,
                   theta = theta,
                   iterations = iter,
                   pseudoalpha = pseudoalpha,
                   bounds = bounds,
                   null.value = null.value),
              class = "nphtest")
  }


McNemarTestRandom <- function(n10, n01, alpha)
  {
    #### performs the randomized McNemar test
    #### n10, n01 ... counts of (1,0) and (0,1) respectively
    #### returns either 1 (rejection), p (prob of rejection) or 0 (no
    #### rejection)

    k <- n01:(n10 + n01)

    prob <- sum((gamma(n10 + n01 + 1)/(gamma(k + 1) * gamma(n10 + n01 - k +
                                                            1)))/(2^(n10 + n01)))
    #### print(prob)

    res <- 0
    if(prob <= alpha)
      {
        res <- 1
      }
    else
      {
        h <- gamma(n10 + n01 + 1)/(gamma(n01 + 1) * gamma(n10 +1))
        if(prob <= alpha)
          {
            res <- (alpha - prob + h)/h
          }
      }
    return(res)
  }

doMatching <- function(var1, var2, alpha)
  {
    #### takes two vectors and does the matching and testing

    n10 <- sum(var1 > var2)
    n01 <- sum(var1 < var2)

    McNemarTestRandom(n10, n01, alpha)

  }

randTransformation <- function(x, n)
  {
    #### this function is used to generate the random transformations of
    #### the data vector (which is in [0,1] -> {0,1})

    #### rand <- runif(n)
    #### res <- (rand < x)
    #### return(as.numeric(res))
    #### or as a one-liner:
    return(as.numeric(runif(n) < x))

  }
