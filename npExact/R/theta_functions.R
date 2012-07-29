
###--------------------------------------------------
### Functions to derive the optimal theta
###--------------------------------------------------



g1 <- function(k, n, z, alpha) ## Calculates g1 term from formula k z
{
    j <- k:n
    sum1 <- choose(n, j) * (z^j) * (1 - z)^(n - j)
    t1 <- as.numeric(alpha - sum(sum1) >= 0)
    t2 <- 1 - as.numeric(alpha - sum(sum1) >= 0)
    j <- (k + 1):n
    sum2 <- choose(n, j) * (z^j) * (1 - z)^(n - j)
    t3 <- as.numeric(alpha - sum(sum2) >= 0)
    t4 <- (alpha - sum(sum2))/(choose(n, k) * (z^k) * (1 - z)^(n -
        k))
    res <- t1 + t2 * t3 * t4
    res
}

g2 <- function(miu, n, z, alpha) ## Calculates g2 term from formula k z uses g1
{
    k <- 0:(n - 1)
    t1 <- sum(sapply(k, function(x) choose(n, x) * (miu^x) *
        (1 - miu)^(n - x) * g1(x, n, z, alpha)), na.rm = TRUE)
    t2 <- miu^n * (as.numeric(alpha - z^n >= 0) + (1 - as.numeric(alpha -
        z^n >= 0)) * alpha/(z^n))
    res <- t1 + t2
    res
}


possible.theta <- function(n, p, alpha) ## Calculates possible values of theta, which are in
## interval (0,1)
{
    k <- 0:n
    j <- lapply(k, function(x) x:n)
    theta <- sapply(j, function(x) (1/alpha) * sum(choose(n,
        x) * p^x * (1 - p)^(n - x)))
    r <- sapply(theta, function(x) if (x < 1 && x > 0)
        x else NA)
    res <- rbind(k, r)
    res <- res[, !is.na(res[2, ])]
    row.names(res) <- c("k", "theta")
    res
}


min_value <- function(n, p, p1, alpha)
  ## Calculates minimum value, for given difference d uses
  ## possible.theta, g2
{
    theta <- possible.theta(n, p, alpha)
    f <- function(x) g2(p1, n, p, alpha * x[2])
    res1 <- apply(theta, 2, function(x) (1 - f(x))/(1 - x[2]))
    res <- min(res1, 1)
    r_theta <- theta[2, which(res == res1)]
    r_theta <- if (length(r_theta) == 0)
        NA else r_theta
    list(theta = r_theta, type2 = res)
}
