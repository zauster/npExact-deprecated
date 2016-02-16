
source("/home/reitero/Arbeit/Rprogramming/npExact/npExact/R/npMeanSingle.R")
context("npMeanSingle")

ones <- rep(1, 20)
zeros <- rep(0, 20)

res <- npMeanSingle(ones, mu = .5)
theta.1 <- res$theta
test_that("correct rejection 1",
          expect_true(res$rejection))

res <- npMeanSingle(zeros, mu = .5)
test_that("correct rejection 1",
          expect_true(res$rejection))

res <- npMeanSingle(ones, mu = .5, alternative = "greater", alpha = 0.025)
test_that("theta equal in two.sided and greater alternative with alpha / 2",
          expect_equal(theta.1, res$theta))

res <- npMeanSingle(ones, mu = .5, alternative = "greater")
test_that("theta unequal in two.sided and greater alternative",
          expect_true(theta.1 != res$theta))


dta <- c(rep(1, 20), rep(.5, 13))
test_that("correct rejection 2",
          expect_true(npMeanSingle(dta, mu = 0.6)$rejection))

test_that("correct not rejection 1",
          expect_false(npMeanSingle(dta, mu = 0.8)$rejection))

lower <- 0
upper <- 1
x <- c(rep(0, 4), rep(.5, 13), rep(1, 20))
x <- (x - lower)/(upper - lower)
n <- length(x)
test_that("transBinomtest, mixed rejection",
          expect_equal({
              set.seed(1); p <- .6; transBinomTest(p, n, 0.03, list(x = x, xp = x - p))
            }, 0.24926764))
test_that("transBinomtest, full rejection",
          expect_equal({
              set.seed(1); p <- .5; transBinomTest(p, n, 0.03, list(x = x, xp = x - p))
            }, 1))
test_that("transBinomtest, no rejection",
          expect_equal({
              set.seed(1); p <- .7; transBinomTest(p, n, 0.03, list(x = x, xp = x - p))
            }, 0))

x <- 1 - c(rep(0, 4), rep(.5, 13), rep(1, 20))
x <- (x - lower)/(upper - lower)
n <- length(x)

test_that("transBinomtest, mixed rejection",
          expect_equal({
              set.seed(1); p <- .175; transBinomTest(p, n, 0.03, list(x = x, xp = x - p))
            }, 0.26506823))

test_that("transBinomtest, full rejection",
          expect_equal({
              set.seed(1); p <- .1; transBinomTest(p, n, 0.03, list(x = x, xp = x - p))
            }, 1))

test_that("transBinomtest, no rejection",
          expect_equal({
              set.seed(1); p <- .2; transBinomTest(p, n, 0.03, list(x = x, xp = x - p))
          }, 0))

set.seed(123)
x <- runif(2)
res <- npMeanSingle(x, mu = .3)
test_that("npMeanSingle, no theta calculation. two-sided",
          expect_true(is.null(res$theta)))

res <- npMeanSingle(x, mu = .3, alternative = "greater")
test_that("npMeanSingle, no theta calculation. greater",
          expect_true(is.null(res$theta)))

res <- npMeanSingle(x, mu = .3, alternative = "less")
test_that("npMeanSingle, no theta calculation. less",
          expect_true(is.null(res$theta)))
