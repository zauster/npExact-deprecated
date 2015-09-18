
source("/home/reitero/Arbeit/Rprogramming/npExact/npExact/R/npMeanSingle.R")
context("npMeanSingle")

test_that("correct rejection 1",
          expect_true(npMeanSingle(rep(1, 20), mu = 0.5)$rejection))

dta <- c(rep(1, 20), rep(.5, 13))

test_that("correct rejection 2",
          expect_true(npMeanSingle(dta, mu = 0.6)$rejection))

test_that("correct not rejection 1",
          expect_false(npMeanSingle(dta, mu = 0.8)$rejection))

## test_that("more iterations",
##           expect_true(npMeanSingle(dta, mu = 0.675)$iterations > 5000))


lower <- 0
upper <- 1
x <- c(rep(0, 4), rep(.5, 13), rep(1, 20))
x <- (x - lower)/(upper - lower)
n <- length(x)
test_that("transBinomtest, mixed rejection",
          expect_equal({
              set.seed(1); p <- .6; transBinomTest(x, p, x - p, n, 0.03)
            }, 0.24926764))
test_that("transBinomtest, full rejection",
          expect_equal({
              set.seed(1); p <- .5; transBinomTest(x, p, x - p, n, 0.03)
            }, 1))
test_that("transBinomtest, no rejection",
          expect_equal({
              set.seed(1); p <- .7; transBinomTest(x, p, x - p, n, 0.03)
            }, 0))

x <- 1 - c(rep(0, 4), rep(.5, 13), rep(1, 20))
x <- (x - lower)/(upper - lower)
n <- length(x)
test_that("transBinomtest, mixed rejection",
          expect_equal({
              set.seed(1); p <- .175; transBinomTest(x, p, x - p, n, 0.03)
            }, 0.26506823))
test_that("transBinomtest, full rejection",
          expect_equal({
              set.seed(1); p <- .1; transBinomTest(x, p, x - p, n, 0.03)
            }, 1))
test_that("transBinomtest, no rejection",
          expect_equal({
              set.seed(1); p <- .2; transBinomTest(x, p, x - p, n, 0.03)
            }, 0))
