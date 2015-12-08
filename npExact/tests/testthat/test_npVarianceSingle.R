
source("/home/reitero/Arbeit/Rprogramming/npExact/npExact/R/npVarianceSingle.R")
context("npVarianceSingle")

set.seed(123)
x <- runif(20)


res <- npVarianceSingle(x, v = .2)
theta.1 <- res$theta
test_that("correct not-rejection 1",
          expect_false(res$rejection))

res <- npVarianceSingle(x, v = .2, alternative = "less", alpha = 0.025)
test_that("theta unequal in two.sided and greater alternative",
          expect_equal(theta.1, res$theta))

res <- npVarianceSingle(x, v = .2, alternative = "greater", alpha = 0.025)
test_that("theta equal in two.sided and greater alternative with alpha / 2",
          expect_true(theta.1 != res$theta))


res <- npVarianceSingle(x, v = .25)
theta.1 <- res$theta
test_that("correct rejection 1",
          expect_true(res$rejection))

res <- npVarianceSingle(x, v = .25, alternative = "less", alpha = 0.025)
test_that("theta unequal in two.sided and greater alternative",
          expect_equal(theta.1, res$theta))


##
## theta calculation
##

set.seed(123)
x <- runif(2)
res <- npVarianceSingle(x, v = 0.01)
test_that("npVarianceSingle, no theta calculation. two-sided",
          expect_true(is.null(res$theta)))

res <- npVarianceSingle(x, v = 0.1, alternative = "greater")
test_that("npVarianceSingle, no theta calculation. greater",
          expect_true(is.null(res$theta)))

res <- npVarianceSingle(x, v = 0.01, alternative = "less")
test_that("npVarianceSingle, no theta calculation. less",
          expect_true(is.null(res$theta)))
