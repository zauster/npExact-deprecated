
source("/home/reitero/Arbeit/Rprogramming/npExact/npExact/R/npStochinUnpaired.R")
context("npStochinUnpaired")

ones <- rep(1, 20)
zeros <- rep(0, 20)

res <- npStochinUnpaired(ones, zeros)
theta.1 <- res$theta
test_that("correct rejection 1",
          expect_true(res$rejection))

res <- npStochinUnpaired(zeros, ones)
test_that("correct rejection 2",
          expect_true(res$rejection))

res <- npStochinUnpaired(zeros, ones, alternative = "less", alpha = 0.025)
test_that("theta equal in both two.sided and greater (with alpha / 2)",
          expect_equal(theta.1, res$theta))

res <- npStochinUnpaired(zeros, ones, alternative = "less")
test_that("theta unequal in both two.sided and greater",
          expect_true(theta.1 != res$theta))


##
## some more tests
## 
set.seed(123)
high <- rnorm(50) + 2
low <- rnorm(50)

##
## high and low
##

## should reject: high > low?
res <- npStochinUnpaired(high, low, d = 0.6,
                          alternative = "greater")
test_that("correct rejection 1",
          expect_true(res$rejection))
test_that("sample estimate 1",
          expect_true(res$stochin.estimate > 0.8))

## should reject: high > low?
res <- npStochinUnpaired(high, low, d = -0.6,
                         alternative = "greater")
test_that("correct rejection 2",
          expect_true(res$rejection))
test_that("sample estimate 2",
          expect_true(res$stochin.estimate > 0.8))

## should NOT reject: high < low?
res <- npStochinUnpaired(high, low, d = 0.6,
                         alternative = "less")
test_that("correct non-rejection 3",
          expect_true(!res$rejection))
test_that("sample estimate 3",
          expect_true(res$stochin.estimate > 0.8))

## should NOT reject: high < low?
res <- npStochinUnpaired(high, low, d = -0.6,
                         alternative = "less")
test_that("correct non-rejection 4",
          expect_true(!res$rejection))
test_that("sample estimate 4",
          expect_true(res$stochin.estimate > 0.8))

##
## low and high
##

## should NOT reject: low > high?
res <- npStochinUnpaired(low, high, d = -0.6,
                         alternative = "greater")
test_that("correct non-rejection 5",
          expect_true(!res$rejection))
test_that("sample estimate 5",
          expect_true(res$stochin.estimate < -0.8))

## should NOT reject: low > high?
res <- npStochinUnpaired(low, high, d = 0.6,
                         alternative = "greater")
test_that("correct non-rejection 6",
          expect_true(!res$rejection))
test_that("sample estimate 6",
          expect_true(res$stochin.estimate < -0.8))

## should reject: low < high?
res <- npStochinUnpaired(low, high, d = -0.6,
                         alternative = "less")
test_that("correct rejection 7",
          expect_true(res$rejection))
test_that("sample estimate 7",
          expect_true(res$stochin.estimate < -0.8))

## should reject: low < high?
res <- npStochinUnpaired(low, high, d = 0.6,
                         alternative = "less")
test_that("correct rejection 8",
          expect_true(res$rejection))
test_that("sample estimate 8",
          expect_true(res$stochin.estimate < -0.8))

