

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
