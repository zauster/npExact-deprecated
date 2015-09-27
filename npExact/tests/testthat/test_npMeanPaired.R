

source("/home/reitero/Arbeit/Rprogramming/npExact/npExact/R/npMeanPaired.R")
context("npMeanPaired")

ones <- rep(1, 20)
zeros <- rep(0, 20)

res <- npMeanPaired(ones, zeros)
theta.1 <- res$theta
test_that("correct rejection 1",
          expect_true(res$rejection))

res <- npMeanPaired(zeros, ones)
test_that("correct rejection 2",
          expect_true(res$rejection))

res <- npMeanPaired(zeros, ones, alternative = "less", alpha = 0.025)
test_that("theta equal in both two.sided and greater",
          expect_equal(theta.1, res$theta))
