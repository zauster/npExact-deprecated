

source("/home/reitero/Arbeit/Rprogramming/npExact/npExact/R/npMeanUnpaired.R")
context("npMeanUnpaired")

ones <- rep(1, 20)
zeros <- rep(0, 20)

res <- npMeanUnpaired(ones, zeros)
theta.1 <- res$theta
test_that("correct rejection 1",
          expect_true(res$rejection))

res <- npMeanUnpaired(zeros, ones)
test_that("correct rejection 2",
          expect_true(res$rejection))

res <- npMeanUnpaired(zeros, ones, alternative = "less", alpha = 0.025)
test_that("theta equal in both two.sided and greater",
          expect_equal(theta.1, res$theta))
