context("GGMtest unit tests")


test_that("GGMtest", {
  set.seed(42)
  n <- 200; p <- 10
  X <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=diag(p))

  edges <- matrix(c(rep(1,p-1),2:p),ncol = 2)
  model <- GGMtest(data = X, edges = edges, nuisance_estimaton = "sqrt-lasso")

  expect_true(model$hyp_max)
  expect_true(model$hyp_sphere)
  expect_length(model$estimates,dim(edges)[1])
  expect_true(model$pvalue_max > model$alpha)
})
