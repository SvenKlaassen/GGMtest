context("GGMtest2 unit tests")

test_that("Compare GGMtest and GGMtest2", {
  set.seed(42)
  n <- 200; p <- 10
  X <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=diag(p))

  edges <- matrix(c(rep(1,p-1),2:p),ncol = 2)
  model1 <- GGMtest(data = X,
                   edges = edges,
                   nuisance_estimaton = "sqrt-lasso",
                   s = 1,
                   exponent = 1)

  model2 <- GGMtest2(data = X,
                     edges = edges,
                     nuisance_estimaton = "sqrt-lasso")

  CR <- create_CR(model = model2,
                  alpha = 0.05,
                  B = 500,
                  s = 1,
                  exp = 1)

  expect_equal(model1$estimates,model2$estimates)
  expect_true(as.logical(model1$hyp_max*CR$hyp_max))
  expect_true(as.logical(model1$hyp_sphere*CR$hyp_sphere))
})

test_that("Test volumes", {
  set.seed(42)
  n <- 200; p <- 100
  X <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=diag(p))

  edges <- matrix(c(rep(1,p-1),2:p),ncol = 2)
  model <- GGMtest2(data = X,
                     edges = edges,
                     nuisance_estimaton = "sqrt-lasso")

  CR <- create_CR(model = model,
                  alpha = 0.05,
                  B = 500,
                  s = 1,
                  exp = 1)

  vol_max <- prod(2*model$sigma_est*CR$quantiles[[1]]/sqrt(n))
  vol_sphere <- prod(2*(model$sigma_est*CR$quantiles[[2]]/sqrt(n)-model$sigma_est*CR$quantiles[[3]]/sqrt(n)))

  expect_equal(CR$vol_max,vol_max)
  expect_equal(CR$vol_sphere,vol_sphere)
  })
