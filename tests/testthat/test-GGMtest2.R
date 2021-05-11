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
                     method = "partialling out",
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

test_cases <- expand.grid(
  method = c("robust", "partialling out"),
  DML_method = c("DML1", "DML2"),
  nuisance_estimation = c("lasso", "post-lasso", "sqrt-lasso"),
  stringsAsFactors = FALSE
)

test_cases["test_name"] = apply(test_cases, 1, paste, collapse = "_")

patrick::with_parameters_test_that("Test different methods work similar over k-fold",
  .cases = test_cases,
  {
    set.seed(42)
    n <- 2000; p <- 10
    Sigma <- diag(p)
    c <- 0.3
    for (i in 2:(p-1)){
      for (j in (i+1):p){
        Sigma[i,j] <- Sigma[j,i] <- c^(j-i)
      }
    }
    X <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma)

    edges <- matrix(c(rep(1,p-1),2:p),ncol = 2)
    model1 <- GGMtest2(data = X,
                       edges = edges,
                       method = method,
                       nuisance_estimaton = nuisance_estimation,
                       DML_method = DML_method,
                       k_fold = 1)
    model2 <- GGMtest2(data = X,
                       edges = edges,
                       method = method,
                       nuisance_estimaton = nuisance_estimation,
                       DML_method = DML_method,
                       k_fold = 5)
    expect_equal(model1$estimates, model2$estimates,tolerance = 1e-1)
    expect_false(all(model1$estimates == model2$estimates))
  }
)

test_cases <- expand.grid(
  method = c("robust", "partialling out"),
  DML_method = c("DML1", "DML2"),
  nuisance_estimation = c("lasso", "post-lasso", "sqrt-lasso"),
  k_fold = c(1,5),
  stringsAsFactors = FALSE
)

test_cases["test_name"] = apply(test_cases, 1, paste, collapse = "_")

patrick::with_parameters_test_that("Hypothesis testing",
  .cases = test_cases,
  {
    set.seed(42)
    n <- 500; p <- 10
    Sigma <- diag(p)
    c <- 0.3
    for (i in 2:(p-1)){
      for (j in (i+1):p){
        Sigma[i,j] <- Sigma[j,i] <- c^(j-i)
      }
    }
    X <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma)

    edges <- matrix(c(rep(1,p-1),2:p),ncol = 2)
    model <- GGMtest2(data = X,
                      edges = edges,
                      method = method,
                      nuisance_estimaton = nuisance_estimation,
                      DML_method = DML_method,
                      k_fold = k_fold)

    CR1 <- create_CR(model = model,
                     alpha = 0.05,
                     B = 500,
                     s = 1,
                     exp = 1)
    CR2 <- create_CR(model = model,
                     alpha = 0.05,
                     B = 500,
                     s = 5,
                     exp = 1)
    CR3 <- create_CR(model = model,
                     alpha = 0.05,
                     B = 500,
                     s = 5,
                     exp = 2)
    CR4 <- create_CR(model = model,
                     alpha = 0.05,
                     B = 500,
                     s = 10,
                     exp = 2)

    expect_true(CR1$hyp_max)
    expect_true(CR2$hyp_max)
    expect_true(CR3$hyp_max)
    expect_true(CR4$hyp_max)
    }
)

