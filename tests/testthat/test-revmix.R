test_that("GRM probabilities are valid and the item M-step gradient is correct", {
  th <- seq(-3, 3, length.out = 11)
  P <- grm_prob(1.3, c(1, 0, -1.2), th)
  expect_equal(rowSums(P), rep(1, 11), tolerance = 1e-10)
  expect_true(all(P > 0))
  expect_equal(d_from_raw(raw_from_d(c(1.5, 0.2, -0.7))), c(1.5, 0.2, -0.7), tolerance = 1e-10)
  set.seed(2)
  r <- matrix(stats::rpois(11 * 4, 20), 11, 4)
  f <- function(p) -sum(r * log(grm_prob(p[1], d_from_raw(p[-1]), th)))
  p0 <- c(1.1, raw_from_d(c(1, 0, -1)))
  num <- vapply(seq_along(p0), function(k) {
    h <- 1e-6; e <- replace(numeric(4), k, h); (f(p0 + e) - f(p0 - e)) / (2 * h)
  }, 0)
  # reproduce analytic gradient used inside mstep_item
  K <- 4; a <- p0[1]; raw <- p0[-1]; d <- d_from_raw(raw)
  z <- a * th + matrix(d, 11, 3, byrow = TRUE); Ps <- stats::plogis(z); W <- Ps * (1 - Ps)
  Pm <- grm_prob(a, d, th); R <- r / Pm; Wf <- cbind(0, W, 0)
  ga <- sum(R * th * (Wf[, 1:4] - Wf[, 2:5])); gd <- colSums(W * (R[, -1] - R[, -4]))
  J <- matrix(0, 3, 3); J[, 1] <- 1; for (m in 2:3) J[m:3, m] <- -exp(raw[m])
  expect_equal(-c(ga, drop(gd %*% J)), num, tolerance = 1e-4)
})

test_that("input checks work", {
  sim <- rm_simulate(n = 50, n_items = 6, n_reversed = 2)
  expect_error(rm_fit(sim$data, sim$items, "zz"), "subset")
  expect_error(rm_fit(sim$data, c(sim$items, "zz"), sim$reversed), "not in")
  expect_error(rm_fit(sim$data, sim$items, sim$reversed, pi = "features"), "item_features")
  expect_error(rm_simulate(n_items = 4, n_reversed = 4), "n_reversed")
})

test_that("simulated data are keyed integers with the documented structure", {
  set.seed(3)
  sim <- rm_simulate(n = 300, n_items = 8, n_reversed = 3, n_cat = 5)
  expect_equal(dim(sim$data), c(300, 9))
  expect_true(all(unlist(sim$data[sim$items]) %in% 0:4))
  expect_length(sim$reversed, 3)
})

test_that("the single-class model and the misreading mixture are estimated", {
  set.seed(4)
  sim <- rm_simulate(n = 600, n_items = 8, n_reversed = 4, prevalence = 0.25, pi = 1)
  f0 <- rm_fit(sim$data, sim$items, sim$reversed, pi = "none", se = FALSE, maxit = 300)
  fa <- rm_fit(sim$data, sim$items, sim$reversed, pi = "all", se = FALSE, nstart = 1, maxit = 400)
  expect_s3_class(fa, "revmix")
  expect_gt(fa$loglik, f0$loglik)
  expect_lt(fa$BIC, f0$BIC)
  expect_true(all(fa$pi$pi == 1))
  expect_equal(nrow(rm_scores(fa)), 600)
  cmp <- rm_compare(f0, fa)
  expect_equal(nrow(cmp), 2)
})

test_that("item-specific misreading and the covariate effect are recovered", {
  skip_on_cran()
  set.seed(11)
  sim <- rm_simulate(n = 3000, n_items = 12, n_reversed = 6, gamma = c(-1.5, -0.8),
                     pi = c(0.95, 0.85, 0.7, 0.6, 0.4, 0.25))
  ff <- rm_fit(sim$data, sim$items, sim$reversed, covariates = ~z, nstart = 1, nquad = 41)
  expect_equal(ff$gamma_table$estimate[2], -0.8, tolerance = 0.2)
  expect_lt(max(abs(ff$pi$pi - sim$pi)), 0.15)
  expect_true(all(is.finite(ff$gamma_table$se)))
  expect_gt(ff$pi$pi[1], ff$pi$pi[6])
})

test_that("feature model, weights, clusters and plausible values run", {
  skip_on_cran()
  set.seed(5)
  sim <- rm_simulate(n = 800, n_items = 8, n_reversed = 4, gamma = c(-1.2, -0.6),
                     pi = c(0.9, 0.9, 0.4, 0.4))
  d <- sim$data
  d$w <- stats::runif(800, 0.5, 2); d$sch <- rep(1:40, each = 20)
  d$pv1 <- d$z * 100 + 500 + stats::rnorm(800, 0, 30); d$pv2 <- d$z * 100 + 500 + stats::rnorm(800, 0, 30)
  feat <- data.frame(hard = c(0, 0, 1, 1))
  fx <- rm_fit(d, sim$items, sim$reversed, pi = "features", item_features = feat, features = ~hard,
               covariates = ~z, weights = "w", cluster = "sch", nstart = 1, nquad = 41)
  expect_equal(nrow(fx$delta_table), 2)
  expect_true(all(is.finite(fx$pi$se)))
  pv <- rm_pv(d, sim$items, sim$reversed, pv = c("pv1", "pv2"), se = TRUE, nstart = 1, nquad = 41)
  expect_s3_class(pv, "revmix_pv")
  expect_equal(nrow(pv$gamma), 2)
})

test_that("bundled data are keyed and complete", {
  expect_equal(nrow(moderation), 4197)
  items <- names(moderation)[grepl("^[A-D][0-9]+$", names(moderation))]
  expect_length(items, 36)
  expect_false(anyNA(moderation[items]))
})
