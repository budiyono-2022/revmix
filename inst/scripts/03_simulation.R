# Simulation study for revmix: parameter recovery, interval coverage,
# Type I error and power of the covariate effect, classification, and the
# consequences of forcing all misreading probabilities to 1.
# Run from the package root:  Rscript inst/scripts/03_simulation.R
suppressMessages(library(revmix))
library(parallel)

cond <- rbind(
  expand.grid(n = c(1000, 3000), J = c(8, 16), pattern = c("homogeneous", "heterogeneous"),
              gamma1 = -0.8, stringsAsFactors = FALSE),
  expand.grid(n = c(1000, 3000), J = 8, pattern = "heterogeneous", gamma1 = 0,
              stringsAsFactors = FALSE))
reps <- as.integer(Sys.getenv("REVMIX_REPS", "200"))
cores <- as.integer(Sys.getenv("REVMIX_CORES", "2"))
out_file <- Sys.getenv("REVMIX_OUT", "sim_results.rds")

pi_vec <- function(nrev, pattern) {
  if (pattern == "homogeneous") rep(0.9, nrev) else seq(0.95, 0.25, length.out = nrev)
}

one_rep <- function(r, cc) {
  p <- cond[cc, ]
  set.seed(1e5 * cc + r)
  nrev <- p$J / 2
  pis <- pi_vec(nrev, p$pattern)
  g <- c(-1.9, p$gamma1)                      # about 15% prevalence at z = 0
  sim <- rm_simulate(n = p$n, n_items = p$J, n_reversed = nrev, gamma = g, pi = pis, mu1 = 0)
  res <- list(ok = FALSE)
  ff <- tryCatch(rm_fit(sim$data, sim$items, sim$reversed, pi = "free", covariates = ~z,
                        nquad = 41, tol = 1e-5, nstart = 1, se = TRUE), error = function(e) NULL)
  fa <- tryCatch(rm_fit(sim$data, sim$items, sim$reversed, pi = "all", covariates = ~z,
                        nquad = 41, tol = 1e-5, nstart = 1, se = FALSE), error = function(e) NULL)
  f0 <- tryCatch(rm_fit(sim$data, sim$items, sim$reversed, pi = "none", nquad = 41, tol = 1e-5,
                        se = FALSE), error = function(e) NULL)
  if (is.null(ff) || is.null(fa) || is.null(f0)) return(res)
  gt <- ff$gamma_table
  ga <- fa$gamma_table
  cls_hat <- as.integer(ff$posterior > 0.5)
  list(ok = TRUE, converged = ff$converged,
       g0 = gt$estimate[1], g1 = gt$estimate[2], g0_se = gt$se[1], g1_se = gt$se[2],
       g1_all = ga$estimate[2], g0_all = ga$estimate[1],
       pi = ff$pi$pi, pi_se = ff$pi$se, pi_true = pis,
       a_err = ff$items$a - sim$a,
       prev = ff$prevalence, prev_all = fa$prevalence, prev_true = mean(sim$class),
       acc = mean(cls_hat == sim$class), entropy = ff$entropy,
       sens = mean(cls_hat[sim$class == 1] == 1), spec = mean(cls_hat[sim$class == 0] == 0),
       bic = c(none = f0$BIC, all = fa$BIC, free = ff$BIC),
       eap_rmse = sqrt(mean((ff$eap - sim$theta)^2)),
       eap0_rmse = sqrt(mean((f0$eap - sim$theta)^2)))
}

results <- if (file.exists(out_file)) readRDS(out_file) else list(cond = cond, res = list())
for (cc in seq_len(nrow(cond))) {
  if (length(results$res) >= cc && !is.null(results$res[[cc]])) next
  results$res[[cc]] <- mclapply(seq_len(reps), one_rep, cc = cc, mc.cores = cores)
  saveRDS(results, out_file)
  cat(cc, format(Sys.time()), "\n")
}
