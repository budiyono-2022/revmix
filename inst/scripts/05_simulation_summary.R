# Summarize the simulation study (output of 03_simulation.R).
# Usage: Rscript inst/scripts/05_simulation_summary.R sim_results.rds [figure.png]
args <- commandArgs(trailingOnly = TRUE)
res <- readRDS(if (length(args)) args[1] else "sim_results.rds")
fig <- if (length(args) > 1) args[2] else "simulation_figure.png"
cond <- res$cond

summ <- lapply(seq_along(res$res), function(cc) {
  rr <- res$res[[cc]]
  if (is.null(rr)) return(NULL)
  ok <- vapply(rr, function(x) isTRUE(x$ok), TRUE)
  r <- rr[ok]
  p <- cond[cc, ]
  get <- function(k) vapply(r, function(x) x[[k]][1], 0)
  g1 <- get("g1"); s1 <- get("g1_se"); g0 <- get("g0")
  fin <- is.finite(s1)
  pis <- do.call(rbind, lapply(r, `[[`, "pi")); pse <- do.call(rbind, lapply(r, `[[`, "pi_se"))
  pit <- r[[1]]$pi_true
  piok <- is.finite(pse)
  cov_pi <- mean((abs(pis - matrix(pit, nrow(pis), ncol(pis), byrow = TRUE)) <= 1.96 * pse)[piok])
  bic <- do.call(rbind, lapply(r, `[[`, "bic"))
  sel <- table(factor(colnames(bic)[apply(bic, 1, which.min)], levels = c("none", "all", "free")))
  prevt <- get("prev_true")
  data.frame(
    cond = cc, n = p$n, J = p$J, pattern = p$pattern, gamma1 = p$gamma1,
    reps = length(rr), ok = sum(ok), converged = mean(vapply(r, function(x) isTRUE(x$converged), TRUE)),
    g0_bias = mean(g0) - (-1.9), g1_bias = mean(g1) - p$gamma1,
    g1_rmse = sqrt(mean((g1 - p$gamma1)^2)),
    g1_sd = sd(g1), g1_mse = mean(s1[fin]),
    g1_cover = mean(abs(g1 - p$gamma1)[fin] <= 1.96 * s1[fin]),
    g1_reject = mean(abs(g1 / s1)[fin] > 1.96),
    g1_all_bias = mean(get("g1_all")) - p$gamma1,
    g0_all_bias = mean(get("g0_all")) - (-1.9),
    pi_bias = mean(colMeans(pis) - pit), pi_mabias = mean(abs(colMeans(pis) - pit)),
    pi_rmse = sqrt(mean((pis - matrix(pit, nrow(pis), ncol(pis), byrow = TRUE))^2)),
    pi_cover = cov_pi,
    prev_true = mean(prevt), prev_bias = mean(get("prev") - prevt),
    prev_all_bias = mean(get("prev_all") - prevt),
    acc = mean(get("acc")), sens = mean(get("sens")), spec = mean(get("spec")),
    entropy = mean(get("entropy")),
    sel_none = sel[["none"]] / sum(sel), sel_all = sel[["all"]] / sum(sel), sel_free = sel[["free"]] / sum(sel),
    eap_rmse = mean(get("eap_rmse")), eap0_rmse = mean(get("eap0_rmse")),
    a_bias = mean(unlist(lapply(r, `[[`, "a_err"))),
    stringsAsFactors = FALSE)
})
tab <- do.call(rbind, summ)
options(width = 200)
print(format(tab, digits = 3), row.names = FALSE)
saveRDS(tab, sub("\\.rds$", "_summary.rds", if (length(args)) args[1] else "sim_results.rds"))

# Item-level recovery of pi in the heterogeneous gamma1 = -0.8 conditions
het <- which(cond$pattern == "heterogeneous" & cond$gamma1 != 0 & seq_len(nrow(cond)) <= length(res$res))
het <- het[!vapply(res$res[het], is.null, TRUE)]
if (length(het)) {
  png(fig, width = 6.5, height = 3.0, units = "in", res = 300, pointsize = 9)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), mgp = c(2.4, 0.7, 0))
  plot(NA, xlim = c(0.2, 1), ylim = c(0, 1.05), xlab = "Population misreading probability",
       ylab = "Estimate", main = "(a) Item-specific misreading", las = 1)
  abline(0, 1, col = "grey60", lty = 2)
  pch <- c(1, 16, 2, 17)
  for (k in seq_along(het)) {
    cc <- het[k]; r <- Filter(function(x) isTRUE(x$ok), res$res[[cc]])
    pis <- do.call(rbind, lapply(r, `[[`, "pi")); pit <- r[[1]]$pi_true
    q <- apply(pis, 2, quantile, c(.05, .95))
    off <- (k - 2.5) * 0.008
    segments(pit + off, q[1, ], pit + off, q[2, ], col = "grey45")
    points(pit + off, colMeans(pis), pch = pch[k], cex = 0.9)
  }
  lab <- sprintf("n = %s, J = %d", format(cond$n[het], big.mark = ","), cond$J[het])
  legend("topleft", lab, pch = pch[seq_along(het)], bty = "n", cex = 0.85)
  # (b) prevalence: free vs all
  rows <- tab[tab$gamma1 != 0, ]
  lab2 <- sprintf("%s\n%s\nJ = %d", ifelse(rows$pattern == "homogeneous", "Hom.", "Het."), format(rows$n, big.mark = ","), rows$J)
  m <- rbind(rows$prev_bias, rows$prev_all_bias)
  par(mar = c(5.5, 4.8, 2, 1), mgp = c(3.4, 0.7, 0))
  bp <- barplot(m, beside = TRUE, col = c("grey25", "grey80"), las = 1, ylim = c(min(m) - 0.035, 0.02),
                ylab = "Bias of estimated prevalence", main = "(b) Prevalence of inconsistent class",
                names.arg = rep("", ncol(m)))
  text(colMeans(bp), par("usr")[3] - 0.004, lab2, xpd = TRUE, cex = 0.7, adj = c(0.5, 1))
  abline(h = 0)
  legend("bottomleft", c("Item-specific misreading", expression("All misread (" * pi * " = 1)")), fill = c("grey25", "grey80"),
         bty = "n", cex = 0.85)
  par(op); dev.off()
}
