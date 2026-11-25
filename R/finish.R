## Parameter packing, standard errors, and output assembly.

par_layout <- function(st) {
  J <- length(st$a)
  lab <- character(0); blk <- character(0)
  for (j in seq_len(J)) {
    lab <- c(lab, paste0("a[", j, "]"), paste0("d[", j, ",", seq_len(st$K[j] - 1L), "]"))
    blk <- c(blk, rep("item", st$K[j]))
  }
  if (st$mixture) {
    lab <- c(lab, paste0("gamma:", colnames(st$Z))); blk <- c(blk, rep("gamma", ncol(st$Z)))
    if (st$pi_type == "free") {
      lab <- c(lab, paste0("logit_pi[", which(st$neg), "]")); blk <- c(blk, rep("pi", sum(st$neg)))
    } else if (st$pi_type == "features") {
      lab <- c(lab, paste0("delta:", colnames(st$Wf))); blk <- c(blk, rep("delta", ncol(st$Wf)))
    }
    if (st$theta_class %in% c("mean", "meanvar")) { lab <- c(lab, "mu1"); blk <- c(blk, "theta") }
    if (st$theta_class == "meanvar") { lab <- c(lab, "log_sigma1"); blk <- c(blk, "theta") }
  }
  list(label = lab, block = blk)
}

pack_par <- function(st) {
  p <- numeric(0)
  for (j in seq_along(st$a)) p <- c(p, st$a[j], raw_from_d(st$d[[j]]))
  if (st$mixture) {
    p <- c(p, st$gamma)
    if (st$pi_type == "free") p <- c(p, stats::qlogis(st$pi[st$neg]))
    if (st$pi_type == "features") p <- c(p, st$delta)
    if (st$theta_class %in% c("mean", "meanvar")) p <- c(p, st$mu1)
    if (st$theta_class == "meanvar") p <- c(p, log(st$sigma1))
  }
  p
}

unpack_par <- function(st, p) {
  k <- 0L
  for (j in seq_along(st$a)) {
    st$a[j] <- p[k + 1L]
    st$d[[j]] <- d_from_raw(p[k + 1L + seq_len(st$K[j] - 1L)])
    k <- k + st$K[j]
  }
  if (st$mixture) {
    q <- length(st$gamma); st$gamma <- p[k + seq_len(q)]; k <- k + q
    if (st$pi_type == "free") {
      q <- sum(st$neg); st$pi[st$neg] <- stats::plogis(p[k + seq_len(q)]); k <- k + q
    }
    if (st$pi_type == "features") {
      q <- length(st$delta); st$delta <- p[k + seq_len(q)]; k <- k + q
      st$pi[st$neg] <- stats::plogis(drop(st$Wf %*% st$delta))
    }
    if (st$theta_class %in% c("mean", "meanvar")) { st$mu1 <- p[k + 1L]; k <- k + 1L }
    if (st$theta_class == "meanvar") { st$sigma1 <- exp(p[k + 1L]); k <- k + 1L }
    st$g1 <- class_prior(st$theta, st$mu1, st$sigma1)
  }
  st
}

# Per-respondent scores by central differences, then (cluster-robust) covariance.
score_vcov <- function(st, clus = NULL, free = NULL) {
  p <- pack_par(st)
  if (is.null(free)) free <- rep(TRUE, length(p))
  idx <- which(free)
  S <- matrix(0, nrow(st$X), length(idx))
  for (k in seq_along(idx)) {
    h <- 1e-4 * max(1, abs(p[idx[k]]))
    pp <- p; pm <- p
    pp[idx[k]] <- p[idx[k]] + h; pm[idx[k]] <- p[idx[k]] - h
    S[, k] <- (person_loglik(unpack_par(st, pp)) - person_loglik(unpack_par(st, pm))) / (2 * h)
  }
  S <- S * st$w
  A <- crossprod(S)
  Ainv <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  V <- if (is.null(clus)) Ainv else {
    Sc <- rowsum(S, clus)
    Ainv %*% crossprod(Sc) %*% Ainv
  }
  out <- matrix(NA_real_, length(p), length(p))
  out[idx, idx] <- V
  out
}

finish_fit <- function(st, items, reversed, n_dropped, clus, se, call) {
  lay <- par_layout(st)
  p <- pack_par(st)
  free <- rep(TRUE, length(p))
  if (st$mixture && st$pi_type == "free") {
    pi_pos <- which(lay$block == "pi")
    bnd <- st$pi[st$neg] < 1e-4 | st$pi[st$neg] > 1 - 1e-4
    free[pi_pos[bnd]] <- FALSE
  }
  V <- if (se) score_vcov(st, clus, free) else matrix(NA_real_, length(p), length(p))
  sev <- sqrt(pmax(diag(V), 0))
  names(p) <- names(sev) <- lay$label
  dimnames(V) <- list(lay$label, lay$label)
  J <- length(items)
  Kmax <- max(st$K)
  dmat <- t(vapply(st$d, function(d) c(d, rep(NA_real_, Kmax - 1L - length(d))), numeric(Kmax - 1L)))
  if (Kmax == 2L) dmat <- matrix(dmat, ncol = 1L)
  colnames(dmat) <- paste0("d", seq_len(Kmax - 1L))
  a_se <- sev[paste0("a[", seq_len(J), "]")]
  item_tab <- data.frame(item = items, reversed = st$neg, a = st$a, a_se = unname(a_se), dmat,
                         row.names = NULL, check.names = FALSE)
  res <- list(call = call, items = item_tab, d = st$d, par = p, se = sev, vcov = V,
              loglik = st$loglik, iterations = st$iterations, converged = st$converged,
              n = nrow(st$X), n_dropped = n_dropped, pi_type = st$pi_type,
              theta_class = st$theta_class, reversed = reversed, item_names = items)
  npar <- length(p) - sum(!free)
  if (!st$mixture) npar <- length(p)
  res$npar <- npar
  res$AIC <- -2 * st$loglik + 2 * npar
  res$BIC <- -2 * st$loglik + log(nrow(st$X)) * npar
  if (st$mixture) {
    gpos <- which(lay$block == "gamma")
    res$gamma <- st$gamma
    res$gamma_table <- coef_table(st$gamma, sev[gpos], colnames(st$Z))
    pis <- st$pi[st$neg]
    if (st$pi_type == "free") {
      lse <- sev[which(lay$block == "pi")]
      pse <- pis * (1 - pis) * lse
    } else if (st$pi_type == "features") {
      dpos <- which(lay$block == "delta")
      G <- st$Wf * (pis * (1 - pis))
      pse <- sqrt(pmax(rowSums((G %*% V[dpos, dpos, drop = FALSE]) * G), 0))
      res$delta <- st$delta
      res$delta_table <- coef_table(st$delta, sev[dpos], colnames(st$Wf))
    } else pse <- rep(NA_real_, length(pis))
    res$pi <- data.frame(item = reversed, pi = pis, se = unname(pse),
                         boundary = pis < 1e-4 | pis > 1 - 1e-4, row.names = NULL)
    res$pi_full <- st$pi
    res$mu1 <- st$mu1; res$sigma1 <- st$sigma1
    res$mu1_se <- if ("mu1" %in% names(sev)) unname(sev["mu1"]) else NA_real_
    post1 <- rowSums(st$post1)
    res$posterior <- post1
    res$prevalence <- sum(st$w * post1) / sum(st$w)
    pp <- cbind(1 - post1, post1)
    ent <- -sum(st$w * rowSums(pp * log(pmax(pp, 1e-300))))
    res$entropy <- 1 - ent / (sum(st$w) * log(2))
    res$classification_accuracy <- sum(st$w * pmax(post1, 1 - post1)) / sum(st$w)
  }
  postall <- st$post0 + st$post1
  res$eap <- drop(postall %*% st$theta)
  res$weights <- st$w
  res$state <- st[c("theta", "K", "neg", "Z", "Wf", "mixture", "pi_type", "theta_class")]
  class(res) <- "revmix"
  res
}

coef_table <- function(est, se, nm) {
  z <- est / se
  data.frame(term = nm, estimate = unname(est), se = unname(se), z = unname(z),
             p = 2 * stats::pnorm(-abs(unname(z))), row.names = NULL)
}

#' @export
print.revmix <- function(x, digits = 3, ...) {
  lab <- switch(x$pi_type, none = "single-class GRM", all = "mixture GRM, all reverse-worded items misread (pi = 1)",
                free = "mixture GRM, item-specific misreading", features = "mixture GRM, misreading explained by item features")
  cat("revmix:", lab, "\n")
  cat(sprintf("N = %d, items = %d (%d reverse-worded); logLik = %.2f, npar = %d, AIC = %.1f, BIC = %.1f%s\n",
              x$n, nrow(x$items), sum(x$items$reversed), x$loglik, x$npar, x$AIC, x$BIC,
              if (x$converged) "" else " [EM did not converge]"))
  if (x$pi_type != "none") {
    cat(sprintf("Inconsistent class: prevalence = %.3f, trait mean = %.3f, entropy = %.3f\n",
                x$prevalence, x$mu1, x$entropy))
    cat("\nClass membership (logit, inconsistent vs consistent):\n")
    print(round_df(x$gamma_table, digits), row.names = FALSE)
    cat("\nMisreading probabilities of reverse-worded items:\n")
    print(round_df(x$pi, digits), row.names = FALSE)
    if (!is.null(x$delta_table)) {
      cat("\nMisreading model (logit):\n")
      print(round_df(x$delta_table, digits), row.names = FALSE)
    }
  }
  invisible(x)
}

round_df <- function(d, digits) {
  num <- vapply(d, is.numeric, TRUE)
  d[num] <- lapply(d[num], round, digits)
  d
}

#' @export
logLik.revmix <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$n, class = "logLik")
}

#' @export
coef.revmix <- function(object, ...) object$par

#' @export
vcov.revmix <- function(object, ...) object$vcov
