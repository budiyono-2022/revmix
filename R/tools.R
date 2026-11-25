#' Fit the model once per plausible value and pool with Rubin's rules
#'
#' When class membership is regressed on a latent proficiency reported as
#' plausible values (as in PIRLS, TIMSS, or PISA), the model is fitted once per
#' plausible value and the estimates are combined with Rubin's rules. The
#' plausible value enters the membership formula through the placeholder
#' variable \code{.pv}.
#'
#' @inheritParams rm_fit
#' @param pv Character vector of plausible-value variables.
#' @param covariates Membership formula that uses \code{.pv}, e.g.
#'   \code{~ .pv}.
#' @param transform Function applied to each plausible value before fitting;
#'   the default rescales an international 500/100 metric to units of 100
#'   points.
#' @param ... Further arguments passed to \code{\link{rm_fit}}.
#' @return An object of class \code{revmix_pv} with the individual fits and
#'   pooled tables for class membership, misreading probabilities, and the
#'   inconsistent-class trait mean.
#' @export
rm_pv <- function(data, items, reversed, pv, covariates = ~ .pv,
                  transform = function(x) (x - 500) / 100, ...) {
  if (!all(pv %in% names(data))) stop("Some `pv` variables are not in `data`.", call. = FALSE)
  fits <- vector("list", length(pv))
  prev <- NULL
  for (m in seq_along(pv)) {
    data$.pv <- transform(data[[pv[m]]])
    fits[[m]] <- rm_fit(data, items, reversed, covariates = covariates, start = prev, ...)
    prev <- fits[[m]]
  }
  pool <- function(est, se) {
    est <- as.matrix(est); se <- as.matrix(se)
    M <- nrow(est)
    qbar <- colMeans(est); W <- colMeans(se^2); B <- apply(est, 2L, stats::var)
    Tv <- W + (1 + 1 / M) * B
    r <- (1 + 1 / M) * B / pmax(W, 1e-12)
    df <- (M - 1) * (1 + 1 / pmax(r, 1e-12))^2
    se_t <- sqrt(Tv)
    t <- qbar / se_t
    data.frame(estimate = qbar, se = se_t, df = df, t = t, p = 2 * stats::pt(-abs(t), df),
               fmi = (r + 2 / (df + 3)) / (r + 1), row.names = NULL)
  }
  g_est <- t(vapply(fits, function(f) f$gamma_table$estimate, numeric(nrow(fits[[1]]$gamma_table))))
  g_se <- t(vapply(fits, function(f) f$gamma_table$se, numeric(nrow(fits[[1]]$gamma_table))))
  gamma <- cbind(term = fits[[1]]$gamma_table$term, pool(g_est, g_se))
  p_est <- t(vapply(fits, function(f) f$pi$pi, numeric(nrow(fits[[1]]$pi))))
  p_se <- t(vapply(fits, function(f) f$pi$se, numeric(nrow(fits[[1]]$pi))))
  pis <- cbind(item = fits[[1]]$pi$item, pool(p_est, p_se)[, c("estimate", "se", "df", "fmi")])
  out <- list(fits = fits, gamma = gamma, pi = pis,
              prevalence = mean(vapply(fits, `[[`, 0, "prevalence")),
              mu1 = mean(vapply(fits, `[[`, 0, "mu1")),
              loglik = vapply(fits, `[[`, 0, "loglik"), pv = pv)
  class(out) <- "revmix_pv"
  out
}

#' @export
print.revmix_pv <- function(x, digits = 3, ...) {
  cat("revmix fits across", length(x$pv), "plausible values (Rubin's rules)\n")
  cat(sprintf("Mean prevalence of the inconsistent class = %.3f; mean trait mean = %.3f\n\n",
              x$prevalence, x$mu1))
  cat("Class membership (logit):\n"); print(round_df(x$gamma, digits), row.names = FALSE)
  cat("\nMisreading probabilities:\n"); print(round_df(x$pi, digits), row.names = FALSE)
  invisible(x)
}

#' Compare fitted revmix models
#'
#' Returns log-likelihood, number of parameters, AIC, and BIC for any number of
#' fits to the same data and, for consecutive pairs, a likelihood-ratio
#' statistic. The chi-square reference distribution is valid only when the
#' restriction lies in the interior of the parameter space (for example,
#' \code{pi = "features"} against \code{pi = "free"}). For restrictions on the
#' boundary, such as \code{pi = "all"} (all \eqn{\pi_j = 1}) against
#' \code{pi = "free"}, or the single-class model against a mixture, use
#' \code{\link{rm_boot_lrt}}.
#'
#' @param ... Fitted \code{revmix} objects, from most to least restricted.
#' @return A data frame.
#' @export
rm_compare <- function(...) {
  fits <- list(...)
  nm <- vapply(fits, function(f) f$pi_type, "")
  tab <- data.frame(model = nm, logLik = vapply(fits, `[[`, 0, "loglik"),
                    npar = vapply(fits, `[[`, 0, "npar"), AIC = vapply(fits, `[[`, 0, "AIC"),
                    BIC = vapply(fits, `[[`, 0, "BIC"), LR = NA_real_, df = NA_real_, p = NA_real_)
  for (k in seq_along(fits)[-1L]) {
    tab$LR[k] <- 2 * (tab$logLik[k] - tab$logLik[k - 1L])
    tab$df[k] <- tab$npar[k] - tab$npar[k - 1L]
    if (tab$df[k] > 0) tab$p[k] <- stats::pchisq(tab$LR[k], tab$df[k], lower.tail = FALSE)
  }
  tab
}

#' Parametric bootstrap likelihood-ratio test
#'
#' Tests a restricted model against a more general one when the restriction is
#' on the boundary of the parameter space, so that the usual chi-square
#' reference does not apply (for example, all misreading probabilities equal
#' to 1, or a single class). Data are simulated from the restricted fit, both
#' models are refitted, and the observed statistic is compared with the
#' bootstrap distribution.
#'
#' @param restricted,general Fitted \code{revmix} objects.
#' @param B Number of bootstrap samples.
#' @param seed Optional seed.
#' @return A list with the observed statistic, the bootstrap statistics, and
#'   the bootstrap p value \eqn{(1 + \#\{LR^* \ge LR\}) / (B + 1)}.
#' @export
rm_boot_lrt <- function(restricted, general, B = 99, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  lr_obs <- 2 * (general$loglik - restricted$loglik)
  stz <- restricted$state
  Z <- stz$Z
  zdf <- as.data.frame(Z[, -1L, drop = FALSE])
  fz <- if (ncol(Z) > 1L) stats::as.formula(paste("~", paste(sprintf("`%s`", colnames(zdf)), collapse = " + "))) else NULL
  lr <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    sim <- rm_generate(restricted)
    dat <- cbind(sim, zdf)
    args <- list(data = dat, items = restricted$item_names, reversed = restricted$reversed,
                 covariates = fz, se = FALSE, nstart = 1, weights = NULL)
    f0 <- tryCatch(do.call(rm_fit, c(args, list(pi = restricted$pi_type, start = restricted))),
                   error = function(e) NULL)
    f1 <- tryCatch(do.call(rm_fit, c(args, list(pi = general$pi_type,
                                                start = if (general$pi_type == "none") NULL else general,
                                                item_features = general$item_features,
                                                features = general$features))),
                   error = function(e) NULL)
    if (!is.null(f0) && !is.null(f1)) lr[b] <- 2 * (f1$loglik - f0$loglik)
  }
  ok <- is.finite(lr)
  list(LR = lr_obs, boot = lr, B = sum(ok), p = (1 + sum(lr[ok] >= lr_obs)) / (sum(ok) + 1))
}

# Simulate a data set from a fitted model (same covariates, unit weights).
rm_generate <- function(fit) {
  stz <- fit$state
  n <- nrow(stz$Z); J <- nrow(fit$items)
  cls <- if (stz$mixture) stats::rbinom(n, 1, stats::plogis(drop(stz$Z %*% fit$gamma))) else rep(0L, n)
  mu <- if (stz$mixture) fit$mu1 else 0; sg <- if (stz$mixture) fit$sigma1 else 1
  theta <- stats::rnorm(n, ifelse(cls == 1, mu, 0), ifelse(cls == 1, sg, 1))
  X <- matrix(NA_integer_, n, J)
  for (j in seq_len(J)) {
    K <- stz$K[j]
    P <- grm_prob(fit$items$a[j], fit$d[[j]], theta)
    u <- stats::runif(n)
    x <- rowSums(u > t(apply(P, 1L, cumsum))[, -K, drop = FALSE])
    if (stz$mixture && stz$neg[j]) {
      mis <- cls == 1 & stats::runif(n) < fit$pi_full[j]
      x[mis] <- K - 1L - x[mis]
    }
    X[, j] <- x
  }
  colnames(X) <- fit$item_names
  as.data.frame(X)
}

#' Person-level results
#'
#' @param fit A fitted \code{revmix} object.
#' @return A data frame with the posterior probability of belonging to the
#'   inconsistent class, the modal class, and the expected a posteriori (EAP)
#'   trait estimate under the mixture.
#' @export
rm_scores <- function(fit) {
  post <- if (is.null(fit$posterior)) rep(0, fit$n) else fit$posterior
  data.frame(p_inconsistent = post, class = as.integer(post > 0.5), eap = fit$eap)
}
