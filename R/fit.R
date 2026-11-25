#' Fit a mixture graded response model for misreading of reverse-worded items
#'
#' Estimates a two-class mixture graded response model (GRM) for a
#' unidimensional rating scale that mixes positively and negatively (reverse)
#' worded items. Responses must already be keyed so that higher categories
#' indicate more of the trait for every item (that is, reverse-worded items have
#' been reverse-scored). In the consistent class every item follows the GRM. In
#' the inconsistent class each reverse-worded item \eqn{j} is misread with
#' probability \eqn{\pi_j}: the respondent then answers as if the item were
#' positively worded, so the keyed response follows the GRM with its categories
#' reversed. Item parameters are shared by the two classes.
#'
#' Three specifications of the misreading probabilities are available.
#' \code{pi = "all"} fixes \eqn{\pi_j = 1}, which is the ordinal analogue of the
#' constrained factor mixture model of Steinmann, Strietholt, and Braeken
#' (2022) in which inconsistent respondents misread every reverse-worded item.
#' \code{pi = "free"} estimates a separate \eqn{\pi_j} for each reverse-worded
#' item. \code{pi = "features"} models \eqn{\mathrm{logit}(\pi_j) = w_j'\delta}
#' from item characteristics supplied in \code{item_features}. The probability
#' of belonging to the inconsistent class follows a logistic regression on the
#' person covariates in \code{covariates} (for example, reading proficiency),
#' estimated in the same step as the measurement model. The trait is standard
#' normal in the consistent class; its mean (and optionally variance) is free in
#' the inconsistent class.
#'
#' Estimation is by marginal maximum likelihood with an EM algorithm on a fixed
#' rectangular quadrature grid. Standard errors are based on the cross-product
#' of per-respondent scores obtained by central differences; with
#' \code{cluster}, cluster-robust standard errors use the score cross-product as
#' the bread of a sandwich estimator. Sampling weights enter as a
#' pseudo-likelihood and are rescaled to sum to the sample size.
#'
#' @param data A data frame.
#' @param items Character vector of item names (keyed, integer categories).
#' @param reversed Character vector of reverse-worded items (a subset of
#'   \code{items}).
#' @param pi One of \code{"free"}, \code{"all"}, \code{"features"}, or
#'   \code{"none"} (a single-class GRM).
#' @param covariates One-sided formula for class membership, e.g.
#'   \code{~ reading}; variables are taken from \code{data}. \code{NULL} gives
#'   an intercept-only model.
#' @param item_features A data frame with one row per reverse-worded item (in
#'   the order of \code{reversed}) used when \code{pi = "features"}, together
#'   with \code{features}, a one-sided formula.
#' @param features One-sided formula for the misreading model.
#' @param theta_class Distribution of the trait in the inconsistent class:
#'   \code{"mean"} (free mean, default), \code{"meanvar"}, or \code{"none"}.
#' @param weights Optional name of a weight variable in \code{data}.
#' @param cluster Optional name of a cluster variable in \code{data} (for
#'   cluster-robust standard errors).
#' @param se Logical; compute standard errors?
#' @param nquad Number of quadrature nodes on \eqn{[-6, 6]}.
#' @param maxit,tol Maximum EM iterations and convergence tolerance on the
#'   change in log-likelihood.
#' @param nstart Number of random starts for the mixture parameters; the
#'   solution with the highest log-likelihood is kept.
#' @param start Optional fitted \code{revmix} object supplying starting values.
#' @param verbose Print progress?
#' @return An object of class \code{revmix}.
#' @references Steinmann, I., Strietholt, R., & Braeken, J. (2022). A
#'   constrained factor mixture analysis model for consistent and inconsistent
#'   respondents to mixed-worded scales. \emph{Psychological Methods, 27}(4),
#'   667--702. \doi{10.1037/met0000392}
#' @examples
#' set.seed(1)
#' sim <- rm_simulate(n = 500, n_items = 6, n_reversed = 3, prevalence = 0.25,
#'                    pi = c(0.9, 0.6, 0.3))
#' fit <- rm_fit(sim$data, sim$items, sim$reversed, se = FALSE, nquad = 21,
#'               nstart = 1, maxit = 200)
#' fit
#' @export
rm_fit <- function(data, items, reversed, pi = c("free", "all", "features", "none"),
                   covariates = NULL, item_features = NULL, features = NULL,
                   theta_class = c("mean", "meanvar", "none"), weights = NULL,
                   cluster = NULL, se = TRUE, nquad = 61, maxit = 2000, tol = 1e-6,
                   nstart = 3, start = NULL, verbose = FALSE) {
  pi <- match.arg(pi)
  theta_class <- match.arg(theta_class)
  cl <- match.call()
  if (!all(items %in% names(data))) stop("Some `items` are not in `data`.", call. = FALSE)
  if (!length(reversed) || !all(reversed %in% items))
    stop("`reversed` must be a non-empty subset of `items`.", call. = FALSE)
  X <- as.matrix(data[, items, drop = FALSE])
  storage.mode(X) <- "double"
  if (any(X != round(X), na.rm = TRUE)) stop("Item responses must be integers.", call. = FALSE)
  mins <- apply(X, 2L, min, na.rm = TRUE)
  X <- sweep(X, 2L, mins, "-")
  K <- apply(X, 2L, max, na.rm = TRUE) + 1L
  if (any(K < 2L)) stop("Every item needs at least two observed categories.", call. = FALSE)
  keep <- rowSums(!is.na(X)) > 0
  n <- nrow(X)
  neg <- items %in% reversed
  w <- norm_weights(if (is.null(weights)) NULL else data[[weights]], n)
  mixture <- pi != "none"
  Z <- if (is.null(covariates)) matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)")) else
    stats::model.matrix(covariates, stats::model.frame(covariates, data, na.action = stats::na.pass))
  if (anyNA(Z)) stop("Missing values in `covariates` are not allowed.", call. = FALSE)
  Wf <- NULL
  if (pi == "features") {
    if (is.null(item_features) || is.null(features))
      stop("`pi = \"features\"` needs `item_features` and `features`.", call. = FALSE)
    if (nrow(item_features) != length(reversed))
      stop("`item_features` needs one row per reversed item.", call. = FALSE)
    Wf <- stats::model.matrix(features, item_features)
  }
  X <- X[keep, , drop = FALSE]; Z <- Z[keep, , drop = FALSE]; w <- w[keep]
  w <- w * length(w) / sum(w)
  clus <- if (is.null(cluster)) NULL else data[[cluster]][keep]
  theta <- seq(-6, 6, length.out = nquad)

  st <- list(X = X, neg = neg, theta = theta, g0 = class_prior(theta, 0, 1),
             g1 = class_prior(theta, 0, 1), mu1 = 0, sigma1 = 1, Z = Z, w = w,
             mixture = FALSE, pi = ifelse(neg, 0, 0), gamma = c(-20, rep(0, ncol(Z) - 1L)),
             delta = NULL, Wf = Wf, K = K, pi_type = pi, theta_class = theta_class)
  # initial item values
  st$a <- rep(1, length(items))
  st$d <- lapply(K, function(k) stats::qlogis(seq(0.9, 0.1, length.out = k - 1L)))
  if (!is.null(start) && inherits(start, "revmix")) {
    st$a <- start$items$a; st$d <- start$d
  } else {
    st <- em_run(st, maxit = 300, tol = tol, verbose = FALSE)
  }
  if (mixture) {
    st$mixture <- TRUE
    base <- st
    best <- NULL
    starts <- max(1L, nstart)
    for (s in seq_len(starts)) {
      st <- base
      if (!is.null(start) && inherits(start, "revmix") && s == 1L && start$pi_type != "none") {
        st$gamma <- start$gamma; st$mu1 <- start$mu1; st$sigma1 <- start$sigma1
        st$pi <- start$pi_full
        st$delta <- start$delta
      } else {
        p0 <- if (s == 1L) 0.10 else stats::runif(1, 0.05, 0.30)
        st$gamma <- c(stats::qlogis(p0), rep(0, ncol(Z) - 1L))
        st$pi <- ifelse(neg, if (pi == "all") 1 else if (s == 1L) 0.8 else stats::runif(sum(neg), 0.4, 0.95), 0)
        if (pi == "features") st$delta <- c(stats::qlogis(0.8), rep(0, ncol(Wf) - 1L))
      }
      st$g1 <- class_prior(theta, st$mu1, st$sigma1)
      st <- em_run(st, maxit = maxit, tol = tol, verbose = verbose)
      if (is.null(best) || st$loglik > best$loglik) best <- st
    }
    st <- best
  } else {
    st <- em_run(st, maxit = maxit, tol = tol, verbose = verbose)
  }
  out <- finish_fit(st, items, reversed, n_dropped = sum(!keep), clus = clus, se = se, call = cl)
  if (mixture && theta_class != "none" && abs(out$mu1) > 2)
    warning("The trait mean of the inconsistent class is more than 2 SD from the consistent class (",
            round(out$mu1, 2), "); the mixture may be absorbing multidimensionality. ",
            "Consider theta_class = \"none\" and compare solutions.", call. = FALSE)
  out$mins <- mins
  out$item_features <- item_features
  out$features <- features
  out
}

# EM iterations ------------------------------------------------------------
em_run <- function(st, maxit, tol, verbose) {
  X <- st$X; theta <- st$theta; J <- ncol(X); n <- nrow(X); w <- st$w
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    # E-step
    P <- lapply(seq_len(J), function(j) grm_prob(st$a[j], st$d[[j]], theta))
    LL0 <- matrix(0, n, length(theta)); LL1 <- LL0
    rho <- vector("list", J)
    for (j in seq_len(J)) {
      l <- item_ll(P[[j]], X[, j])
      LL0 <- LL0 + l
      if (st$mixture && st$neg[j]) {
        m <- item_ll_mix(P[[j]], X[, j], st$pi[j])
        LL1 <- LL1 + m$ll; rho[[j]] <- m$rho
      } else LL1 <- LL1 + l
    }
    if (st$mixture) {
      eta <- stats::plogis(drop(st$Z %*% st$gamma))
      A <- sweep(LL0, 2L, log(st$g0), "+") + log1p(-eta)
      B <- sweep(LL1, 2L, log(st$g1), "+") + log(eta)
      m <- pmax(apply(A, 1L, max), apply(B, 1L, max))
      eA <- exp(A - m); eB <- exp(B - m)
      L <- rowSums(eA) + rowSums(eB)
      post0 <- eA / L; post1 <- eB / L
      li <- m + log(L)
    } else {
      A <- sweep(LL0, 2L, log(st$g0), "+")
      m <- apply(A, 1L, max); eA <- exp(A - m); L <- rowSums(eA)
      post0 <- eA / L; post1 <- 0 * post0
      li <- m + log(L)
    }
    ll <- sum(w * li)
    if (verbose && it %% 25 == 0) message(sprintf("iter %d  loglik %.4f", it, ll))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
    # M-step: items
    wp0 <- post0 * w; wp1 <- post1 * w
    for (j in seq_len(J)) {
      x <- X[, j]; ok <- !is.na(x); Kj <- st$K[j]
      if (st$mixture && st$neg[j]) {
        Mn <- wp0[ok, , drop = FALSE] + wp1[ok, , drop = FALSE] * (1 - rho[[j]][ok, , drop = FALSE])
        Mr <- wp1[ok, , drop = FALSE] * rho[[j]][ok, , drop = FALSE]
        r <- cat_sums(Mn, x[ok], Kj) + cat_sums(Mr, Kj - 1L - x[ok], Kj)
      } else {
        Mn <- wp0[ok, , drop = FALSE] + wp1[ok, , drop = FALSE]
        r <- cat_sums(Mn, x[ok], Kj)
      }
      upd <- mstep_item(r, theta, st$a[j], st$d[[j]])
      st$a[j] <- upd$a; st$d[[j]] <- upd$d
    }
    if (st$mixture) {
      p1 <- pmin(pmax(rowSums(post1), 0), 1)
      # class membership
      g <- suppressWarnings(stats::glm.fit(st$Z, p1, weights = w, start = st$gamma,
                                           family = stats::quasibinomial()))
      st$gamma <- g$coefficients
      # misreading probabilities
      if (st$pi_type %in% c("free", "features")) {
        idx <- which(st$neg)
        Aj <- Bj <- numeric(length(idx))
        for (k in seq_along(idx)) {
          j <- idx[k]; ok <- !is.na(X[, j])
          Aj[k] <- sum(wp1[ok, , drop = FALSE] * rho[[j]][ok, , drop = FALSE])
          Bj[k] <- sum(wp1[ok, , drop = FALSE])
        }
        if (st$pi_type == "free") {
          st$pi[idx] <- pmin(pmax(Aj / Bj, 1e-6), 1 - 1e-6)
        } else {
          gf <- suppressWarnings(stats::glm.fit(st$Wf, pmin(pmax(Aj / Bj, 0), 1), weights = Bj, start = st$delta,
                                                family = stats::quasibinomial()))
          st$delta <- gf$coefficients
          st$pi[idx] <- pmin(pmax(stats::plogis(drop(st$Wf %*% st$delta)), 1e-6), 1 - 1e-6)
        }
      }
      # trait distribution in the inconsistent class
      if (st$theta_class != "none") {
        s1 <- sum(wp1)
        st$mu1 <- sum(wp1 %*% theta) / s1
        if (st$theta_class == "meanvar")
          st$sigma1 <- sqrt(max(sum(wp1 %*% (theta - st$mu1)^2) / s1, 0.05))
        st$g1 <- class_prior(theta, st$mu1, st$sigma1)
      }
    }
  }
  st$loglik <- ll
  st$iterations <- it
  st$converged <- it < maxit
  st$post0 <- post0; st$post1 <- post1
  st$li <- li
  st
}
