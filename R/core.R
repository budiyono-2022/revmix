## Internal building blocks: GRM probabilities, parameter packing, likelihood.

# Category probabilities of the graded response model on a quadrature grid.
# a: slope; d: K-1 decreasing intercepts; theta: nodes. Returns Q x K matrix.
grm_prob <- function(a, d, theta) {
  Q <- length(theta)
  K <- length(d) + 1L
  z <- a * theta + matrix(d, Q, K - 1L, byrow = TRUE)
  Ps <- cbind(1, stats::plogis(z), 0)
  P <- Ps[, seq_len(K), drop = FALSE] - Ps[, -1L, drop = FALSE]
  pmax(P, 1e-300)
}

# Intercepts from the unconstrained parameterization (c1, log increments).
d_from_raw <- function(raw) {
  if (length(raw) == 1L) return(raw)
  raw[1L] - c(0, cumsum(exp(raw[-1L])))
}
raw_from_d <- function(d) {
  if (length(d) == 1L) return(d)
  inc <- pmax(-diff(d), 1e-4)
  c(d[1L], log(inc))
}

# Expected-count M-step for one GRM item. r: Q x K matrix of expected counts.
mstep_item <- function(r, theta, a, d) {
  K <- ncol(r)
  par0 <- c(a, raw_from_d(d))
  nll <- function(p) {
    P <- grm_prob(p[1L], d_from_raw(p[-1L]), theta)
    -sum(r * log(P))
  }
  gr <- function(p) {
    a <- p[1L]; raw <- p[-1L]; d <- d_from_raw(raw)
    Q <- length(theta)
    z <- a * theta + matrix(d, Q, K - 1L, byrow = TRUE)
    Pstar <- stats::plogis(z)
    W <- Pstar * (1 - Pstar)
    P <- grm_prob(a, d, theta)
    R <- r / P
    Wfull <- cbind(0, W, 0)
    ga <- sum(R * theta * (Wfull[, seq_len(K), drop = FALSE] - Wfull[, -1L, drop = FALSE]))
    gd <- colSums(W * (R[, -1L, drop = FALSE] - R[, -K, drop = FALSE]))
    J <- matrix(0, K - 1L, K - 1L)
    J[, 1L] <- 1
    if (K > 2L) for (m in 2:(K - 1L)) J[m:(K - 1L), m] <- -exp(raw[m])
    -c(ga, drop(gd %*% J))
  }
  fit <- stats::nlminb(par0, nll, gr, control = list(iter.max = 50, eval.max = 80))
  list(a = fit$par[1L], d = d_from_raw(fit$par[-1L]))
}

# Rescale weights to sum to the number of respondents.
norm_weights <- function(w, n) {
  if (is.null(w)) return(rep(1, n))
  if (length(w) != n || any(!is.finite(w)) || any(w < 0))
    stop("`weights` must be finite, non-negative and of length nrow(data).", call. = FALSE)
  w * n / sum(w)
}

# log-sum-exp by row
rlse <- function(M) {
  m <- apply(M, 1L, max)
  m + log(rowSums(exp(M - m)))
}

# Log-likelihood contributions of one item on the grid: N x Q matrix.
item_ll <- function(P, x) {
  out <- matrix(0, length(x), nrow(P))
  ok <- !is.na(x)
  out[ok, ] <- t(log(P))[x[ok] + 1L, , drop = FALSE]
  out
}

# Class-1 contribution for a negatively worded item: log((1-pi) P(x) + pi P(K-1-x)).
item_ll_mix <- function(P, x, pi) {
  K <- ncol(P)
  out <- matrix(0, length(x), nrow(P))
  rho <- matrix(0, length(x), nrow(P))
  ok <- !is.na(x)
  Pn <- t(P)[x[ok] + 1L, , drop = FALSE]
  Pr <- t(P)[K - x[ok], , drop = FALSE]
  mix <- (1 - pi) * Pn + pi * Pr
  out[ok, ] <- log(mix)
  rho[ok, ] <- pi * Pr / mix
  list(ll = out, rho = rho)
}

# Full marginal log-likelihood per person (used for scores and final values).
person_loglik <- function(st) {
  X <- st$X; J <- ncol(X); theta <- st$theta
  LL0 <- matrix(0, nrow(X), length(theta)); LL1 <- LL0
  for (j in seq_len(J)) {
    P <- grm_prob(st$a[j], st$d[[j]], theta)
    l <- item_ll(P, X[, j])
    LL0 <- LL0 + l
    if (st$neg[j] && st$mixture) {
      LL1 <- LL1 + item_ll_mix(P, X[, j], st$pi[j])$ll
    } else LL1 <- LL1 + l
  }
  lg0 <- log(st$g0)
  if (!st$mixture) return(rlse(sweep(LL0, 2L, lg0, "+")))
  eta <- stats::plogis(drop(st$Z %*% st$gamma))
  lg1 <- log(st$g1)
  A <- log1p(-eta) + rlse(sweep(LL0, 2L, lg0, "+"))
  B <- log(eta) + rlse(sweep(LL1, 2L, lg1, "+"))
  m <- pmax(A, B)
  m + log(exp(A - m) + exp(B - m))
}

class_prior <- function(theta, mu, sigma) {
  g <- stats::dnorm(theta, mu, sigma)
  g / sum(g)
}

# Column sums of M by response category: returns Q x K matrix.
cat_sums <- function(M, x, K) {
  rs <- rowsum(M, x, reorder = TRUE)
  out <- matrix(0, K, ncol(M))
  out[as.integer(rownames(rs)) + 1L, ] <- rs
  t(out)
}
