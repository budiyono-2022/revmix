#' Simulate keyed responses from the misreading mixture GRM
#'
#' Generates rating-scale responses (already keyed so that higher means more
#' of the trait) from the model fitted by \code{\link{rm_fit}}. Reverse-worded
#' items are the last \code{n_reversed} items.
#'
#' @param n Number of respondents.
#' @param n_items,n_reversed Number of items and of reverse-worded items.
#' @param n_cat Number of response categories.
#' @param a Slopes (scalar or vector); default draws from U(1, 2).
#' @param prevalence Marginal prevalence of the inconsistent class when
#'   \code{gamma} is \code{NULL}.
#' @param gamma Optional coefficients (intercept, slope) of a logistic model of
#'   class membership on a standard normal covariate \code{z}, which is
#'   returned in the data.
#' @param pi Misreading probabilities of the reverse-worded items (recycled).
#' @param mu1 Trait mean in the inconsistent class.
#' @return A list with \code{data}, \code{items}, \code{reversed}, the true
#'   \code{class}, \code{theta}, and generating parameters.
#' @examples
#' sim <- rm_simulate(n = 200)
#' head(sim$data)
#' @export
rm_simulate <- function(n = 1000, n_items = 12, n_reversed = 6, n_cat = 4, a = NULL,
                        prevalence = 0.15, gamma = NULL, pi = 0.8, mu1 = 0) {
  if (n_reversed < 1 || n_reversed >= n_items) stop("Need 1 <= n_reversed < n_items.", call. = FALSE)
  if (is.null(a)) a <- stats::runif(n_items, 1, 2)
  a <- rep_len(a, n_items)
  d <- lapply(seq_len(n_items), function(j) sort(stats::qnorm(seq_len(n_cat - 1L) / n_cat) * 1.6 +
                                                    stats::runif(1, -0.6, 0.6), decreasing = TRUE))
  z <- stats::rnorm(n)
  eta <- if (is.null(gamma)) rep(prevalence, n) else stats::plogis(gamma[1] + gamma[2] * z)
  cls <- stats::rbinom(n, 1, eta)
  theta <- stats::rnorm(n, ifelse(cls == 1, mu1, 0), 1)
  neg <- c(rep(FALSE, n_items - n_reversed), rep(TRUE, n_reversed))
  pi_full <- numeric(n_items); pi_full[neg] <- rep_len(pi, n_reversed)
  X <- matrix(NA_integer_, n, n_items)
  for (j in seq_len(n_items)) {
    P <- grm_prob(a[j], d[[j]], theta)
    u <- stats::runif(n)
    x <- rowSums(u > t(apply(P, 1L, cumsum))[, -n_cat, drop = FALSE])
    if (neg[j]) {
      mis <- cls == 1 & stats::runif(n) < pi_full[j]
      x[mis] <- n_cat - 1L - x[mis]
    }
    X[, j] <- x
  }
  items <- sprintf("i%02d", seq_len(n_items))
  colnames(X) <- items
  list(data = data.frame(X, z = z), items = items, reversed = items[neg], class = cls,
       theta = theta, a = a, d = d, pi = pi_full[neg], gamma = gamma,
       prevalence = mean(eta), mu1 = mu1)
}
