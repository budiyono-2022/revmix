# Cross-checks with mirt (Chalmers, 2012):
#  (a) the single-class GRM of revmix reproduces the mirt log-likelihood;
#  (b) the misreading mixture is compared with the usual competing account,
#      a bifactor GRM with an orthogonal wording factor on the reverse-worded items.
# Covariate-free models are compared so that both likelihoods refer to p(y).
# Run from the package root:  Rscript inst/scripts/04_method_factor.R
suppressMessages({library(revmix); library(mirt); library(haven)})

compare_accounts <- function(d, items, rev, weights = NULL, label) {
  keep <- rowSums(!is.na(d[items])) > 0
  d <- d[keep, ]
  w <- if (is.null(weights)) NULL else d[[weights]] / mean(d[[weights]])
  y <- as.data.frame(lapply(d[items], function(x) x - min(x, na.rm = TRUE)))
  m_grm <- mirt(y, 1, itemtype = "graded", survey.weights = w, verbose = FALSE,
                technical = list(NCYCLES = 3000))
  spec <- ifelse(items %in% rev, 1, NA)
  m_bif <- bfactor(y, spec, itemtype = "graded", survey.weights = w, verbose = FALSE,
                   quadpts = 41, technical = list(NCYCLES = 3000))
  r0 <- rm_fit(d, items, rev, pi = "none", weights = weights, se = FALSE)
  rf <- rm_fit(d, items, rev, pi = "free", weights = weights, se = FALSE, nstart = 4)
  ra <- rm_fit(d, items, rev, pi = "all", weights = weights, se = FALSE, nstart = 4)
  n <- sum(keep)
  bic <- function(ll, k) -2 * ll + log(n) * k
  out <- data.frame(
    data = label,
    model = c("GRM (mirt)", "GRM (revmix)", "Bifactor GRM, wording factor (mirt)",
              "Mixture, all misread (revmix)", "Mixture, item-specific (revmix)"),
    logLik = c(extract.mirt(m_grm, "logLik"), r0$loglik, extract.mirt(m_bif, "logLik"),
               ra$loglik, rf$loglik),
    npar = c(extract.mirt(m_grm, "nest"), r0$npar, extract.mirt(m_bif, "nest"), ra$npar, rf$npar))
  out$BIC <- bic(out$logLik, out$npar)
  out
}

## PIRLS 2011 Indonesia (requires data-raw/asgidnr3.sav)
g <- as.data.frame(zap_labels(read_sav(Sys.getenv("PIRLS_ASG", "data-raw/asgidnr3.sav"))))
key <- function(v, rev) { x <- as.numeric(v); x[!x %in% 1:4] <- NA; if (rev) x else 5 - x }
conf <- c(A = FALSE, B = FALSE, C = TRUE, D = FALSE, E = TRUE, F = FALSE, G = TRUE)
for (k in names(conf)) g[[paste0("SC", k)]] <- key(g[[paste0("ASBR08", k)]], conf[[k]])
items_p <- paste0("SC", names(conf)); rev_p <- paste0("SC", names(conf)[conf])
set.seed(2011)
tp <- compare_accounts(g, items_p, rev_p, weights = "TOTWGT", label = "PIRLS")

## Religious moderation data
d <- moderation
items_m <- grep("^[A-D][0-9]+$", names(d), value = TRUE)
rev_m <- c("A2", "A3", "A6", "A7", "A10", "B12", "B14", "B16", "B19", "C21", "C22", "C25",
           "C26", "C29", "C30", "D32", "D34", "D36")
set.seed(2026)
tm <- compare_accounts(d, items_m, rev_m, label = "Moderation")

res <- rbind(tp, tm)
print(res, digits = 8)
saveRDS(res, "method_factor_results.rds")
