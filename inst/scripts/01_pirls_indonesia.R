# PIRLS 2011 Indonesia: "Students Confident in Reading" (ASBR08A-G).
# Requires the IEA file asgidnr3.sav (not redistributed). Run from the package root
# with the file in data-raw/:  Rscript inst/scripts/01_pirls_indonesia.R
suppressMessages({library(revmix); library(haven)})
g <- as.data.frame(zap_labels(read_sav(Sys.getenv("PIRLS_ASG", "data-raw/asgidnr3.sav"))))
key <- function(v, rev) { x <- as.numeric(v); x[!x %in% 1:4] <- NA; if (rev) x else 5 - x }
conf <- c(A = FALSE, B = FALSE, C = TRUE, D = FALSE, E = TRUE, F = FALSE, G = TRUE)
for (k in names(conf)) g[[paste0("SC", k)]] <- key(g[[paste0("ASBR08", k)]], conf[[k]])
items <- paste0("SC", names(conf)); rev <- paste0("SC", names(conf)[conf])
pv <- sprintf("ASRREA0%d", 1:5)
cat("N =", nrow(g), "; respondents with no item response:", sum(rowSums(!is.na(g[items])) == 0), "\n")
cat("Weighted mean reading score:", round(weighted.mean(rowMeans(g[pv]), g$TOTWGT), 1), "\n")
r <- cor(g[items], use = "pairwise")
cat("Mean r positive-positive:", round(mean(r[!conf, !conf][upper.tri(r[!conf, !conf])]), 3),
    " negative-negative:", round(mean(r[conf, conf][upper.tri(r[conf, conf])]), 3),
    " positive-negative:", round(mean(r[!conf, conf]), 3), "\n")
set.seed(2011)
f0 <- rm_fit(g, items, rev, pi = "none", weights = "TOTWGT", se = FALSE)
pa <- rm_pv(g, items, rev, pv = pv, pi = "all", theta_class = "none", weights = "TOTWGT",
            cluster = "IDSCHOOL", nstart = 4)
pf <- rm_pv(g, items, rev, pv = pv, pi = "free", theta_class = "none", weights = "TOTWGT",
            cluster = "IDSCHOOL", nstart = 4)
fm <- suppressWarnings(rm_fit(transform(g, .pv = (ASRREA01 - 500) / 100), items, rev, pi = "free",
                              covariates = ~ .pv, theta_class = "mean", weights = "TOTWGT",
                              se = FALSE, nstart = 6))
print(pa); print(pf)
cmp <- rm_compare(f0, pa$fits[[1]], pf$fits[[1]])
print(cmp)
cat("Free class mean (PV1):", round(fm$loglik, 1), "mu1 =", round(fm$mu1, 2), "\n")
saveRDS(list(f0 = f0, pa = pa, pf = pf, fm = fm, cmp = cmp, r = r), "pirls_results.rds")
