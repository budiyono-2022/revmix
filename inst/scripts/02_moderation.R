# Bundled religious moderation data (Pratama et al., 2026).
# Run from the package root:  Rscript inst/scripts/02_moderation.R
suppressMessages(library(revmix))
d <- moderation
items <- grep("^[A-D][0-9]+$", names(d), value = TRUE)
rev <- c("A2", "A3", "A6", "A7", "A10", "B12", "B14", "B16", "B19", "C21", "C22", "C25",
         "C26", "C29", "C30", "D32", "D34", "D36")
d$female <- as.numeric(d$gender == "female")
d$islamic <- as.numeric(d$school_type == "islamic")
d$vocational <- as.numeric(d$school_type == "vocational")
d$muslim <- as.numeric(d$religion == "Islam")
cov <- ~ female + islamic + vocational + muslim
set.seed(2026)
m0 <- rm_fit(d, items, rev, pi = "none", se = FALSE)
ma <- rm_fit(d, items, rev, pi = "all", covariates = cov, nstart = 4)
mf <- rm_fit(d, items, rev, pi = "free", covariates = cov, nstart = 4)
feat <- data.frame(item = rev, dimension = factor(substr(rev, 1, 1)))
mx <- rm_fit(d, items, rev, pi = "features", item_features = feat, features = ~ dimension,
             covariates = cov, nstart = 1, start = mf)
print(ma); print(mf); print(mx)
print(rm_compare(m0, ma, mf)); print(rm_compare(mx, mf))
sc <- rm_scores(mf)
cat("Respondents classified inconsistent:", round(mean(sc$class), 3), "\n")
cat("Correlation of EAP (mixture) with EAP (single class):", round(cor(sc$eap, m0$eap), 3), "\n")
saveRDS(list(m0 = m0, ma = ma, mf = mf, mx = mx), "moderation_results.rds")
