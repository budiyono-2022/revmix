# revmix

Mixture graded response models for misreading of reverse-worded items.

Rating scales that mix positively and negatively worded items often contain a
group of respondents who answer some reverse-worded items as if they were
positively worded, for example because they misread a negation or lack the
reading proficiency to process it. `revmix` estimates a two-class mixture
graded response model (GRM) in which:

- the **consistent class** answers every item according to a GRM;
- in the **inconsistent class**, each reverse-worded item *j* is misread with
  probability π_j, in which case the keyed response follows the same GRM with
  its categories reversed;
- π_j can be fixed to 1 (`pi = "all"`, the ordinal analogue of the constrained
  factor mixture model of Steinmann, Strietholt, & Braeken, 2022), estimated
  separately for each item (`pi = "free"`), or explained by item features
  such as wording type or content facet (`pi = "features"`);
- class membership is regressed on person covariates (for example, reading
  proficiency) in the same estimation step, with support for plausible values
  (Rubin's rules), sampling weights, and cluster-robust standard errors.

Estimation is by marginal maximum likelihood with an EM algorithm.

## Installation

```r
# install.packages("remotes")
remotes::install_github("budiyono-2022/revmix")
```

## Example

```r
library(revmix)
d <- moderation                      # bundled data, N = 4,197, 36 items
items <- grep("^[A-D][0-9]+$", names(d), value = TRUE)
rev <- c("A2","A3","A6","A7","A10","B12","B14","B16","B19",
         "C21","C22","C25","C26","C29","C30","D32","D34","D36")
d$female <- as.numeric(d$gender == "female")

f0 <- rm_fit(d, items, rev, pi = "none", se = FALSE)          # single-class GRM
fa <- rm_fit(d, items, rev, pi = "all",  covariates = ~ female)
ff <- rm_fit(d, items, rev, pi = "free", covariates = ~ female)
rm_compare(f0, fa, ff)
ff$pi                                 # item-specific misreading probabilities
head(rm_scores(ff))                   # posterior class probability and EAP
```

Main functions:

| Function | Purpose |
|---|---|
| `rm_fit()` | Fit the single-class GRM or a misreading mixture |
| `rm_pv()` | Fit across plausible values and pool with Rubin's rules |
| `rm_compare()` | Information criteria and likelihood-ratio tests |
| `rm_boot_lrt()` | Parametric bootstrap LRT for boundary restrictions |
| `rm_scores()` | Posterior class probabilities, classification, EAP scores |
| `rm_simulate()` | Generate data from the model |

## Reproducing the article

Scripts in `inst/scripts/` reproduce the empirical examples (PIRLS 2011
Indonesia, which requires the IEA data file, and the bundled religious
moderation data), the cross-checks against `mirt`, and the simulation study.

## Reference

Steinmann, I., Strietholt, R., & Braeken, J. (2022). A constrained factor
mixture analysis model for consistent and inconsistent respondents to
mixed-worded scales. *Psychological Methods, 27*(4), 667–702.
https://doi.org/10.1037/met0000392
