## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(fig.align = "center")


## ---- message = FALSE---------------------------------------------------------
library(assertthat)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(viridis)
library(nloptr)


## -----------------------------------------------------------------------------
simulate_confusion_matrix <- function(
  n, p, reference_sensitivity, reference_specificity, test_sensitivity,
  test_specificity) {
  # assert that arguments are valid
  assert_that(is.count(n), is.number(p), p >= 0, p <= 1,
    is.number(reference_sensitivity), is.number(reference_specificity),
    reference_sensitivity >= 0, reference_sensitivity <= 1,
    reference_specificity >= 0, reference_specificity <= 1,
    is.number(test_sensitivity), is.number(test_specificity),
    test_sensitivity >= 0, test_specificity <= 1,
    test_sensitivity >= 0, test_specificity <= 1)
  # rename variables to follow Staquet (eqns 1--4)
  x <- n * p
  y <- n * (1 - p)
  sr <- reference_sensitivity
  sn <- test_sensitivity
  spr <- reference_specificity
  spn <- test_specificity
  ma <- (x * sr * sn) + (y * (1 - spn) * (1 - spr))
  mb <- (x * (1 - sn) * sr) + (y * spn * (1 - spr))
  mc <- (x * sn * (1 - sr)) + (y * spr * (1 - spn))
  md <- (x * (1 - sr) * (1 - sn)) + (y * spr * spn)
  # create confusion matrix following Staquet
  cm <- matrix(c(ma, mb, mc, md), ncol = 2,
         dimnames = list(
           c("test_presence", "test_absence"),
           c("reference_presence", "reference_absence")))
  # round values in matrix
  cm[] <- round(cm)
  # check for valid values
  assert_that(all(c(cm) >= 0), sum(cm) >= (0.9 * n), sum(cm) <= (1.1 * n))
  # return matrix
  cm
}


## -----------------------------------------------------------------------------
# simulate confusion matrix where species has a 50% prevalence, and
# both reference data and model have high sensitivity and specificity values
simulate_confusion_matrix(
  n = 1000, p = 0.5,
  reference_sensitivity = 0.9, reference_specificity = 0.9,
  test_sensitivity = 0.89, test_specificity = 0.89)

# simulate confusion matrix where species has a 10% prevalence, and
# both reference data and model have high sensitivity and specificity values
simulate_confusion_matrix(
  n = 1000, p = 0.1,
  reference_sensitivity = 0.9, reference_specificity = 0.9,
  test_sensitivity = 0.89,  test_specificity = 0.89)

# simulate confusion matrix where species has a 10% prevalence,
# the reference data has high sensitivity and specificity, and the model
# has low high sensitivity and specificity
simulate_confusion_matrix(
  n = 1000, p = 0.1,
  reference_sensitivity = 0.9, reference_specificity = 0.9,
  test_sensitivity = 0.6, test_specificity = 0.6)

# simulate confusion matrix where species has a 10% prevalence,
# and both the reference data and model have low sensitivity and specificity
simulate_confusion_matrix(
  n = 1000, p = 0.1,
  reference_sensitivity = 0.7, reference_specificity = 0.7,
  test_sensitivity = 0.6, test_specificity = 0.6)


## -----------------------------------------------------------------------------
staquet_performance <- function(
  x, reference_sensitivity, reference_specificity) {
  # assert that arguments are valid
  assert_that(is.matrix(x), ncol(x) == 2, nrow(x) == 2,
    is.number(reference_sensitivity), is.number(reference_specificity),
    reference_sensitivity >= 0, reference_sensitivity <= 1,
    reference_specificity >= 0, reference_specificity <= 1)
  # rename variables to follow Staquet (eqns 1--9)
  a <- x[1, 1]
  b <- x[2, 1]
  c <- x[1, 2]
  d <- x[2, 2]
  n <- a + b + c + d
  s1r <- reference_sensitivity
  s2r <- reference_specificity
  # estimate correct sensitivity and specificity
  cor_sens <- ceiling(((a + c) * s2r) - c) / ceiling((n * (s2r - 1)) + (a + b))
  cor_spec <- ceiling(((b + d) * s1r) - b) / ceiling((n * s1r) - (a + b))
  # return result
  c(sensitivity = cor_sens, specificity = cor_spec)
}

haverford_performance <- function(
  x, reference_sensitivity, reference_specificity) {
  # assert that arguments are valid
  assert_that(is.matrix(x), ncol(x) == 2, nrow(x) == 2,
    is.number(reference_sensitivity), is.number(reference_specificity),
    reference_sensitivity >= 0, reference_sensitivity <= 1,
    reference_specificity >= 0, reference_specificity <= 1)
  # calculate naive estimates for sensitivity and specificity
  # (assuming that reference data are perfect)
  naive_sens <- x[1, 1] / sum(x[, 1])
  naive_spec <- x[2, 2] / sum(x[, 2])
  # apply a correction to these naive estimates based on the
  # known sensitivity and specificity of the reference data
  cor_sens <- (naive_sens * reference_sensitivity) +
              ((1 - naive_spec) * (1 - reference_sensitivity))
  cor_spec <- (naive_spec * reference_specificity) +
              ((1 - naive_sens) * (1 - reference_specificity))
  # return result
  c(sensitivity = cor_sens, specificity = cor_spec)
}

maxlik_performance <- function(
  x, reference_sensitivity, reference_specificity) {
  # assert that arguments are valid
  assert_that(is.matrix(x), ncol(x) == 2, nrow(x) == 2,
    is.number(reference_sensitivity), is.number(reference_specificity),
    reference_sensitivity >= 0, reference_sensitivity <= 1,
    reference_specificity >= 0, reference_specificity <= 1)
  # redefine variables to follow source
  seR <- reference_sensitivity
  spR <- reference_specificity
  s00 <- x[2, 2]
  s10 <- x[2, 1]
  s01 <- x[1, 2]
  s11 <- x[1, 1]
  # define maximum likelihood function
  f <- function(x) {
    prec <- 1e+4 # precision for calculations
    seT <- x[1]
    spT <- x[2]
    Pi <- x[3]
    t1 <- (seR * seT * Pi) + ((1 - spR) * (1 - spT) * (1 - Pi))
    t2 <- (seR * (1 - seT) * Pi) + ((1 - spR) * spT * (1 - Pi))
    t3 <- ((1 - seR) * seT * Pi) + (spR * (1 - spT) * (1 - Pi))
    t4 <- ((1 - seR) * (1 - seT) * Pi) + (spR * spT * (1 - Pi))
    # calculate negative log likelihood
    # use Rmpfr package to account for really small numbers in calculations
    -as.numeric(sum(log(c(Rmpfr::mpfr(t1, prec) ^ Rmpfr::mpfr(s11, prec),
                          Rmpfr::mpfr(t2, prec) ^ Rmpfr::mpfr(s10, prec),
                          Rmpfr::mpfr(t3, prec) ^ Rmpfr::mpfr(s01, prec),
                          Rmpfr::mpfr(t4, prec) ^ Rmpfr::mpfr(s00, prec)))))
  }
  # estimate test sensitivity and specificity by minimizing negative loglik
  res <- nloptr::bobyqa(c(0.9, 0.9, 0.5), f,
                        lower = rep(1e-5, 3), upper = rep(1 - 1e-5, 3))
  # return result
  c(sensitivity = res$par[[1]], specificity = res$par[[2]])
}


## -----------------------------------------------------------------------------
# create confusion matrix
conf_matrix <- matrix(c(120, 180, 220, 480), ncol = 2, nrow = 2,
                      dimnames = list(
                        c("test_presence", "test_absence"),
                        c("reference_presence", "reference_absence")))
print(conf_matrix)

# define reference sensitivity and specificity
ref_sens <- 0.75
ref_spec <- 0.75

# estimate model sensitivity and specificity using Staquet method
staquet_performance(conf_matrix, ref_sens, ref_spec)

# estimate model sensitivity and specificity using Haverford method
haverford_performance(conf_matrix, ref_sens, ref_spec)

# estimate model sensitivity and specificity using maximum likelihood method
maxlik_performance(conf_matrix, ref_sens, ref_spec)


## -----------------------------------------------------------------------------
# set parameters for simulation
number_obs <- 1000
prevalence <- 0.05
ref_sens <- 0.77
ref_spec <- 0.999
test_sens <- seq(0, 1, length.out = 100)
test_spec <- 0.9

# simulate confusion matrices
sens_cm <- lapply(test_sens, simulate_confusion_matrix, n = number_obs,
  p = prevalence, reference_sensitivity = ref_sens,
  reference_specificity = ref_spec, test_specificity = test_spec)

# calculate corrected sensitivity and specificity using Staquet method
sens_staquet_results <- sapply(sens_cm, staquet_performance,
  reference_sensitivity = ref_sens, reference_specificity = ref_spec)

# calculate corrected sensitivity and specificity using Haverford method
sens_haverford_results <- sapply(sens_cm, haverford_performance,
  reference_sensitivity = ref_sens, reference_specificity = ref_spec)

# calculate corrected sensitivity and specificity using max likelihood method
sens_maxlik_results <- sapply(sens_cm, maxlik_performance,
  reference_sensitivity = ref_sens, reference_specificity = ref_spec)

# compile results
sens_results <-
  rbind(sens_staquet_results, sens_haverford_results, sens_maxlik_results)  %>%
  t() %>%
  as.data.frame() %>%
  setNames(c("staquet_sensitivity", "staquet_specificity",
             "haverford_sensitivity", "haverford_specificity",
             "maxlik_sensitivity", "maxlik_specificity")) %>%
  mutate(true_sensitivity = test_sens)

# plot estimates of sensitivity
p1 <-
  ggplot(sens_results,
         aes(x = true_sensitivity, y = staquet_sensitivity)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  ylab("Staquet method for estimating test sensitivity") +
  xlab("True test sensitivity") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
p2 <-
  ggplot(sens_results,
         aes(x = true_sensitivity, y = haverford_sensitivity)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  ylab("Haverford method for estimating test sensitivity") +
  xlab("True test sensitivity") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
p3 <-
  ggplot(sens_results,
         aes(x = true_sensitivity, y = maxlik_sensitivity)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  ylab("Maximum likelihood method for estimating test sensitivity") +
  xlab("True test sensitivity") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
grid.arrange(p1, p2, p3, nrow = 1)


## -----------------------------------------------------------------------------
# set parameters for simulation
number_obs <- 1000
prevalence <- 0.05
ref_sens <- 0.77
ref_spec <- 0.999
test_sens <- 0.75
test_spec <- seq(0, 1, length.out = 100)

# simulate confusion matrices
spec_cm <- lapply(test_spec, simulate_confusion_matrix, n = number_obs,
  p = prevalence, reference_sensitivity = ref_sens,
  reference_specificity = ref_spec, test_sensitivity = test_sens)

# calculate corrected sensitivity and specificity using Staquet method
spec_staquet_results <- sapply(spec_cm, staquet_performance,
  reference_sensitivity = ref_sens, reference_specificity = ref_spec)

# calculate corrected sensitivity and specificity using Haverford method
spec_haverford_results <- sapply(spec_cm, haverford_performance,
  reference_sensitivity = ref_sens, reference_specificity = ref_spec)

# calculate corrected sensitivity and specificity using max likelihood method
spec_maxlik_results <- sapply(spec_cm, maxlik_performance,
  reference_sensitivity = ref_sens, reference_specificity = ref_spec)

# compile results
spec_results <-
  rbind(spec_staquet_results, spec_haverford_results, spec_maxlik_results)  %>%
  t() %>%
  as.data.frame() %>%
  setNames(c("staquet_sensitivity", "staquet_specificity",
             "haverford_sensitivity", "haverford_specificity",
             "maxlik_sensitivity", "maxlik_specificity")) %>%
  mutate(true_specificity = test_spec)

# plot estimates of specificity
p1 <-
  ggplot(spec_results,
         aes(x = true_specificity, y = staquet_specificity)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  ylab("Staquet method estimate for estimating test specificity") +
  xlab("True test specificity") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
p2 <-
  ggplot(spec_results,
         aes(x = true_specificity, y = haverford_specificity)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  ylab("Haverford method for estimating test specificity") +
  xlab("True test specificity") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
p3 <-
  ggplot(spec_results,
         aes(x = true_specificity, y = maxlik_specificity)) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  ylab("Maximum likelihood method for estimating test specificity") +
  xlab("True test specificity") +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
grid.arrange(p1, p2, p3, nrow = 1)


## -----------------------------------------------------------------------------
# create confusion matrix
conf_matrix <- matrix(c(1005, 134, 147, 1668), ncol = 2, nrow = 2,
                      dimnames = list(
                        c("test_presence", "test_absence"),
                        c("reference_presence", "reference_absence")))
print(conf_matrix)

# define reference sensitivity and specificity
ref_sens <- 0.98
ref_spec <- 0.85

# estimate model sensitivity and specificity using Staquet method
staquet_performance(conf_matrix, ref_sens, ref_spec)

# estimate model sensitivity and specificity using Haverford method
haverford_performance(conf_matrix, ref_sens, ref_spec)

# estimate model sensitivity and specificity using maximum likelihood method
maxlik_performance(conf_matrix, ref_sens, ref_spec)


## ---- results = "hide"--------------------------------------------------------
# generate different combinations of parameters
model_params <- expand.grid(model_sens = seq(1e-5, 1 - 1e-5, length.out = 25),
                            model_spec = seq(1e-5, 1 - 1e-5, length.out = 25))

# calculate log-likelihood for each parameter combination
model_params$loglik <-
  plyr::laply(seq_len(nrow(model_params)), .progress = "text", function(i) {
  # set parameters
  seR <- ref_sens
  spR <- ref_spec
  seT <- model_params[[1]][[i]]
  spT <- model_params[[2]][[i]]
  s00 <- conf_matrix[2, 2]
  s10 <- conf_matrix[2, 1]
  s01 <- conf_matrix[1, 2]
  s11 <- conf_matrix[1, 1]
  # define maximum likelihood function
  f <- function(x) {
    prec <- 1e+4 # precision for calculations
    Pi <- x[1]
    t1 <- (seR * seT * Pi) + ((1 - spR) * (1 - spT) * (1 - Pi))
    t2 <- (seR * (1 - seT) * Pi) + ((1 - spR) * spT * (1 - Pi))
    t3 <- ((1 - seR) * seT * Pi) + (spR * (1 - spT) * (1 - Pi))
    t4 <- ((1 - seR) * (1 - seT) * Pi) + (spR * spT * (1 - Pi))
    # calculate negative log likelihood
    # use Rmpfr package to account for really small numbers in calculations
    -as.numeric(sum(log(c(Rmpfr::mpfr(t1, prec) ^ Rmpfr::mpfr(s11, prec),
                          Rmpfr::mpfr(t2, prec) ^ Rmpfr::mpfr(s10, prec),
                          Rmpfr::mpfr(t3, prec) ^ Rmpfr::mpfr(s01, prec),
                          Rmpfr::mpfr(t4, prec) ^ Rmpfr::mpfr(s00, prec)))))
  }
  # find loglik for best supported prevalence given
  # test sensitivity
  res <- nloptr::bobyqa(c(0.5), f, lower = c(1e-5), upper = c(1 - 1e-5))
  # return result
  -1 * res$value
})


## -----------------------------------------------------------------------------
# create plot
ggplot(data = model_params,
       aes(x = model_sens, y = model_spec, fill = loglik)) +
geom_tile() +
xlab("Estimate of model sensitivity") +
xlab("Estimate of model specificity") +
scale_fill_viridis("Log Likelihood") +
theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
sessionInfo()

