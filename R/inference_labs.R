# Reusable calculations for the teaching laboratories in Chapters 7--9.
# The functions use only base R so the book can render in a clean session.

binomial_lr_deviance <- function(p, y, n) {
  phat <- y / n
  success_term <- if (y == 0) 0 else y * log(phat / p)
  failure_term <- if (y == n) 0 else (n - y) * log((1 - phat) / (1 - p))
  2 * (success_term + failure_term)
}

binomial_intervals <- function(y, n, conf.level = 0.95) {
  stopifnot(length(y) == 1L, length(n) == 1L, n > 0, y >= 0, y <= n)
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  phat <- y / n

  # Wald interval, truncated to the binomial parameter space.
  wald_half <- z * sqrt(phat * (1 - phat) / n)
  wald <- pmax(0, pmin(1, c(phat - wald_half, phat + wald_half)))

  # Wilson score interval obtained by inverting the score test.
  wilson_center <- (phat + z^2 / (2 * n)) / (1 + z^2 / n)
  wilson_half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) /
    (1 + z^2 / n)
  wilson <- c(
    max(0, wilson_center - wilson_half),
    min(1, wilson_center + wilson_half)
  )

  # Equal-tail Clopper--Pearson interval from beta quantiles.
  exact <- c(
    if (y == 0) 0 else stats::qbeta(alpha / 2, y, n - y + 1),
    if (y == n) 1 else stats::qbeta(1 - alpha / 2, y + 1, n - y)
  )

  # Likelihood-ratio interval: solve the deviance cutoff on either side of the MLE.
  cutoff <- stats::qchisq(conf.level, df = 1)
  objective <- function(p) binomial_lr_deviance(p, y, n) - cutoff
  eps <- sqrt(.Machine$double.eps)
  lr_lower <- if (y == 0) 0 else stats::uniroot(objective, c(eps, phat))$root
  lr_upper <- if (y == n) 1 else stats::uniroot(objective, c(phat, 1 - eps))$root

  intervals <- rbind(
    Wald = wald,
    Wilson = wilson,
    `Likelihood ratio` = c(lr_lower, lr_upper),
    `Clopper-Pearson` = exact
  )
  colnames(intervals) <- c("lower", "upper")
  intervals
}

binomial_interval_coverage <- function(n, p_grid, conf.level = 0.95) {
  interval_array <- array(
    NA_real_,
    dim = c(n + 1L, 4L, 2L),
    dimnames = list(y = 0:n, method = NULL, endpoint = c("lower", "upper"))
  )
  for (y in 0:n) interval_array[y + 1L, , ] <- binomial_intervals(y, n, conf.level)
  dimnames(interval_array)$method <- rownames(binomial_intervals(0, n, conf.level))

  coverage <- matrix(NA_real_, nrow = length(p_grid), ncol = 4L)
  mean_length <- matrix(NA_real_, nrow = length(p_grid), ncol = 4L)
  for (i in seq_along(p_grid)) {
    probabilities <- stats::dbinom(0:n, size = n, prob = p_grid[i])
    for (j in seq_len(4L)) {
      covered <- interval_array[, j, "lower"] <= p_grid[i] &
        p_grid[i] <= interval_array[, j, "upper"]
      coverage[i, j] <- sum(probabilities * covered)
      mean_length[i, j] <- sum(
        probabilities * (interval_array[, j, "upper"] - interval_array[, j, "lower"])
      )
    }
  }

  colnames(coverage) <- dimnames(interval_array)$method
  colnames(mean_length) <- dimnames(interval_array)$method
  list(p = p_grid, coverage = coverage, mean_length = mean_length)
}
