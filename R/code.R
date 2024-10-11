doBinarization <- function(score.data) {
  # fit bimodal distribution
  fit <- normalmixEM(score.data[[1]], k = 2, maxit = 1000, eps = 1e-6)
  mean1 <- fit$mu[1]
  mean2 <- fit$mu[2]
  sd1 <- fit$sigma[1]
  sd2 <- fit$sigma[2]

  # compute the intersection of two normal distributions
  x_range <- seq(mean1, mean2, length.out = 1000)
  norm1 <- function(x) dnorm(x, mean = mean1, sd = sd1)
  norm2 <- function(x) dnorm(x, mean = mean2, sd = sd2)
  y1 <- norm1(x_range)
  y2 <- norm2(x_range)
  diff <- abs(y1 - y2)
  idx <- which.min(diff)
  intersection <- x_range[idx]

  threshold <- (mean1 + mean2) / 2
  score.data$label <- ifelse(score.data[[1]] > threshold, "positive", "negative")
  return(score.data)
}
