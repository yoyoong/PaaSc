doBinarization2 <- function(score.data) {
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

doNormalization <- function(object, score.data, background.file = NULL, pathway.geneset = NULL) {
  background_genesets <- readGMT(background.file)
  all_gene <- unique(unlist(background_genesets, use.names = FALSE))
  repeat_num <- 100

  for (pathway in names(pathway.geneset)) {
    random_gene_num <- length(unlist(pathway.geneset))
    random_score_data <- matrix(0, nrow = nrow(score.data), ncol = repeat_num,
                                dimnames = list(rownames(score.data), 1:repeat_num))
    for (i in 1:repeat_num) {
      random_geneset <- sample(all_gene, random_gene_num)
      random_geneset = list(random_geneset = random_geneset)
      gene_rate <- getGeneRate(background.file, pathway.geneset = random_geneset)
      regression_data <- doRegression(object, gene.rate = gene_rate)
      score_data <- computeScore(object, regression.data = regression_data, pvalue = 0.05, weight = FALSE)
      random_score_data[, i] <- score_data[[1]]
    }

    random_score_data_mean <- apply(random_score_data, 1, mean)
    random_score_data_sd <- apply(random_score_data, 1, sd)

    score.data[pathway] <- (score.data[pathway] - random_score_data_mean) / random_score_data_sd
    message(paste(pathway, "process end."))
  }

  return(score.data)
}
