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

