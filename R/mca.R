

RunMCA <- function(X, nmcs = 50, features = NULL, reduction.name = "mca", slot = "data", assay = DefaultAssay(X), ...) {
  InitAssay <- DefaultAssay(X)
  DefaultAssay(X) <- assay
  data_matrix <- as.matrix(GetAssayData(X, slot = slot))
  MCA <- RunMCA.matrix(X = data_matrix, nmcs = nmcs, features = features)
  geneEmb <- MCA$featuresCoordinates
  cellEmb <- MCA$cellsCoordinates
  stdev <- MCA$stdev
  X <- setDimMCSlot.Seurat(X = X, cellEmb = cellEmb, geneEmb = geneEmb, stdev = stdev, reduction.name = reduction.name)
  DefaultAssay(X) <- InitAssay
  return(X)
}

RunMCA.matrix <- function(X, nmcs = 50, features = NULL, reduction.name = "MCA", ...) {
  # preprocessing matrix
  # ----------------------------------------------------
  if (!is.null(features)) {
    X <- X[features, ]
  }
  X <- as.matrix(X)
  X <- X[rowVars(X, useNames=FALSE) != 0, ]
  X <- X[str_length(rownames(X)) > 0, ]
  X <- X[!duplicated(rownames(X)), ]
  cellsN <- colnames(X)
  featuresN <- rownames(X)
  tic()
  message("Computing Fuzzy Matrix")
  MCAPrepRes <- MCAStep1(X)
  toc()
  message("Computing SVD")
  tic()
  SVD <- irlba::irlba(A = MCAPrepRes$Z, nv = nmcs + 1, nu = 1)[seq(3)]
  toc()
  message("Computing Coordinates")
  tic()
  MCA <- MCAStep2(Z = MCAPrepRes$Z, V = SVD$v[, -1], Dc = MCAPrepRes$Dc)
  component <- paste0(reduction.name, "_", seq(ncol(MCA$cellsCoordinates)))
  colnames(MCA$cellsCoordinates) <- component
  colnames(MCA$featuresCoordinates) <- component
  rownames(MCA$cellsCoordinates) <- cellsN
  rownames(MCA$featuresCoordinates) <- featuresN
  MCA$stdev <- SVD$d[-1]
  class(MCA) <- "MCA"
  toc()
  return(MCA)
}

setDimMCSlot.Seurat <- function(X, cellEmb, geneEmb, stdev = NULL, reduction.name = "mca", assay = DefaultAssay(X), ...) {
  colnames(cellEmb) <- paste0(reduction.name, "_", seq(ncol(cellEmb)))
  colnames(geneEmb) <- paste0(reduction.name, "_", seq(ncol(geneEmb)))
  DimReducObject <- CreateDimReducObject(embeddings = cellEmb, loadings = geneEmb, key = paste0(reduction.name, "_"), assay = assay)
  X@reductions[[reduction.name]] <- DimReducObject
  if (!is.null(stdev)) {
    X@reductions[[reduction.name]]@stdev <- sqrt(stdev)
  }
  X@reductions[[reduction.name]]@misc[["mca.flag"]] <- TRUE
  return(X)
}
