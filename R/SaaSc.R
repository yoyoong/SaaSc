


#' Title
#'
#' @param object
#' @param features
#' @param assay
#' @param slot
#' @param method
#' @param avc
#'
#' @return
#' @export
#'
#' @examples
rebuildMatrix <- function(object, features = NULL, assay = "RNA", slot = "data", method = NULL, avc = 0.8) {
  if (avc > 1 | avc < 0) {
    stop("The value of avc should be 0 to 1.")
  }

  if (is.null(method)) {
    matrix <- GetAssayData(object, assay = assay, slot = slot)
  } else if (method == "MCA") { # rebuild by MCA
    object <- NormalizeData(object)
    object <- RunMCA(object, nmcs = 50, features = features, assay = assay, slot = slot)
    loadings <- object@reductions$mca@feature.loadings
    embeddings <- object@reductions$mca@cell.embeddings
    stdev <- object@reductions$mca@stdev
    stdev.cumsum <- cumsum(stdev) / sum(stdev)
    topk <- max(min(which(stdev.cumsum >= avc)), 10)
    matrix <- loadings[, 1:topk] %*% t(embeddings[, 1:topk])
    matrix <- log(abs(matrix) + 1) * sign(matrix)
  }

  # Scale data
  matrix <- matrix - apply(matrix, 1, mean)
  matrix <- Matrix(data = matrix, sparse = TRUE)
  object[['SaaSc']] <- CreateAssayObject(data = matrix)
  return(object)
}


#' Get gene occurrence rate in background and signature genesets
#'
#' @param background.geneset Background genesets file name (in GMT format) or R list object.
#' @param signature.geneset Signature geneset, a R object of list which include tags (signature name) and values (gene list).
#' @param mode Signature geneset mode, "single" means each geneset is a signature, "multiple" means all genesets form a signature.
#'
#' @return A data frame include background and signatures gene rate.
#' @export
#'
#' @examples
getGeneRate <- function(background.geneset = NULL, signature.geneset = NULL, mode = "single") {
  # get background genesets
  if (class(background.geneset) == "list") {
    background_genesets <- background.geneset
  } else {
    if (!file.exists(background.geneset)) {
      stop(paste0("Cannot find file: ", background.geneset))
    }
    background_genesets <- readGMT(background.geneset)
  }

  # get signature geneset
  if (class(signature.geneset) != "list") {
    signature_geneset <- list(signature.geneset)
  } else {
    signature_geneset <- signature.geneset
  }
  if (is.null(names(signature_geneset))) {
    names(signature_geneset) <- paste0("signature", 1:length(signature_geneset))
  }

  # get all gene name
  all_genesets <- removeDuplicates(c(background_genesets, signature_geneset))
  all_gene <- unique(unlist(all_genesets, use.names = FALSE))

  background_matrix <- do.call('cbind', lapply(background_genesets, function(x) all_gene %in% x))
  background_col <- apply(background_matrix, 1, mean)
  signature_matrix <- do.call('cbind', lapply(signature_geneset, function(x) ifelse(all_gene %in% x, 1, 0)))
  if (mode == "single") {
    gene_rate = cbind(background_col, signature_matrix)
    colnames(gene_rate) <- c('background', colnames(signature_matrix))
  } else if (mode == "multiple") {
    signature_col <- apply(signature_matrix, 1, mean)
    gene_rate = cbind(background_col, signature_col)
    colnames(gene_rate) <- c('background', 'signature')
  }

  rownames(gene_rate) <- all_gene
  return (data.frame(gene_rate))
}

readGMT = function(gmt_file_path){
  chunks <- strsplit(readLines(gmt_file_path), '\t')
  genesets <- lapply(chunks, function(x) setdiff(x[3:length(x)], ''))
  names(genesets) <- sapply(chunks, function(x) x[1])
  return (genesets)
}

removeDuplicates <- function(lst) {
  unique_lst <- list()
  for (name in names(lst)) {
    if (!(name %in% names(unique_lst))) {
      unique_lst[[name]] <- lst[[name]]
    }
  }
  return(unique_lst)
}


#' Title
#'
#' @param object
#' @param model.file
#' @param celltype
#' @param cytokine
#' @param lambda
#' @param num.permutations
#' @param test.method
#'
#' @return
#' @export
#'
#' @examples
computeSignaling <- function(object, model.file = NULL, celltype = NULL, cytokine = NULL,
                             lambda = 10000, num.permutations = 1000, test.method = "two-sided") {
  if (is.null(model.file)){
    stop("Please input the model.file args")
  } else {
    if (!file.exists(model.file)) {
      stop(paste0("Cannot find file: ", model.file))
    }
  }
  model.data <- read.table(model.file, row.names = 1, header = TRUE, check.names = FALSE, sep = '\t')

  if (!is.null(cytokine)){
    model.cytokine <- colnames(model.data)
    if (all(cytokine %in% model.cytokine)) {
      use.cytokine <- cytokine
    } else if (any(cytokine %in% model.cytokine)) {
      use.cytokine <- intersect(cytokine, model.cytokine)
      unvalid.cytokine <- setdiff(cytokine, use.cytokine)
      message(paste("The following cytokine not exist in model file:", unvalid.cytokine))
    } else {
      stop("cytokine not exist in model file.")
    }
    model.data <- model.data[, use.cytokine]
  }

  object.celltype <- unique(object@meta.data$celltype)
  if (is.null(celltype)){ # support multiple celltype
    use.celltype <- object.celltype
  } else {
    use.celltype <- NULL
    if (all(celltype %in% object.celltype)) {
      use.celltype <- celltype
    } else if (any(celltype %in% object.celltype)) {
      use.celltype <- intersect(celltype, object.celltype)
      unvalid.celltype <- setdiff(celltype, use.celltype)
      message(paste("The following celltype not exist in Seurat object:", unvalid.celltype))
    } else {
      stop("celltype not exist in Seurat object.")
    }
  }

  use.gene <- intersect(rownames(model.data), rownames(object))
  use.cell = colnames(object[, object@meta.data$celltype %in% use.celltype])
  use.expr <- object[['SaaSc']]$data[use.gene, use.cell]
  use.model.data <- model.data[use.gene, ]

  cytokine.name <- colnames(model.data)
  cytokine.name.valid <- make.names(cytokine.name)

  # Build matrix for ridge regression
  X <- as.matrix(use.model.data)
  Y <- as.matrix(use.expr)
  # Calculate ridge regression coefficients
  tmp1 <- t(X) %*% X
  if (det(tmp1) == 0) {
    stop("t(X) %*% X == 0, can not do ridge regression...")
  }
  tmp2 <- solve(tmp1 + lambda * diag(ncol(X))) %*% t(X)
  beta <- tmp2 %*% Y

  if (!is.null(num.permutations) & num.permutations > 0) { # Do permutation
    step = max(1, floor(num.permutations / 10))
    average.matrix <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    average.sq.matrix <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    pvalue.matrix <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    for (i in 1:num.permutations) {
      if (i %% step == 0) {
        message(paste0("Process ", 100 * i / num.permutations, "%"))
      }

      Y.rand <- Y[sample(1:nrow(Y)), , drop = FALSE]
      beta.rand <- tmp2 %*% Y.rand

      if (test.method == "two-sided") {
        pvalue.matrix = pvalue.matrix + (beta >= abs(beta.rand))
      } else if (test.method == "greater") {
        pvalue.matrix = pvalue.matrix + (beta >= beta.rand)
      } else if (test.method == "less") {
        pvalue.matrix = pvalue.matrix + (beta <= beta.rand)
      }

      average.matrix <- average.matrix + beta.rand
      average.sq.matrix <- average.sq.matrix + beta.rand * beta.rand
    }

    average.matrix <- average.matrix / num.permutations
    average.sq.matrix <- average.sq.matrix / num.permutations
    pvalue.matrix <- pvalue.matrix / num.permutations
    std.matrix <- sqrt(average.sq.matrix - average.matrix * average.matrix)
    zscore.matrix <- (beta - average.matrix) / std.matrix
    zscore.matrix[is.na(zscore.matrix)] <- 0
    result <- t(zscore.matrix)
  } else {
    result <- t(beta)
  }

  return(result)
}



#' Title
#'
#' @param object
#' @param gene.rate
#' @param celltype
#' @param signature
#' @param cell.threshold
#'
#' @return
#' @export
#'
#' @examples
computeResponse <- function(object, gene.rate = NULL, celltype = NULL, signature = NULL, cell.threshold = 100) {
  if (is.null(gene.rate)){
    stop("Please input the gene.rate args")
  }

  if (is.null(object@meta.data$celltype)){
    stop("Celltype annotation not exist in Seurat object.")
  }
  object.celltype <- unique(object@meta.data$celltype)
  if (is.null(celltype)){
    stop("Please input the celltype args.")
  } else {
    if (!(celltype %in% object.celltype)) {
      stop("celltype not exist in Seurat object.")
    }
  }

  if (is.null(signature)) {
    signature <- paste0(setdiff(colnames(gene.rate), "background"), ".response")
  } else {
    if (length(celltype) != length(signature)) {
      stop("The number of celltypes not consistent with the signature number.")
    }
  }

  use.gene <- intersect(rownames(gene.rate), rownames(object))
  use.cell = colnames(object[, object@meta.data$celltype == celltype])
  if (length(use.cell) < cell.threshold) {
    # message(paste("Cell number of celltype", cell.type,  "less than 100, continue..."))
    # break
    stop("Cell number of celltype less than 100.")
  }
  use.expr <- object[['SaaSc']]$data[use.gene, use.cell]
  use.gene.rate <- gene.rate[use.gene, ]

  # Build matrix for regression
  X <- cbind(1, as.matrix(use.gene.rate))
  Y <- as.matrix(use.expr)

  # Do batch regression
  regression.result <- batchRegression(X, Y)

  result <- as.matrix(regression.result[3, ])
  colnames(result) <- signature
  return(result)
}



batchRegression <- function(X, Y) {
  tmp <- t(X) %*% X
  if (det(tmp) == 0) {
    stop("t(X) %*% X == 0, can not do regression...")
  }
  # Calculate regression coefficients
  beta <- solve(tmp) %*% t(X) %*% Y
  # Calculate residuals
  res <- Y - X %*% beta
  # Calculate the variance of the residuals
  dof <- nrow(Y) - ncol(X) + 1 # degree of freedom
  sigma_squared <- apply(res, MARGIN = 2, FUN = function(x) {
    return (sum(x ^ 2) / dof)
  })
  # Calculate the standard deviation of the coefficients
  XtX_inv <- solve(t(X) %*% X)
  se_beta <- sapply(sigma_squared, FUN = function(x) {
    return (sqrt(diag(x * XtX_inv)))
  })
  # Calculate the t-values
  t_values <- beta / se_beta
  return(t_values)
}


#' Title
#'
#' @param object
#' @param response.data
#' @param signaling.data
#' @param signature
#' @param cytokine
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
computeCorrelation <- function(object, response.data = NULL, signaling.data = NULL,
                               signature = NULL, cytokine = NULL, threshold = 100) {
  if (!setequal(rownames(response.data), rownames(signaling.data))) {
    stop("Cell names of response.data and signaling.data are not consistent.")
  }
  cell.name <- rownames(response.data)

  if (is.null(signature)){
    # stop("Please input the signature args.")
    valid.signature <- colnames(response.data)
  } else {
    valid.signature <- NULL
    if (all(signature %in% colnames(response.data))) {
      valid.signature <- signature
    } else if (any(signature %in% colnames(response.data))) {
      valid.signature <- intersect(signature, colnames(response.data))
      unvalid.signature <- setdiff(signature, valid.signature)
      message(paste0("The following signature not exist in response.data", unvalid.signature))
    } else {
      stop("signature not exist in response.data.")
    }
  }

  if (is.null(cytokine)){
    # stop("Please input the cytokine args.")
    valid.cytokine <- colnames(signaling.data)
  } else {
    valid.cytokine <- NULL
    if (all(cytokine %in% colnames(signaling.data))) {
      valid.cytokine <- cytokine
    } else if (any(cytokine %in% colnames(signaling.data))) {
      valid.cytokine <- intersect(cytokine, colnames(signaling.data))
      unvalid.cytokine <- setdiff(cytokine, valid.cytokine)
      message(paste0("The following cytokine not exist in signaling.data", unvalid.cytokine))
    } else {
      stop("cytokine not exist in signaling.data.")
    }
  }

  result.list = list()
  sample.names <- object@meta.data$sample[colnames(object) %in% cell.name]
  for (use.sample in unique(sample.names)) {
    use.sample.cell.name <- cell.name[sample.names %in% use.sample]
    nCell <- length(use.sample.cell.name)
    if (nCell >= threshold) {
      # message(paste("Processing sample:", use.sample))
      for (use.signature in valid.signature) {
        for (use.cytokine in valid.cytokine) {
              response.subset <- as.numeric(response.data[use.sample.cell.name, use.signature])
              signaling.subset <- as.numeric(signaling.data[use.sample.cell.name, use.cytokine])
              cor.test.result <- cor.test(response.subset, signaling.subset)
              new.list <- list(sample = use.sample, nCell = nCell,signature = use.signature, cytokine = use.cytokine,
                               correlation = cor.test.result$estimate[1], pvalue = cor.test.result$p.value)
              result.list[[length(result.list) + 1]] <- new.list

          }
      }
    } else {
      message(paste("Cell count of sample", use.sample, "less than threshold, continue..."))
    }
  }

  result <- data.frame(sample = sapply(result.list, function(x) x$sample),
                       nCell = sapply(result.list, function(x) x$nCell),
                       signature = sapply(result.list, function(x) x$signature),
                       cytokine = sapply(result.list, function(x) x$cytokine),
                       correlation = sapply(result.list, function(x) x$correlation),
                       pvalue = sapply(result.list, function(x) x$pvalue))
  return(result)
}


#' Title
#'
#' @param object
#' @param response.data
#' @param signaling.data
#' @param signature
#' @param cytokine
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
doInteraction <- function(object, response.data = NULL, signaling.data = NULL,
                          signature = NULL, cytokine = NULL, threshold = 100) {
  if (!setequal(rownames(response.data), rownames(signaling.data))) {
    stop("Cell names of response.data and signaling.data are not consistent.")
  }
  cell.name <- rownames(response.data)

  if (is.null(signature)){
    stop("Please input the signature args.")
  } else {
    valid.signature <- NULL
    if (all(signature %in% colnames(response.data))) {
      valid.signature <- signature
    } else if (any(signature %in% colnames(response.data))) {
      valid.signature <- intersect(signature, colnames(response.data))
      unvalid.signature <- setdiff(signature, valid.signature)
      message(paste0("The following signature not exist in response.data", unvalid.signature))
    } else {
      stop("signature not exist in response.data.")
    }
  }

  if (is.null(cytokine)){
    stop("Please input the cytokine args.")
  } else {
    valid.cytokine <- NULL
    if (all(cytokine %in% colnames(signaling.data))) {
      valid.cytokine <- cytokine
    } else if (any(cytokine %in% colnames(signaling.data))) {
      valid.cytokine <- intersect(cytokine, colnames(signaling.data))
      unvalid.cytokine <- setdiff(cytokine, valid.cytokine)
      message(paste0("The following cytokine not exist in signaling.data", unvalid.cytokine))
    } else {
      stop("cytokine not exist in signaling.data.")
    }
  }

  result.list = list()
  sample.names <- object@meta.data$sample[colnames(object) %in% cell.name]
  for (use.sample in unique(sample.names)) {
    use.sample.cell.name <- cell.name[sample.names %in% use.sample]
    nCell <- length(use.sample.cell.name)
    if (nCell >= threshold) {
      # message(paste("Processing sample:", use.sample))
      for (use.signature in valid.signature) {
        for (use.cytokine in valid.cytokine) {
          use.gene <- rownames(object[['SaaSc']]$data)
          response.subset <- as.numeric(response.data[use.sample.cell.name, use.signature])
          signaling.subset <- as.numeric(signaling.data[use.sample.cell.name, use.cytokine])
          expr.subset <- object[['SaaSc']]$data[use.gene, use.sample.cell.name]

          regression.result <- apply(expr.subset, 1, function(row) {
            signaling.expr <- signaling.subset * row
            X <- cbind(1, as.matrix(signaling.subset), as.matrix(row), as.matrix(signaling.expr))
            Y <- as.matrix(response.subset)
            lm.result <- lm(Y ~ X)
            tvalue <- summary(lm.result)$coefficients[4, 3]
            pvalue <- summary(lm.result)$coefficients[4, 4]
            return(c(tvalue, pvalue))
          })

          new.list <- lapply(use.gene, function(gene) {
            return(list(sample = use.sample, signature = use.signature, cytokine = use.cytokine, gene = gene,
                        interaction = regression.result[1, gene], pvalue = regression.result[2, gene]))
          })
          result.list <- c(result.list, new.list)
        }
      }
    } else {
      message(paste("Cell count of sample", use.sample, "less than threshold, continue..."))
    }
  }

  result <- data.frame(sample = sapply(result.list, function(x) x$sample),
                       signature = sapply(result.list, function(x) x$signature),
                       cytokine = sapply(result.list, function(x) x$cytokine),
                       gene = sapply(result.list, function(x) x$gene),
                       interaction = sapply(result.list, function(x) x$interaction),
                       pvalue = sapply(result.list, function(x) x$pvalue))
  return(result)
}


#' Title
#'
#' @param interaction.dataset
#' @param signature.cytokine
#' @param pvalue
#' @param method
#'
#' @return
#' @export
#'
#' @examples
getTresSignature <- function(interaction.dataset = NULL, signature.cytokine = NULL, pvalue = 0.05, method = "median") {
  if (is.null(interaction.dataset)) {
    stop("Please input the interaction.dataset args.")
  } else {
    if (class(interaction.dataset) == "character") {
      if (!file.exists(interaction.dataset)) {
        stop(paste0("Cannot find file: ", interaction.dataset))
      }
    } else if (class(interaction.dataset) == "list") {
      interaction.dataset <- do.call(rbind, interaction.dataset)
    }
  }

  all.cytokine <- unique(interaction.dataset$cytokine)
  if (is.null(signature.cytokine)){
    stop("Please input the signature.cytokine args.")
  } else {
    valid.cytokine <- NULL
    if (all(signature.cytokine %in% all.cytokine)) {
      valid.cytokine <- signature.cytokine
    } else if (any(cytokine %in% all.cytokine)) {
      valid.cytokine <- intersect(signature.cytokine, all.cytokine)
      unvalid.cytokine <- setdiff(signature.cytokine, valid.cytokine)
      message(paste0("The following signature.cytokine not exist in interaction.data", unvalid.cytokine))
    } else {
      stop("signature.cytokine not exist in interaction.data.")
    }
  }

  if (is.null(method)){
    stop("Please input the method args.")
  } else {
    if (method != "median" | method != "mean") {
      stop("Input method must be median or mean.")
    }
  }

  # if (is.null(filter)){
  #   stop("Please input the filter args.")
  # } else {
  #   if (filter != "all" | filter != "any") {
  #     stop("Input filter must be all or any.")
  #   }
  # }

  interaction.data <- interaction.dataset
  interaction.data <- interaction.data[interaction.data$cytokine == signature.cytokine, ] # filter cytokine
  interaction.data <- interaction.data[interaction.data$pvalue < pvalue, ] # filter pvalue

  interaction.data.grouped <- split(interaction.data, interaction.data$gene)
  if (method == "median") {
    result <- sapply(interaction.data.grouped, function(x) median(x$interaction))
  } else if (method == "mean") {
    result <- sapply(interaction.data.grouped, function(x) mean(x$interaction))
  }
  result <- data.frame(Tres.signature = result)

  return (result)
}


