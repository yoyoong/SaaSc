


#' Rebuild the gene expression matrix
#'
#' @param object A seurat object, must include the celltype in meta.data assay.
#' @param features Character vector of feature names. If not specified all features will be taken.
#' @param assay Name of assay MCA is being run on, default set to "RNA".
#' @param slot The slot to pull expression data for Seurat, default set to "data".
#' @param method Dimensionality reduction and reconstruction method. The choices are NULL(do not rebuild), "mca" , "nmf" and "pca".
#' @param dim The number of selected dimensions, 1 <= dim <= 50. Enter only one of avc and dim.
#' @param avc Accumulative variance contribution threshold, 0 <= avc <= 1. Enter only one of avc and dim.
#'
#' @return A new seurat object include the 'rebuild.RNA' assay that store the reconstructed matrix.
#' @export
#'
#' @examples
rebuildMatrix <- function(object, features = NULL, assay = "RNA", slot = "data", method = NULL, dim = NULL, avc = NULL) {
  if (!is.null(avc) & !is.null(dim)) {
    stop("Cannot enter avc and dim at the same time.")
  }
  if (is.null(avc) & is.null(dim)) {
    stop("Cannot enter avc and dim at the same time.")
  }
  if (!is.null(avc)) {
    if (avc > 1.0 | avc < 0.0) {
      stop("The value of avc should be 0 to 1.")
    }
  }

  object <- NormalizeData(object)
  if (!is.null(methods)) {
    for (use.method in method){
      if (use.method == "mca") { # rebuild by MCA
        object <- RunMCA(object, nmcs = 50, features = features, assay = assay, slot = slot)
      } else if (use.method == "nmf") { # rebuild by NMF
        if (!is.null(features)) {
          object <- object[features, ]
        }
        nmf.object <- nmf(object[[assay]][slot], k = 50, tol = 1e-05, maxit = 500)
        loadings <- nmf.object$w
        rownames(loadings) <- rownames(object)
        embeddings <- t(nmf.object$h)
        rownames(embeddings) <- Cells(object)
        stdev <- nmf.object$d

        DimReducObject <- CreateDimReducObject(embeddings = embeddings, loadings = loadings,
                                               key = paste0(use.method, "_"), assay = assay)
        object@reductions[[use.method]] <- DimReducObject
        object@reductions[[use.method]]@stdev <- stdev
        object@reductions[[use.method]]@misc[[paste0(use.method, ".flag")]] <- TRUE
      } else if (use.method == "pca") { # rebuild by PCA
        if (is.null(features)) {
          features <- rownames(object)
        }
        object <- RunPCA(object, npcs = 50, features = features, assay = assay)
      } else {
        print(paste0("The value of methods include invalid item: ", method))
      }
    }

    object@reductions[['active.method']] = method
    object@reductions[['active.dim']] = dim
    object@reductions[['active.avc']] = avc
  }

  return(object)
}


#' Get gene occurrence rate in background and signature genesets
#'
#' @param background.geneset Background genesets file name (in GMT format) or R list object.
#' @param signature.geneset Signature geneset, a R object of list which include tags (signature name) and values (gene list).
#' @param mode Signature geneset mode, "single" means each geneset is a signature, "multiple" means all genesets form a signature.
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


#' Get rebuilt matrix
#'
#' @param object A Seurat object processed by rebuildMatrix().
#'
#' @return A matrix of rebuilt matrix by reconstruction method.
#' @export
#'
#' @examples
getRebuildMatrix <- function(object) {
  active.method = object@reductions[['active.method']]

  if (is.null(active.method)) {
    rebuild.GEM = GetAssayData(object, assay = "RNA", slot = "data")
  } else {
    dim = object@reductions[['active.dim']]
    avc = object@reductions[['active.avc']]
    if (is.null(object@reductions[[active.method]])) {
      stop(paste("Cannot find", active.method, "reductions object in object. Run rebuildMatrix() at first."))
    }
    loadings <- object@reductions[[active.method]]@feature.loadings
    embeddings <- object@reductions[[active.method]]@cell.embeddings
    stdev <- object@reductions[[active.method]]@stdev
    stdev.cumsum <- cumsum(stdev) / sum(stdev)
    if (!is.null(avc)) {
      topk <- min(which(stdev.cumsum >= avc))
    } else if (!is.null(dim)) {
      topk <- dim
    } else {
      topk <- 50
    }
    rebuild.GEM <- loadings[, 1:topk] %*% t(embeddings[, 1:topk])
    rownames(rebuild.GEM) <- rownames(loadings)
    colnames(rebuild.GEM) <- rownames(embeddings)
    rebuild.GEM <- log(abs(rebuild.GEM) + 1) * sign(rebuild.GEM)
  }

  rebuild.GEM <- rebuild.GEM - apply(rebuild.GEM, 1, mean)
  return(rebuild.GEM)
}


#' Compute signaling
#'
#' @param object A Seurat object processed by rebuildMatrix().
#' @param model.file Cytokine effects on genes quantitative matrix file, which row is gene and column is cytokine.
#' @param celltype Celltypes to be calculated.
#' @param cytokine Cytokines to be calculated.
#' @param lambda Lambda values of ridge regression, default set to 10000.
#' @param num.permutations Number of permutations to use, if greater than 0, compute empirical t-values using a permutation test. Default set to 1000.
#' @param test.method Significance test method. The choices are "two-sided", "greater" and "less", default is "two-sided".
#'
#' @return A data frame store the cytokine signaling of single cells.
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

  if (is.null(celltype)){ # support multiple celltype
    use.cell = colnames(object)
  } else {
    object.celltype <- unique(object@meta.data$celltype)
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
    use.cell = colnames(object[, object@meta.data$celltype %in% use.celltype])
  }

  rebuild.GEM = getRebuildMatrix(object)
  use.gene <- intersect(rownames(model.data), rownames(rebuild.GEM))
  use.expr <- rebuild.GEM[use.gene, use.cell]
  use.model.data <- model.data[use.gene, ]

  cytokine.name <- colnames(model.data)
  cytokine.name.valid <- make.names(cytokine.name)

  # Build matrix for ridge regression
  X <- as.matrix(use.model.data)
  Y <- as.matrix(use.expr)

  # Calculate ridge regression coefficients
  result <- ridgeRegression(X, Y, scale = TRUE, lambda = lambda, num.permutations = num.permutations, test.method = test.method)

  return(result)
}


ridgeRegression <- function(X, Y, scale = TRUE, lambda = 0, num.permutations = NULL, test.method = "two-sided") {
  if (scale) {
    X <- scale(X)
    Y <- scale(Y)
  }

  tmp1 <- crossprod(X)
  # if (det(tmp1) == 0) {
  #   stop("t(X) %*% X == 0, can not do regression...")
  # }
  # Calculate regression coefficients
  tmp2 <- solve(tmp1 + lambda * diag(ncol(X))) %*% t(X)
  beta <- tmp2 %*% Y

  if (is.null(num.permutations) || num.permutations <= 0) {
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
    result <- t(t_values)
  } else { # Do permutation
    step = max(1, floor(num.permutations / 10))
    average.matrix <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    average.sq.matrix <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    pvalue.matrix <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
    for (i in 1:num.permutations) {
      if (i %% step == 0) {
        message(paste0("Process ", 100 * i / num.permutations, "%"))
      }

      Y.rand <- Y[sample(1:nrow(Y)), ]
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
  }

  return(result)
}


#' Compute response
#'
#' @param object A Seurat object processed by rebuildMatrix().
#' @param gene.rate Gene rate matrix generated by getGeneRate().
#' @param celltype A cell type to be calculated.
#' @param signature Signature name corresponding to the cell type.
#'
#' @return A data frame store the signature response of single cells.
#' @export
#'
#' @examples
computeResponse <- function(object, gene.rate = NULL, celltype = NULL, signature = NULL) {
  if (is.null(gene.rate)){
    stop("Please input the gene.rate args")
  }

  if (is.null(celltype)){
    use.cell = colnames(object)
  } else {
    if (!(celltype %in% object@meta.data$celltype)) {
      stop("celltype not exist in Seurat object.")
    }
    use.cell = colnames(object[, object@meta.data$celltype == celltype])
  }

  rebuild.GEM = getRebuildMatrix(object)
  use.gene <- intersect(rownames(gene.rate), rownames(rebuild.GEM))
  use.expr <- rebuild.GEM[use.gene, use.cell]
  use.gene.rate <- gene.rate[use.gene, ]

  # Build matrix for regression
  X <- as.matrix(use.gene.rate)
  Y <- as.matrix(use.expr)

  # Do batch regression
  regression.result <- ridgeRegression(X, Y, scale = TRUE, lambda = 0, num.permutations = 0)
  result <- as.matrix(regression.result[, 2])
  colnames(result) <- signature
  return(result)
}


#' Compute correlation between signaling and response
#'
#' @param object A Seurat object processed by rebuildMatrix, and include the sample in meta.data assay.
#' @param response.data Response data frame generated by computeResponse().
#' @param signaling.data Signaling data frame generated by computeSignaling().
#' @param signature Signature to be calculated.
#' @param cytokine Cytokine to be calculated.
#' @param threshold The cell number threshold for samples. Default set to 100.
#'
#' @return A data frame include 6 columns: sample, nCell, signature, cytokine, correlation, qvalue.
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


#' Do interaction of Tres model
#'
#' @param object A Seurat object processed by rebuildMatrix().
#' @param response.data Response data frame generated by computeResponse().
#' @param signaling.data Signaling data frame generated by computeSignaling().
#' @param signature Signature to be calculated.
#' @param cytokine Cytokine to be calculated.
#' @param threshold The cell number threshold for samples. Default set to 100.
#'
#' @return A data frame include 6 columns: sample, nCell, signature, gene, t, qvalue.
#' @export
#'
#' @examples
doInteraction <- function(object, response.data = NULL, signaling.data = NULL, signature = NULL, cytokine = NULL, threshold = 100) {
  if (!setequal(rownames(response.data), rownames(signaling.data))) {
    stop("Cell names of response.data and signaling.data are not consistent.")
  }
  cell.name <- intersect(rownames(response.data), colnames(object[["RNA"]]$data))

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
  rebuild.GEM = getRebuildMatrix(object)
  sample.names <- object@meta.data$sample[colnames(object) %in% cell.name]
  for (use.sample in unique(sample.names)) {
    use.sample.cell.name <- cell.name[sample.names %in% use.sample]
    nCell <- length(use.sample.cell.name)
    if (nCell >= threshold) {
      for (use.signature in valid.signature) {
        for (use.cytokine in valid.cytokine) {
          use.gene <- rownames(rebuild.GEM)
          response.subset <- response.data[use.sample.cell.name, use.signature]
          signaling.subset <- signaling.data[use.sample.cell.name, use.cytokine]
          expr.subset <- rebuild.GEM[use.gene, use.sample.cell.name]

          regression.result <- apply(expr.subset, 1, function(row) {
            signaling.expr <- signaling.subset * row
            X <- cbind(1, as.matrix(signaling.subset), as.matrix(row), as.matrix(signaling.expr))
            colnames(X) <- c('const', 'signaling', 'expr', 'interaction')
            Y <- as.matrix(response.subset)
            tryCatch({
              lm.result <- lm(Y ~ X)
              tvalue <- summary(lm.result)$coefficients['Xinteraction', 3]
              pvalue <- summary(lm.result)$coefficients['Xinteraction', 4]
              return(c(tvalue, pvalue))
            }, error = function(e) {
              tvalue <- NA
              pvalue <- NA
              return(c(tvalue, pvalue))
            })
          })
          tvalue <- regression.result[1, ]
          pvalue <- regression.result[2, ]
          qvalue <- p.adjust(pvalue, method = "BH")

          new.list <- lapply(use.gene, function(gene) {
            return(list(sample = use.sample, signature = use.signature, cytokine = use.cytokine, gene = gene,
                        t = tvalue[gene], pvalue = pvalue[gene], qvalue = qvalue[gene]))
          })
          result.list <- c(result.list, new.list)
          # message(paste("Process sample:", use.sample, ", signature:", use.signature, ", cytokine:", use.cytokine, "end."))
        }
        # message(paste("Process sample:", use.sample, ", signature:", use.signature, "end."))
      }
      message(paste("Process sample:", use.sample, "end."))
    } else {
      message(paste("Cell count of sample", use.sample, "less than threshold, continue..."))
    }
  }

  result <- data.frame(sample = sapply(result.list, function(x) x$sample),
                       signature = sapply(result.list, function(x) x$signature),
                       cytokine = sapply(result.list, function(x) x$cytokine),
                       gene = sapply(result.list, function(x) x$gene),
                       t = sapply(result.list, function(x) x$t),
                       qvalue = sapply(result.list, function(x) x$qvalue))
  return(result)
}


#' Get Tres signature
#'
#' @param interaction.dataset Interaction dataset of all training samples, which generated by computeResponse and aggregated.
#' @param signature.cytokine Cytokine to be calculated.
#' @param method Method for summarizing Tres signatures when the number of signature.cytokine more than one. The choices include "medium" and "mean", default set to "medium."
#' @param qvalue q value threshold, default is 0.05.
#' @param cutoff Sample proportion threshold less than q value, default set to 0.5.
#'
#' @return Tres signature data frame include Tres score of valid genes.
#' @export
#'
#' @examples
getTresSignature <- function(interaction.dataset = NULL, signature.cytokine = NULL,
                             qvalue = 0.05, cutoff = 0.5, method = "median") {
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
    if (method != "median" & method != "mean") {
      stop("Input method must be median or mean.")
    }
  }

  if (!is.null(valid.cytokine)) { # filter cytokine
    interaction.dataset <- interaction.dataset[interaction.dataset$cytokine %in% valid.cytokine, ]
  }


  if (!is.null(cutoff)){
    # group by gene and calculate the total number of samples for each gene
    gene_sample_counts <- aggregate(sample ~ gene, data = interaction.dataset, FUN = length)
    colnames(gene_sample_counts) <- c("gene", "total_samples")
    # filter qvalue and calculate the filter number of samples for each gene
    filtered_samples <- interaction.dataset[interaction.dataset$qvalue <= qvalue, ]
    gene_filtered_counts <- aggregate(sample ~ gene, data = filtered_samples, FUN = length)
    colnames(gene_filtered_counts) <- c("gene", "filtered_samples")
    # merge the total sample and filter number of samples
    merged_data <- merge(gene_sample_counts, gene_filtered_counts, by = "gene", all = TRUE)
    merged_data$filtered_samples[is.na(merged_data$filtered_samples)] <- 0
    merged_data$proportion <- merged_data$filtered_samples / merged_data$total_samples
    # get the gene fit the cutoff
    filtered_gene <- subset(merged_data, proportion >= cutoff)$gene

    interaction.dataset <- interaction.dataset[interaction.dataset$gene %in% filtered_gene, ]
  }
  interaction.data <- interaction.dataset
  # interaction.data <- interaction.dataset[interaction.dataset$qvalue <= qvalue, ]

  if (method == "median") {
    aggregated_by_cytokine <- aggregate(t ~ gene + sample, data = interaction.data, FUN = median)
    aggregated_by_sample <- aggregate(t ~ gene, data = aggregated_by_cytokine, FUN = median)
    result <- aggregated_by_sample
  } else if (method == "mean") {
    aggregated_by_cytokine <- aggregate(t ~ gene + sample, data = interaction.data, FUN = mean)
    aggregated_by_sample <- aggregate(t ~ gene, data = aggregated_by_cytokine, FUN = mean)
    result <- aggregated_by_sample
  }

  colnames(result) <- c("gene", "Tres.score")
  return (result)
}


