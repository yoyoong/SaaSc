


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
    topk <- max(which(stdev.cumsum >= avc)[1], 10)
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
    # 先不支持多个celltype
    # use.celltype <- NULL
    # if (all(celltype %in% object.celltype)) {
    #   use.celltype <- celltype
    # } else if (any(celltype %in% object.celltype)) {
    #   use.celltype <- intersect(celltype, object.celltype)
    #   unvalid.celltype <- setdiff(celltype, use.celltype)
    #   message(paste0("The following celltype not exist in Seurat object", unvalid.celltype))
    # } else {
    #   stop("celltype not exist in Seurat object.")
    # }
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
  use.expr <- object[['SaaSc']]$data[use.gene, ]
  use.gene.rate <- gene.rate[use.gene, ]

  result <- matrix(0, nrow = length(use.cell), ncol = 1, dimnames = list(use.cell, signature))
  formula <- paste("expr", "~", paste(colnames(use.gene.rate), collapse = " + "))
  for (cell in use.cell) {
    use.cell.expr <- use.expr[, cell]
    lm_data <- data.frame(use.gene.rate, use.cell.expr)
    names(lm_data) <- c(colnames(use.gene.rate), "expr")
    lm <- lm(formula, data = lm_data)

    signature.tag <- setdiff(colnames(use.gene.rate), "background")
    tvalue <- summary(lm)$coefficients[signature.tag, 3]
    pvalue <- summary(lm)$coefficients[signature.tag, 4]

    result[cell, signature] <- tvalue
  }


  return(result)
}


#' Title
#'
#' @param object
#' @param model.file
#' @param celltype
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
computeSignaling <- function(object, model.file = NULL, celltype = NULL, lambda = 10000) {
  if (is.null(model.file)){
    stop("Please input the model.file args")
  } else {
    if (!file.exists(model.file)) {
      stop(paste0("Cannot find file: ", model.file))
    }
  }
  model.data <- read.table(model.file, row.names = 1, header = TRUE, check.names = FALSE, sep = '\t')

  if (is.null(celltype)){
    stop("Please input the celltype args.")
  }

  use.gene <- intersect(rownames(model.data), rownames(object))
  use.cell = colnames(object[, object@meta.data$celltype == celltype])
  use.expr <- object[['SaaSc']]$data[use.gene, ]
  use.model.data <- model.data[use.gene, ]

  cytokine.name <- colnames(model.data)
  result <- matrix(0, nrow = length(use.cell), ncol = ncol(use.model.data),
                             dimnames = list(use.cell, cytokine.name))
  cytokine.name.valid <- make.names(cytokine.name)
  formula <- paste("expr", "~", paste(cytokine.name.valid, collapse = " + "))
  for (cell in use.cell) {
    use.cell.expr <- use.expr[, cell]
    ridge.data <- data.frame(use.model.data, use.cell.expr)
    names(ridge.data) <- c(cytokine.name.valid, "expr")
    lm.ridge <- lm.ridge(formula, lambda = lambda, data = ridge.data, model = TRUE)

    result[cell, cytokine.name] <- lm.ridge$coef[cytokine.name.valid]
  }

  return(result)
}


computeCorrelation <- function(object, response.data = NULL, signaling.data = NULL,
                               signature = NULL, cytokine = NULL,  threshold = 50) {
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
  for (use.signature in valid.signature) {
    for (use.cytokine in valid.cytokine) {
      sample.names <- object@meta.data$sample[colnames(object) %in% cell.name]
      for (use.sample in unique(sample.names)) {
        use.sample.cell.name <- cell.name[sample.names %in% use.sample]
        nCell <- length(use.sample.cell.name)
        if (nCell >= threshold) {
          response.subset <- as.numeric(response.data[use.sample.cell.name, use.signature])
          signaling.subset <- as.numeric(signaling.data[use.sample.cell.name, use.cytokine])
          cor.test.result <- cor.test(response.subset, signaling.subset)
          new.list <- list(signature = use.signature, cytokine = use.cytokine, sample = use.sample, nCell = nCell,
                           correlation = cor.test.result$estimate[1], pvalue = cor.test.result$p.value)
          result.list[[length(result.list) + 1]] <- new.list
        }
      }
    }
  }

  result <- data.frame(signature = sapply(result.list, function(x) x$signature),
                    cytokine = sapply(result.list, function(x) x$cytokine),
                    sample = sapply(result.list, function(x) x$sample),
                    nCell = sapply(result.list, function(x) x$nCell),
                    correlation = sapply(result.list, function(x) x$correlation),
                    pvalue = sapply(result.list, function(x) x$pvalue))
  return(result)
}



