#' Modify fgsea::fgseaLabel to use user input fold change estimates
#' Original authors Gennady Korotkevich, Vladimir Sukhov, Nikolay Budin, Alexey Sergushichev
#' See https://bioconductor.org/packages/release/bioc/html/fgsea.html
#'
#' @param pathways List of gene sets to check.
#' @param dat Voom object with gene expression and sample metadata. Used in label randomization only'
#' @param label Numeric vector of labels for the correlation score of the same length as the number of columns in `mat`
#' @param nperm Number of permutations to do. Minimal possible nominal p-value is about 1/nperm
#' @param rand_var Character string specifying the variable to randomize in the label method
#' @param rand_est Data frame with random variable estimates. From prior run of BIGsea in the estimates slot
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of `gseaParam` before calculation of GSEA enrichment scores.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @param estimates Named vector of numeric fold change estimates
#' @param ... Additional parameters for kimma::kmFit needed if rand = "label"
#'
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#'
#' @export
#'
#' @importFrom Matrix invPerm
#' @importFrom fgsea calcGseaStat
#' @importFrom BiocParallel bplapply

fgseaLabel2 <- function(pathways, dat, label,
                        nperm, rand_var, rand_est,
                        minSize=1, maxSize=nrow(dat$E)-1,
                        nproc=0,
                        gseaParam=1,
                        BPPARAM=NULL,
                        estimates=NULL, ...) {
  g <- category <- subcategory <- gene <- variable <- estimate <- NULL

  # Define minimum pathway size at least 1
  minSize <- max(minSize, 1)

  # Determine pathways with genes of interest
  pathwaysFiltered <- lapply(pathways, function(p) {
    unique(stats::na.omit(fastmatch::fmatch(p, names(estimates))))
  })
  pathwaysSizes <- sapply(pathwaysFiltered, length)

  # Filter pathways by min and max sizes
  toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
  m <- length(toKeep)

  # Empty result if no genes in pass filter pathways
  if (m == 0) {
    return(tibble::tibble(
      group=g, gs_cat=category, gs_subcat=subcategory,
      pathway="No overlap of query genes and specified database."))
  }

  # Filter pathways to those with genes of interest
  pathwaysFiltered <- pathwaysFiltered[toKeep]
  pathwaysSizes <- pathwaysSizes[toKeep]

  # Order original model estimates
  ranksOrder <- order(estimates, decreasing=TRUE)
  ranksOrderInv <- Matrix::invPerm(ranksOrder)
  stats <- estimates[ranksOrder]

  pathwaysReordered <- lapply(pathwaysFiltered, function(x) ranksOrderInv[x])

  # Calculate GSEA for original model estiamtes
  gseaStatRes <- do.call(rbind,
                         lapply(pathwaysReordered, fgsea::calcGseaStat,
                                stats=stats,
                                returnLeadingEdge=TRUE))

  # Extract and format GSEA results
  leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
  pathwayScores <- unlist(gseaStatRes[, "res"])

  universe <- seq_along(stats)

  # Make random variables
  n_samples <- ncol(dat$E)
  rand_org <- as.numeric(as.factor(unlist(dat$targets[,label])))
  n_ones <- length(rand_org[rand_org==1])

  combinations <- utils::combn(n_samples, n_ones)

  # Convert combinations to assignments of 1 and 2
  rand_labels <- apply(combinations, 2, function(idx) {
    rand_labels <- rep(2, n_samples)  # start with all 2s
    rand_labels[idx] <- 1  # assign 1 to the positions in the combination
    return(rand_labels)
  })

  # Select total randomizations to run
  if(nperm > ncol(rand_labels)){
    nperm <- ncol(rand_labels)
    message(paste0("nperm is greater than total unique randomizations. Setting nperm to ", nperm, "."))
  } else{
    set.seed(42)
    rand_select <- sample(c(1:ncol(rand_labels)), size=nperm, replace=FALSE)
    rand_labels <- rand_labels[,rand_select]
  }

  # Setup parallel processing
  if(nperm>=100){ granularity <- 100 } else{ granularity <- 1 }
  permPerProc <- rep(granularity, floor(nperm / granularity))
  if (nperm - sum(permPerProc) > 0) {
    permPerProc <- c(permPerProc, nperm - sum(permPerProc))
  }
  set.seed(42)
  seeds <- sample.int(10^9, length(permPerProc))
  BPPARAM <- fgsea:::setUpBPPARAM(nproc, BPPARAM)

  # Run models for randomized labels
  if(is.null(rand_est)){
    message("Label randomization runs multiple kimma models. This may take a long time.\n")
    # Hold kimma results
    rand_est <- data.frame(matrix(nrow=nrow(dat$E), ncol=0))
    rownames(rand_est) <- rownames(dat$E)

    for(i in 1:ncol(rand_labels)){
      print(paste("kimma permutation",i))
      #Extract group labels
      dat$targets$random <- as.factor(rand_labels[,i])

      # parse kimma parameters in ...
      dots <- list(...)
      if(!"kin" %in% names(dots)){dots$kin <- NULL}
      if(!"patientID" %in% names(dots)){dots$patientID <- "ptID"}
      if(!"libraryID" %in% names(dots)){dots$libraryID <- "libID"}
      if(!"model" %in% names(dots)){stop("model required for label randomization.")}
      if(!"use_weights" %in% names(dots)){dots$use_weights <- TRUE}
      if(!"run_lm" %in% names(dots)){dots$run_lm <- FALSE}
      if(!"run_lme" %in% names(dots)){dots$run_lme <- FALSE}
      if(!"run_lmerel" %in% names(dots)){dots$run_lmerel <- FALSE}
      if(!"p_method" %in% names(dots)){dots$p_method <- "BH"}

      # reformat model to use random variable
      model_rand <- gsub(rand_var,"random",dots$model)
      # Run kimma model to calculate estimates
      suppressMessages(fit <- kimma::kmFit(dat = dat,
                                           kin = dots$kin,
                                           patientID = dots$patientID,
                                           libraryID = dots$libraryID,
                                           model=model_rand,
                                           use_weights = dots$use_weights,
                                           run_lm = dots$run_lm,
                                           run_lme = dots$run_lme,
                                           run_lmerel = dots$run_lmerel,
                                           metrics = FALSE,
                                           run_contrast = FALSE,
                                           contrast_var = NULL,
                                           processors = nproc,
                                           p_method = dots$p_method
      ))

      #Extract estimates
      rand_est[,i] <- fit[[1]] %>%
        dplyr::arrange(match(gene, rownames(rand_est))) %>%
        dplyr::filter(variable=="random2") %>%
        dplyr::pull(estimate)
    }
  } else{
    rand_est <- rand_est
  }

  # Run GSEA
  counts <- BiocParallel::bplapply(seq_along(permPerProc), function(j) {
    nperm1 <- permPerProc[j]

    randEsPs <- lapply(seq_len(nperm1), function(k) {
      rank_est_temp <- rand_est[, k]
      rank_ord_temp <- sort.list(rank_est_temp, decreasing=TRUE)
      gene_temp <- Matrix::invPerm(rank_ord_temp)
      stats_temp <- rank_est_temp[rank_ord_temp]

      randEsP <- fgsea:::calcGseaStatBatch(
        stats = stats_temp,
        selectedStats = pathwaysFiltered,
        geneRanks = gene_temp,
        gseaParam = gseaParam)
    })

    randEsPs <- do.call(cbind, randEsPs)

    leEs <- apply(sweep(randEsPs, MARGIN = 1, pathwayScores, `<=`), 1, sum)
    geEs <- apply(sweep(randEsPs, MARGIN = 1, pathwayScores, `>=`), 1, sum)
    leZero <- apply(randEsPs <= 0, 1, sum)
    geZero <- apply(randEsPs >= 0, 1, sum)
    leZeroSum <- apply(pmin(randEsPs, 0), 1, sum)
    geZeroSum <- apply(pmax(randEsPs, 0), 1, sum)

    data.table::data.table(pathway=seq_len(m),
                           leEs=leEs, geEs=geEs,
                           leZero=leZero, geZero=geZero,
                           leZeroSum=leZeroSum, geZeroSum=geZeroSum
    )
  }, BPPARAM=BPPARAM)

  counts <- data.table::rbindlist(counts)

  # Getting rid of check NOTEs
  leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
  pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
  nMoreExtreme=nGeEs=nLeEs=size=NULL
  leadingEdge=NULL
  .="damn notes"

  pvals <- counts %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(pval = min((1+sum(leEs)) / (1 + sum(leZero)),
                                (1+sum(geEs)) / (1 + sum(geZero))),
                     leZeroMean = sum(leZeroSum) / sum(leZero),
                     geZeroMean = sum(geZeroSum) / sum(geZero),
                     nLeEs=sum(leEs),
                     nGeEs=sum(geEs)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(padj = stats::p.adjust(pval, method="BH"),
                  ES = pathwayScores[pathway],
                  NES = ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean)),
                  nMoreExtreme = ifelse(ES > 0, nGeEs, nLeEs),
                  size = pathwaysSizes,
                  pathway= names(pathwaysFiltered),
                  leadingEdge=leadingEdges) %>%
    dplyr::select(-leZeroMean, -geZeroMean, -nLeEs, -nGeEs)

  result <- list()
  result[["pvals"]] <- pvals
  result[["estimates"]] <- rand_est
  return(result)
}

