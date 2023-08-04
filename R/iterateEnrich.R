#' Iterative Hypergeometric Enrichment
#' @description
#' Iteratively run flexEnrich on a list of features containing multiple annotations. Will select a random annotation for a given feature and run enrichment n times. Returns p-values and k/K ratios for each iteration, as well as a summative data frame with median p-values, median k/K ratios, and FDR-correction of the median p-values.
#'
#'
#' @param anno_df Data frame containing features (such as methylation probes, SNP IDs, ATAC-seq peaks) and corresponding gene IDs for those features. Features with multiple annotations (ie rows) are iteratively sampled to 1 annotation per feature for enrichment
#' @param anno_featCol Character string. Name of the column containing the feature IDs (e.g. methylation probes, SNP IDs, ATAC-seq peaks). Default is the first column of anno_df
#' @param anno_annotationCol Character string. Name of the column containing the gene IDs (e.g. HGNC symbols, Entrez IDs, Ensembl IDs). Default is the second column of anno_df
#' @param niter Integer. Number of iterations. Default is 100.
#' @param p_adjust Character string. Method for p-value adjustment from p.adjust.methods(). Default is "fdr"
#' @param ID Character string. Type of gene annotation used in anno_df. One of SYMBOL, ENTREZ, ENSEMBL. Default is "SYMBOL"
#' @param species Character string denoting species of interest. "human" or "mouse." Default is "human"
#' @param category Character string denoting Broad gene set database
#' @param subcategory Character string denoting Broad gene set sub-database
#' \tabular{rrrrr}{
#'  \strong{category} \tab    \strong{subcategory}\cr
#'  C1  \tab      \cr
#'  C2  \tab      CGP\cr
#'  C2  \tab      CP\cr
#'  C2  \tab      CP:BIOCARTA\cr
#'  C2  \tab      CP:KEGG\cr
#'  C2  \tab      CP:PID\cr
#'  C2  \tab      CP:REACTOME\cr
#'  C2  \tab      CP:WIKIPATHWAYS\cr
#'  C3  \tab      \cr
#'  C5  \tab      GO:BP\cr
#'  C5  \tab      GO:CC\cr
#'  C5  \tab      GO:MF\cr
#'  C5  \tab      HPO\cr
#'  C6  \tab      \cr
#'  H   \tab      \cr
#'  }
#' @param db If not using Broad databases, a data frame with gene ontology including gene set name (column 1: gs_name) and gene ID (column2: gene_symbol, entrez_gene, or ensembl_gene as matches your anno_df annotations)
#' @param custom_bg Custom background. Formatted as a vector of gene IDs.
#' @param protein_coding Logical. Do you want to limit the background to only protein-coding genes? Default is TRUE
#' @param minOverlap Minimum overlap between a gene set and your list of query genes for hypergeometric enrichment to be calculated. Default is 1. For the iterative function, the only valid value is 1 at the moment.
#' @param minGeneSetSize Maximum overlap between a gene set and your list of query genes for hypergeometric enrichment to be calculated. Default is 10.
#' @param maxGeneSetSize Maximum size of a reference gene set for hypergeometric enrichment to be calculated. Default is 1e10
#' @param print_genes TRUE or FALSE. Do you want the results to include a list of genes that overlap between any given gene set and your query genes. Default is TRUE. Leaving this parameter as TRUE will make the function run slowly for very large datasets.
#' @param ncores Number of cores for parallel processing. Default is 1
#' @author Madison Cox
#'
#' @return List of data frames including
#'    - summary: median and min/max p-values, median k/K values and adjusted p-values. Note that medians treat iterations with no overlap in a gene set as p-value = 1 and k/K = 0
#'    - p_iterations: p-values for each iteration,
#'    - k/K_iterations:  k/K values for each iteration
#'
#' @importFrom foreach %dopar%
#' @importFrom stats median
#' @export
#'
#' @examples
#' df <- data.frame(feat = paste0("probe", rep(1:25, 4)),
#'                  annotation = names(example.gene.list[[1]]))
#'
#' iterateEnrich(anno_df = df,
#'               anno_featCol = "feat",
#'               anno_annotationCol = "annotation",
#'               niter = 5, ID = "ENSEMBL", category = "H")

iterateEnrich <- function(anno_df = NULL,
                          anno_featCol = NULL,
                          anno_annotationCol = NULL,
                          niter = 100,
                          p_adjust = "fdr",
                          ID = "SYMBOL",
                          species = "human",
                          category = NULL,
                          subcategory = NULL,
                          db = NULL,
                          custom_bg = NULL,
                          protein_coding = TRUE,
                          minOverlap = 1,
                          minGeneSetSize = 10,
                          maxGeneSetSize = 10000,
                          print_genes = TRUE,
                          ncores = 1){
  gs_cat <- gs_subcat <- pathway <- `k/K` <- pvalue <- genes <- max_pval <- median_kK <- median_pval <- min_pval <- NULL


  if(minOverlap > 1) {
    stop("Sorry, at this time iterative p-values can only be generated for a minOverlap of 1.")
  }
  #Set colnames if not provided
  if(is.null(anno_featCol)) { anno_featCol <- colnames(anno_df[1])}
  if(is.null(anno_annotationCol)) { anno_annotationCol <- colnames(anno_df[2])}

  ###### Parallel ######
  #setup parallel processors
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    #Use 2 in CRAN/Travis/AppVeyor
    processors.to.use <- 2
  } else if (is.null(ncores)){
    #Use 2 less than total if not user defined
    processors.to.use <- parallel::detectCores()-2
    if(processors.to.use == 0){
      stop("Error processors: Default resulted in 0. Please correct.")}
  } else {
    #Use user defined number
    processors.to.use <- ncores
  }

  cl <- parallel::makeCluster(processors.to.use)

  ### get pathway names for base data frame ###
  if(!is.null(category)){
    db.format <- msigdbr::msigdbr(species, category)
    #Subset subcategory if selected
    if(!is.null(subcategory)){
      if(length(db.format$gs_subcat[which(db.format$gs_subcat == subcategory)]) == 0){
        stop("Entered subcategory does not exist in MSigDB")
      }
      else{
        db.format <- db.format %>%
          dplyr::filter(grepl(paste0("^",subcategory), gs_subcat))
      }
    }
  } else if(!is.null(db)){
    db.format <- db %>%
      dplyr::mutate(gs_cat = "custom", gs_subcat=NA)
    #Name columns to match mgsibdbr
    colnames(db.format)[1] <- "gs_name"
    if(ID == "SYMBOL"){
      colnames(db.format)[2] <- "gene_symbol"
    } else if(ID == "ENSEMBL"){
      colnames(db.format)[2] <- "ensembl_gene"
    } else if(ID == "ENTREZ"){
      colnames(db.format)[2] <- "entrez_gene"
    } else{
      stop("Please like ID from SYMBOL, ENSEMBL, or ENTREZ.")
    }
  } else {
    stop("Please provide gene set information as Broad category/subcategory or in a data frame as db.")
  }

  base_df <- data.frame("pathway" = unique(db.format$gs_name))
  ### iterate enrichment with individual annotations ###

  parallel::clusterExport(cl, c("ensembl.human.db.pc", "ensembl.human.db.full", "ensembl.mouse.db.pc",
                                "entrez.human.db.pc", "entrez.human.db.full", "entrez.mouse.db.pc",
                                "symbol.human.db.pc", "symbol.human.db.full", "symbol.mouse.db.pc"),
                          envir = environment())
  doParallel::registerDoParallel(cl)

  iter_list <- foreach::foreach(i = 1:niter,
                                .packages = c("dplyr", "doParallel","msigdbr","stats","tibble", "foreach"),
                                .export = c("flexEnrich"), .noexport = c("ensembl.human.db.pc",
                                                                         "ensembl.human.db.full",
                                                                         "ensembl.mouse.db.pc",

                                                                         "entrez.human.db.pc",
                                                                         "entrez.human.db.full",
                                                                         "entrez.mouse.db.pc",

                                                                         "symbol.human.db.pc",
                                                                         "symbol.human.db.full",
                                                                         "symbol.mouse.db.pc")

  ) %dopar% {


    gl <- list()
    gv <- c()
    ### Get single annotation for each feature ###
    for(feat in unique(unlist(anno_df[,anno_featCol]))){
      syms <- anno_df
      syms <- syms[,anno_annotationCol][which(syms[,anno_featCol] == feat)]
      syms <- unlist(syms)

      set.seed(42)
      sym.pick <- sample(syms,1) # choose random symbol from list of options

      gv <- c(gv, sym.pick)
    }

    gl[["genes"]] <- gv

    ### run enrichemt on reduced gene list ###
    prof <- flexEnrich(gene_list = gl,
                       ID = ID,
                       species = species,
                       category = category,
                       subcategory = subcategory,
                       db = db,
                       custom_bg = custom_bg,
                       protein_coding = protein_coding,
                       minOverlap = minOverlap,
                       minGeneSetSize = minGeneSetSize,
                       maxGeneSetSize = maxGeneSetSize,
                       print_genes = print_genes)
    if(print_genes == TRUE){
      prof2<- prof %>%
        dplyr::select(pathway, pvalue, `k/K`, genes) %>%
        dplyr::rename_with(~paste0("pvalue_", i),
                           "pvalue",
                           recycle0 = TRUE) %>%
        dplyr::rename_with(~paste0("k/K_", i),
                           "k/K",
                           recycle0 = TRUE) %>%
        dplyr::rename_with(~paste0("genes_", i),
                           "genes",
                           recycle0 = TRUE)

    } else{
      prof2<- prof %>%
        dplyr::select(pathway, pvalue, `k/K`) %>%
        dplyr::rename_with(~paste0("pvalue_", i),
                           "pvalue",
                           recycle0 = TRUE) %>%
        dplyr::rename_with(~paste0("k/K_", i),
                           "k/K",
                           recycle0 = TRUE)
    }

    iter_list <- prof2

  }
  parallel::stopCluster(cl)

  for(i in 1:length(iter_list)){
    df <- as.data.frame(iter_list[[i]])
    base_df <- base_df %>%
      dplyr::full_join(df, by = c("pathway" = "pathway"))
  }

  df_p <- base_df %>%
    dplyr::select(pathway, dplyr::starts_with("pvalue"))

  df_k <- base_df %>%
    dplyr::select(pathway, dplyr::starts_with("k/K"))

  ### get median p-values, perform correction, format result ###

  # remove rows with all NA values (no enrichment in any iteration)
  df_p <- df_p[rowSums(is.na(df_p[2:(niter+1)])) != (ncol(df_p)-1), ]
  df_k <- df_k[rowSums(is.na(df_k[2:(niter+1)])) != (ncol(df_k)-1), ]

  if(print_genes == TRUE){
    df_genes <- base_df %>%
      dplyr::select(pathway, dplyr::starts_with("genes")) %>%
      dplyr::filter(pathway %in% df_p$pathway)

    results_vec_genes <- rep(NA, nrow(df_genes))
    for(i in 1:nrow(df_genes)){
      genes_in_row <- c()
      for(j in 2:ncol(df_genes)){
        if(!is.na(df_genes[i,j])){
          el <- stringr::str_split_1(df_genes[i,j], ";")
        } else{
          el <- NA
        }
        genes_in_row <- unique(c(genes_in_row, el))
      }
      results_vec_genes[i] <- paste0(genes_in_row, collapse = ";")
    }
  }


  results_df_k <- replace(df_k, is.na(df_k),  0)# NA overlaps recoded as 0
  results_df_k <- results_df_k %>%
    dplyr::rowwise() %>%
    dplyr::mutate(median_kK = stats::median(dplyr::c_across(dplyr::where(~is.numeric(.x))), na.rm=FALSE)) %>% # calculate medians across iterations
    dplyr::ungroup() %>%
    dplyr::select(pathway, median_kK)

  results_df_summary <- replace(df_p, is.na(df_p), 1) # NA p-values recoded as 1
  results_df_summary <- results_df_summary %>%
    dplyr::rowwise() %>%
    dplyr::mutate(median_pval = stats::median(dplyr::c_across(dplyr::where(~is.numeric(.x))), na.rm=FALSE)) %>% # calculate medians across iterations
    dplyr::mutate(max_pval = max(dplyr::c_across(dplyr::where(~is.numeric(.x))), na.rm = FALSE),
                  min_pval = min(dplyr::c_across(dplyr::where(~is.numeric(.x))), na.rm = FALSE)) %>%
    dplyr::ungroup() %>% #undo rowwise
    dplyr::select(pathway, median_pval, min_pval, max_pval) %>%
    dplyr::left_join(results_df_k, by = c("pathway" = "pathway"))
  results_df2_p <- df_p %>%
    dplyr::rowwise() %>%
    dplyr::mutate(median_pval = stats::median(dplyr::c_across(dplyr::where(~is.numeric(.x))), na.rm=FALSE)) %>% # calculate medians across iterations
    dplyr::ungroup()

  results_df2_k <- df_k %>%
    dplyr::rowwise() %>%
    dplyr::mutate(median_kK = stats::median(dplyr::c_across(dplyr::where(~is.numeric(.x))), na.rm=FALSE)) %>% # calculate medians across iterations
    dplyr::ungroup()

  # p adjustment
  results_df_summary[,toupper(p.adjust.method)] <- stats::p.adjust(results_df_summary$median_pval, method = p.adjust.method)

  results_df_summary <- results_df_summary %>%
    dplyr::rename(`k/K` = median_kK) %>%
    dplyr::mutate("group" = deparse(substitute(anno_df)), .before = pathway)

  if(print_genes == TRUE){
    results_df_summary <- results_df_summary %>%
      dplyr::mutate(overlap_in_any_iteration = results_vec_genes)
  }



  results <- list("summary" = results_df_summary, "p_iterations" = results_df2_p, "k/K_iterations" = results_df2_k)

  return(results)
}
