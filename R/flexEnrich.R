#' Flexible hypergeometric enrichment
#' @description
#' Run hypergeometric enrichment on multiple lists using phyper. Permits gene symbols, Ensembl IDs, or Entrez IDs as gene identifiers and allows a background of protein-coding mouse genes or either whole-genome or protein-coding genes for human data. Also allows for custom backgrounds. Performs hypergeometric enrichment against a Broad Gene Set database or a custom database.
#'
#' @param gene_list Named list object with gene IDs
#' @param gene_df Data frame including variable/module groups (column 1: group), and gene name (column2: gene). Can be used instead of gene_list
#' @param ID Character string for type of ID used in gene_list. One of SYMBOL, ENTREZ, ENSEMBL. Default is "SYMBOL"
#' @param species Character string denoting species of interest. "human" or "mouse" Default is "human"
#' @param category Character string denoting Broad gene set database
#' @param subcategory Character string denoting Broad gene set sub-database
#' @param db Custom database
#' @param custom_bg Custom background. Formatted as a vector of gene IDs.
#' @param protein_coding TRUE or FALSE: do you want to limit the background to only protein-coding genes? Default is TRUE
#' @param minOverlap Minimum overlap between a gene set and your list of query genes for hypergeometric enrichment to be calculated. Default is 1.
#' @param minGeneSetSize Maximum overlap between a gene set and your list of query genes for hypergeometric enrichment to be calculated. Default is 10.
#' @param maxGeneSetSize Maximum size of a reference gene set for hypergeometric enrichment to be calculated. Default is 1e10
#' @param print_genes TRUE or FALSE. Do you want the results to include a list of genes that overlap between any given gene set and your query genes. Default is TRUE. Leaving this parameter as TRUE will make the function run slowly for very large datasets.
#' @author Madison Cox
#'
#' @return Data frame enrichments including pathway, significance, and k/K ratios
#' @export
#'
#' @examples
#' gene_list <- list(HRV1 = names(example.gene.list[[1]]),
#'                   HRV2 = names(example.gene.list[[2]]))
#' flexEnrich(gene_list = gene_list, category = "H", ID = "ENSEMBL")

flexEnrich <- function(gene_list = NULL,
                       gene_df = NULL,
                       ID = "SYMBOL",
                       species = "human",
                       category = NULL,
                       subcategory = NULL,
                       db = NULL,
                       custom_bg = NULL,
                       protein_coding = TRUE,
                       minOverlap = 1,
                       minGeneSetSize = 10,
                       maxGeneSetSize = 1e10,
                       print_genes = TRUE){

  pathway_ID <- gs_exact_source <- FDR <- gs_name <- n <- db.format <- group <- n_query_genes <- n_background_genes <- gs_cat <- gs_subcat <- pathway <- n_pathway_genes <- n_query_genes_in_pathway <- `k/K` <- pvalue <- genes <- ensembl_gene <-  gene_symbol <- entrez_gene <- geneID <- NULL

  ##### Database #####
  #Load gene ontology

  if(!is.null(category)){
    db.format <- msigdbr::msigdbr(species, category)
    # remove gene sets that are too small or too large
    good_pw <- db.format %>%
      dplyr::group_by(gs_name) %>%
      dplyr::count() %>%
      dplyr::filter(n > minGeneSetSize,
                    n < maxGeneSetSize) %>%
      dplyr::pull(gs_name)
    db.format <- db.format %>%
      dplyr::filter(gs_name %in% good_pw)
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
      stop("Please use ID from SYMBOL, ENSEMBL, or ENTREZ.")
    }

    # remove gene sets that are too small or too large
    good_pw <- db.format %>%
      dplyr::group_by(gs_name) %>%
      dplyr::count() %>%
      dplyr::filter(n > minGeneSetSize,
                    n < maxGeneSetSize) %>%
      dplyr::pull(gs_name)
    db.format <- db.format %>%
      dplyr::filter(gs_name %in% good_pw)
  } else {
    stop("Please provide gene set information as Broad category/subcategory or in a data frame as db.")
  }


  ##### Get gene ID #####
  if(ID == "SYMBOL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, gene_symbol, gs_cat, gs_subcat) %>%
      dplyr::mutate(geneID = gene_symbol)
  } else if(ID == "ENSEMBL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, ensembl_gene, gs_cat, gs_subcat) %>%
      dplyr::mutate(geneID = ensembl_gene)
  } else if(ID == "ENTREZ"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, entrez_gene, gs_cat, gs_subcat) %>%
      dplyr::mutate(geneID = entrez_gene)
  } else{
    stop("Please use ID from SYMBOL, ENSEMBL, or ENTREZ.")
  }


  ##### Format data #####
  if(!is.null(gene_df)){
    gene_list_format <- list()
    col1 <- colnames(gene_df)[1]
    col2 <- colnames(gene_df)[2]
    for(g in unique(unlist(gene_df[,1]))){
      gene_list_format[[g]] <- gene_df %>%
        dplyr::filter(get(col1) == g) %>%
        dplyr::pull(get(col2)) %>% unique()
    }
  } else if(!is.null(gene_list)){
    gene_list_format <- gene_list

  } else{
    stop("Please provide either gene_list or gene_df.")
  }


  ##### Get background gene list #####
  if(!is.null(custom_bg)){
    if(length(base::intersect(unlist(gene_list_format), custom_bg)) == 0){
      stop("None of your input gene IDs are present in your background. Are you sure you used the same ID format?")
    }
    else{bg <- unique(custom_bg)}
  } else{
    if(!species %in% c("mouse","human")){
      stop("Please enter either 'human' or 'mouse' for species.")
    } else {
      if(ID == "SYMBOL"){
        if(protein_coding == TRUE){
          if(species == "human"){
            bg <- symbol.human.db.pc
          }
          else if(species == "mouse"){
            bg <- symbol.mouse.db.pc
          }
        }
        else{
          if(species == "human"){
            bg <- symbol.human.db.full
          }
          else if(species == "mouse"){
            stop("At this time, only protein-coding backgrounds are available for mouse genes. Please use 'protein_coding = TRUE'.")
          }
        }
      }
      else if(ID == "ENSEMBL"){
        if(protein_coding == TRUE){
          if(species == "human"){
            bg <- ensembl.human.db.pc
          }
          else if(species == "mouse"){
            bg <- ensembl.mouse.db.pc
          }
        }
        else{
          if(species == "human"){
            bg <- ensembl.human.db.full
          }
          else if(species == "mouse"){
            stop("At this time, only protein-coding backgrounds are available for mouse genes. Please use 'protein_coding = TRUE'.")
          }
        }
      }
      else if(ID == "ENTREZ"){
        if(protein_coding == TRUE){
          if(species == "human"){
            bg <- entrez.human.db.pc
          }
          else if(species == "mouse"){
            bg <- entrez.mouse.db.pc
          }
        }
        else{
          if(species == "human"){
            bg <- entrez.human.db.full
          }
          else if(species == "mouse"){
            stop("At this time, only protein-coding backgrounds are available for mouse genes. Please use 'protein_coding = TRUE'.")
          }
        }
      }
      else{stop("Please enter a valid ID. Options are SYMBOL, ENSEMBL, or ENTREZ")}
    }
  }

  ##### Get background value for stats::phyper #####
  n_background_genes <- length(unique(bg))


  ##### Loop through groups in gene_list_format #####
  all.results <- list()

  for(g in names(gene_list_format)){
    n_query_genes <- length(unique(gene_list_format[[g]][which(gene_list_format[[g]] %in% bg)]))
    query <- gene_list_format[[g]]
    n_genesets <- length(unique(db.format2$gs_name))

    # Blank holders
    pvals <- c()
    set_sizes <- c()
    overlaps <- c()
    set_names <- c()
    kK_ratios <- c()
    genes_in_overlap <- list()

    # Loop through gene sets
    for(s in unique(db.format2$gs_name)){
      set <- db.format2$geneID[which(db.format2$gs_name == s)]
      set_size = length(unique(set))
      query_in_set <- length(unique(base::intersect(set, query)))
      kK_ratio <- query_in_set/set_size
      if(print_genes){
        query_genes_in_set <- unique(base::intersect(set, query))
        query_genes_in_set <- list(sort(unique(base::intersect(set, query))))
      }

      ## print progress ##
      i <- which(unique(db.format2$gs_name) == s)
      if(i/1000 == round(i/1000)){
        print(paste0(i, " out of ", length(unique(db.format2$gs_name)), "gene sets complete"))
      }

      set_sizes <- c(set_sizes, set_size)
      overlaps <- c(overlaps, query_in_set)
      set_names <- c(set_names, s)
      kK_ratios <- c(kK_ratios, kK_ratio)

      if(print_genes) {
        if(!is.null(query_genes_in_set)){
          genes_in_overlap[[i]] <- query_genes_in_set

        } else{
          genes_in_overlap[[i]] <- NA
        }
      }

      if(maxGeneSetSize < set_size | set_size < minGeneSetSize | query_in_set < minOverlap){
        pvals <- c(pvals, NA)
      } else{
        p <- stats::phyper(query_in_set - 1, set_size, (n_background_genes - set_size), n_query_genes, lower.tail = F)
        pvals <- c(pvals, p)
      }
    }

    nrep <- length(unique(db.format2$gs_name))

    res.temp <- tibble::tibble("group" = rep(g, nrep),
                               "n_query_genes" = rep(n_query_genes, nrep),
                               "n_background_genes" = rep(n_background_genes, nrep),
                               "gs_cat" = rep(category, nrep),
                               "gs_subcat" = rep(subcategory, nrep),
                               "pathway" = set_names,
                               "n_pathway_genes" = set_sizes,
                               "n_query_genes_in_pathway" = overlaps,
                               "k/K" = kK_ratios,
                               "pvalue" = pvals #,
                               #     "genes" = genes_in_overlap
    )
    if(print_genes){
      res.temp$genes <- genes_in_overlap
    }


    res.temp <- res.temp %>%
      dplyr::filter(n_pathway_genes >= minGeneSetSize,
                    n_pathway_genes <= maxGeneSetSize,
                    n_query_genes_in_pathway >= minOverlap,
                    !is.na(pvalue))


    res.temp$FDR <- stats::p.adjust(res.temp$pvalue, method = "fdr")
    res.temp <- res.temp %>%
      dplyr::relocate(FDR, .after = pvalue)

    if(!is.null(category)){
      if(category == "C5"){
        res.temp <- res.temp %>%
          dplyr::left_join(dplyr::select(db.format, c("gs_name", "gs_exact_source")), by = c("pathway" = "gs_name")) %>%
          dplyr::rename(pathway_ID = gs_exact_source) %>%
          dplyr::relocate(pathway_ID, .after = pathway)
      }
    }

    all.results[[g]] <- res.temp
  }


  ##### Save results #####
  #combine list of df results
  results.all.df <- dplyr::bind_rows(all.results)
  return(results.all.df)
}
