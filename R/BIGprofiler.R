#' clusterProfiler of significant genes in Broad gene sets
#'
#' @param gene_list Named list object with gene symbols
#' @param gene_df Data frame including variable/module groups (column 1: group) and gene name (column2: gene). Can be used instead of gene_list
#' @param ID Character string for type of ID used in gene_list. One of SYMBOL, ENTREZ, ENSEMBL. Default is SYMBOL
#' @param species Character string denoting species of interest. Default is human
#' @param category Character string denoting Broad gene set database. See https://www.gsea-msigdb.org/gsea/msigdb/
#' @param subcategory Character string denoting Broad gene set sub-database. See https://www.gsea-msigdb.org/gsea/msigdb/
#' @param db If not using Broad databases, a data frame with gene ontology
#' including gene set name (column 1: gs_name) and gene ID (column2:
#' gene_symbol, entrez_gene, or ensembl_gene as matches your gene_list names)
#' @param minGSSize Numeric. Minimum size of genes annotated for testing (minimum pathway
#'   size). Default \code{10}. See [clusterProfiler::enricher()]
#' @param maxGSSize Numeric. Maximum size of genes annotated for testing (maximum pathway
#'   size). Default \code{500}. See [clusterProfiler::enricher()]
#'
#' @return Data frame of enrichments including pathway, significance, and genes
#'   in pathway
#' @export
#'
#' @examples
#' #Get gene names for enrichment
#' gene_list <- list(HRV1 = "notReal",
#'                   HRV2 = names(example.gene.list[[2]]))
#' BIGprofiler(gene_list, ID="ENSEMBL", category="H")
#'
#' # Use gene_df
#' gene_df <- data.frame(gs_name = c(rep("HRV1", 100), rep("HRV2",100)),
#'                       gene = c(names(example.gene.list[[1]]),
#'                                names(example.gene.list[[2]])))
#' BIGprofiler(gene_df=gene_df, ID="ENSEMBL", category="H")
#'
#' #Use custom data base
#' db <- data.frame(module = c(rep("module1",10), rep("module2",10)),
#'                  symbol = sample(gene_list[[2]], 20))
#'
#' BIGprofiler(gene_list, ID="ENSEMBL", db=db)

BIGprofiler <- function(gene_list = NULL, gene_df = NULL, ID = "SYMBOL",
                        species = "human",
                        category = NULL, subcategory = NULL,
                        db = NULL,
                        minGSSize = 10,
                        maxGSSize = 500){
  db_join <- pathway_GOID <- gs_exact_source <- BgRatio <- Description <- FDR <- GeneRatio <- ensembl_gene <- entrez_gene <- geneID <- gene_symbol <- genes <- group <- group_in_cat.subcat <- group_in_pathway <- gs_cat <- gs_name <- gs_subcat <- `k/K` <- p.adjust <- pathway <- pval <- pvalue <- qvalue <- size_cat.subcat <- size_group <- size_pathway <- gene_list_overlap <- NULL

  ##### Database #####
  #Load gene ontology
  if(!is.null(category)){
    #Check that category exists in msigdb
    all_cat <- msigdbr::msigdbr_collections() %>%
      dplyr::pull(gs_cat) %>% unique()
    if(!category %in% all_cat){
      stop("Category does not exist. Use msigdbr::msigdbr_collections() to see options.") }

    db.format <- msigdbr::msigdbr(species, category)
    #Subset subcategory if selected
    if(!is.null(subcategory)){
      #Check that subcategory exists in msigdb
      all_subcat <- msigdbr::msigdbr_collections() %>%
        dplyr::pull(gs_subcat) %>% unique()
      if(!subcategory %in% all_subcat){
        stop("Subcategory does not exist. Use msigdbr::msigdbr_collections() to see options.") }

      db.format <- db.format %>%
        dplyr::filter(grepl(paste0("^",subcategory), gs_subcat))
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

  #Get gene ID
  if(ID == "SYMBOL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, gene_symbol, gs_cat, gs_subcat)
  } else if(ID == "ENSEMBL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, ensembl_gene, gs_cat, gs_subcat)
  } else if(ID == "ENTREZ"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, entrez_gene, gs_cat, gs_subcat)
  } else{
    stop("Please like ID from SYMBOL, ENSEMBL, or ENTREZ.")
  }

  ##### Format data ####
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

  ##### Loop through gene df #####
  #Blank holders
  all.results <- list()

  for(g in names(gene_list_format)){
    print(g)

    #Check that genes exist in db
    gene_list_overlap <- base::intersect(gene_list_format[[g]],
                                         unlist(db.format2[,2]))
    if(length(gene_list_overlap)>0){
      #run enrichment on gene list
      enrich.result <- clusterProfiler::enricher(gene=gene_list_overlap,
                                                 TERM2GENE=db.format2[,1:2],
                                                 minGSSize = minGSSize,
                                                 maxGSSize = maxGSSize)

      #handle no enrichment results
      if(is.null(enrich.result)){
        if(!is.null(db)){
          result.clean <- data.frame(
            group=g,
            gs_cat="custom",
            pathway="No enriched terms")
        } else{
          result.clean <- data.frame(
            group=g,
            gs_cat=category,
            gs_subcat=subcategory,
            pathway="No enriched terms")
        }

      } else{
        #Format category labels
        db.species.clean <- db.format2 %>%
          dplyr::distinct(gs_cat, gs_subcat, gs_name) %>%
          dplyr::rename(pathway=gs_name)

        #Format results
        #Format gene column to vector
        result.clean <- enrich.result@result %>%
          tibble::remove_rownames() %>%
          dplyr::rename(pathway=Description, FDR=p.adjust, pval=pvalue) %>%
          dplyr::mutate(genes = strsplit(geneID, split="/")) %>%
          #Extract values from ratios
          tidyr::separate_wider_delim(BgRatio, names=c("size_pathway","size_cat.subcat"), delim="/") %>%
          tidyr::separate_wider_delim(GeneRatio, names=c("group_in_pathway",
                                            "group_in_cat.subcat"),
                          delim="/") %>%
          dplyr::mutate_at(dplyr::vars("size_pathway","size_cat.subcat",
                                       "group_in_pathway","group_in_cat.subcat"),
                           as.numeric) %>%
          #Calculate k/K
          dplyr::mutate("k/K"=group_in_pathway/size_pathway) %>%

          #Add ID columns for database names
          dplyr::left_join(db.species.clean, by = "pathway") %>%
          #Add columns for group info
          dplyr::mutate(group=g, size_group = length(gene_list_format[[g]])) %>%
          #Reorder variables
          dplyr::select(group, size_group,
                        gs_cat, gs_subcat, size_cat.subcat,
                        group_in_cat.subcat,
                        pathway, size_pathway, group_in_pathway, `k/K`,
                        pval, FDR, qvalue, genes) %>%
          dplyr::arrange(FDR)

        # add GO term reference ID to results
        if(!is.null(category)){
          if(category == "C5"){
            db_join <- db.format %>%
              dplyr::select(c("gs_name", "gs_exact_source")) %>%
              dplyr::distinct()
            result.clean <- result.clean %>%
              dplyr::left_join(db_join, by = c("pathway" = "gs_name")) %>%
              dplyr::rename(pathway_GOID = gs_exact_source) %>%
              dplyr::relocate(pathway_GOID, .after = pathway)
          }
        }

        #Run enrich and save to results list
        all.results[[g]] <- result.clean
    }
    }else {
      all.results[[g]] <- tibble::tibble(
        group=g, size_group = length(gene_list_format[[g]]),
        gs_cat=category, gs_subcat=subcategory,
        size_cat.subcat=NA, group_in_cat.subcat=0,
        pathway="No overlap of query genes and specified database.")
    }
    }

  ##### Save results #####
  #combine list of df results
  results.all.df <- dplyr::bind_rows(all.results)
  return(results.all.df)
}
