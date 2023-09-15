#' @title BIGsea
#' @description Run gene set enrichment analysis (GSEA) on multiple list using fgsea. More information about GSEA and geneset databases check https://www.gsea-msigdb.org/gsea/msigdb
#' @param gene_list Named list object with named numeric vectors of gene symbols and logFC
#' @param gene_df Data frame including variable/module groups (column 1: group), gene name (column2: gene), and log fold change (column 3: logFC). Can be used instead of gene_list
#' @param ID Character string for type of ID used in gene_list. One of SYMBOL, ENTREZ, ENSEMBL. Default is "SYMBOL"
#' @param nperm Numeric permutations for P-value calculations. Default is 1000
#' @param species Character string denoting species of interest. Default is "human"
#' @param category Character string denoting Broad gene set database
#' @param subcategory Character string denoting Broad gene set sub-database \cr
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
#' @param db If not using Broad databases, a data frame with gene ontology including gene set name (column 1: gs_name) and gene ID (column2: gene_symbol, entrez_gene, or ensembl_gene as matches your gene_list names)
#'
#' @return Data frame of enrichments including pathway, significance, and leading edge genes
#' @export
#'
#' @examples
#' BIGsea(example.gene.list, category="H", ID="ENSEMBL")
#' BIGsea(example.gene.list, category="C2", subcategory="CP", ID="ENSEMBL")
#'
#' #Use gene_df
#' gene_df <- data.frame(gs_name = c(rep("HRV1", 100), rep("HRV2",100)),
#'                       gene = c(names(example.gene.list[[1]]),
#'                                names(example.gene.list[[2]])),
#'                      logFC = c(example.gene.list[[1]],
#'                                example.gene.list[[2]]))
#' BIGsea(gene_df=gene_df, category="H", ID="ENSEMBL")
#'
#' #Use custom data base
#' db <- data.frame(module = c(rep("module1",10), rep("module2",10)),
#'                  symbol = sample(names(example.gene.list[[1]]), 20))
#' BIGsea(example.gene.list, ID="ENSEMBL", db=db)
#'

BIGsea <- function(gene_list = NULL, gene_df = NULL,
                   nperm=1000, species="human", ID="SYMBOL",
                   category = NULL, subcategory = NULL, db = NULL){
  pathway_ID <- gs_exact_source <- ensembl_gene <- entrez_gene <- gene_symbol <- group <- gs_name <- gs_subcat <- padj <- pathway <- col1 <- NULL
  #Blank list to hold results
  all.results <- list()

  #### Database ####
  #Load gene ontology
  if(!is.null(category)){
    db.format <- msigdbr::msigdbr(species, category)
    #Subset subcategory if selected
    if(!is.null(subcategory)){
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
      dplyr::group_by(gs_name) %>%
      dplyr::summarise(genes = list(gene_symbol)) %>%
      dplyr::ungroup() %>%  tibble::column_to_rownames("gs_name")
  } else if(ID == "ENSEMBL"){
    db.format2 <- db.format %>%
      dplyr::group_by(gs_name) %>%
      dplyr::summarise(genes = list(ensembl_gene)) %>%
      dplyr::ungroup() %>%  tibble::column_to_rownames("gs_name")
  } else if(ID == "ENTREZ"){
    db.format2 <- db.format %>%
      dplyr::group_by(gs_name) %>%
      dplyr::summarise(genes = list(entrez_gene)) %>%
      dplyr::ungroup() %>%  tibble::column_to_rownames("gs_name")
  } else{
    stop("Please like ID from SYMBOL, ENSEMBL, or ENTREZ.")
  }

  #Format database to list
  db.ls <- lapply(as.list(1:dim(db.format2)[1]),
                  function(x) unlist(db.format2[x[1],]))
  names(db.ls) <- rownames(db.format2)

  #### Format data ####
  if(!is.null(gene_df)){
    gene_list_format <- list()
    col1 <- colnames(gene_df)[1]
    for(g in unique(unlist(gene_df[,1]))){
      temp <- gene_df %>%
        dplyr::filter(get(col1) == g) %>%
        dplyr::distinct()
    gene_vec <- unlist(temp[,3])
    names(gene_vec) <- unlist(temp[,2])

    gene_list_format[[g]] <- gene_vec
    }
  } else if(!is.null(gene_list)){
    gene_list_format <- gene_list
  } else{
    stop("Please provide either gene_list or gene_df.")
  }


  #### Loop ####
  #Loop through each list in the gene_list_format object
  for(g in names(gene_list_format)){
    message(g)
    #Extract 1 gene list
    genes.temp <- gene_list_format[[g]]
    #Order by fold change
    genes.temp <- sort(genes.temp, decreasing = TRUE)

    #### FGSEA ####
    #Set score type based on fold change
    if(min(genes.temp) < 0 & max(genes.temp) > 0){
      scoreType <- "std"
    } else if(max(genes.temp) <= 0){
      scoreType <- "neg"
    } else if(min(genes.temp) >= 0){
      scoreType <- "pos"
    } else{
      stop("Could not determine score type from fold changes.")
    }

    #Run GSEA with fgsea
    fg.result <- fgsea::fgseaSimple(pathways = db.ls,
                                    stats = genes.temp,
                                    nperm=nperm,
                                    #eps=0,
                                    scoreType=scoreType) %>%
      as.data.frame() %>%
      dplyr::rename(FDR=padj) %>%
      dplyr::mutate(group=g, gs_cat=category, gs_subcat=subcategory, .before=1)

    # add GO term reference ID to results
    if(!is.null(category)){
      if(category == "C5"){
        fg.result <- fg.result %>%
          dplyr::left_join(dplyr::select(db.format, c("gs_name", "gs_exact_source")), by = c("pathway" = "gs_name")) %>%
          dplyr::rename(pathway_ID = gs_exact_source) %>%
          dplyr::relocate(pathway_ID, .after = pathway)
      }
    }


    #### Save ####
    all.results[[g]] <- fg.result
  }

  #### Format output ####
  #Unlist results into 1 df
  all.results.df <- dplyr::bind_rows(all.results)

  return(all.results.df)
}

