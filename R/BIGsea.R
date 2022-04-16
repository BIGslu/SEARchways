#' @title BIGsea
#' @description Run gene set enrichment analysis (GSEA) on multiple list using fgsea. More information about GSEA and geneset databases check https://www.gsea-msigdb.org/gsea/msigdb
#' @param gene_list Named list object with named numeric vectors of gene symbols and logFC
#' @param ID Character string for type of ID used in gene_list. One of SYMBOL, ENTREZ, ENSEMBL. Default is SYMBOL
#' @param nperm Numeric permutations for P-value calculations. Default is 1000
#' @param species Character string denoting species of interest. Default is human
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
#' @param db If not using Broad databases, a data frame with gene ontology including category (gs_cat), subcategory (gs_subcat), gene set name (gs_name), and gene ID (gene_symbol, entrez_gene, or ensembl_gene as matches your gene_list names)
#'
#' @return Output of a list of list with GSEA results
#' @export
#'
#' @examples
#' BIGsea(example_gene_list, category="H", ID="ENSEMBL")
#' BIGsea(example_gene_list, category="C2", subcategory="CP", ID="ENSEMBL")

BIGsea <- function(gene_list, nperm=1000, species="human", ID="SYMBOL",
                   category = NULL, subcategory = NULL, db = NULL){
  ensembl_gene <- entrez_gene <- gene_symbol <- group <- gs_name <- gs_subcat <- padj <- pathway <- NULL
  #Blank list to hold results
  all.results <- list()

  #### Data ####
  #Load gene ontology
  if(!is.null(category)){
    db.format <- msigdbr::msigdbr(species, category)
    #Subset subcategory if selected
    if(!is.null(subcategory)){
      db.format <- db.format %>%
        dplyr::filter(grepl(paste0("^",subcategory), gs_subcat))
    }
  } else if(!is.null(db)){
    db.format <- db
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

  #### Loop ####
  #Loop through each list in the gene_list object
  for(genes in names(gene_list)){
    message(genes)
    #Extract 1 gene list
    genes.temp <- gene_list[[genes]]
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
      dplyr::mutate(gs_cat=category, gs_subcat=subcategory, .before=1)

    #### Save ####
    all.results[[genes]] <- fg.result
  }

  #### Format output ####
  #Unlist results into 1 df
  all.results.df <- do.call(rbind.data.frame, all.results) %>%
    tibble::rownames_to_column("group") %>%
    dplyr::mutate(group = gsub("[.][0-9]{0,4}","",group)) %>%
    dplyr::mutate(pathway = sub("[A-Z]*_","",pathway)) %>%
    dplyr::mutate(pathway = gsub("_"," ",pathway))

  return(all.results.df)
}
