

# THIS PART IS FROM THE ENRICHR R FUNCTION https://github.com/wjawaid/enrichR

# ##' onLoad hook to setup package options
# ##'
# ##' onLoad hook to setup package options and to check connection to website
# ##' @title onLoad hook to setup package options
# ##' @param libname (Required). Library name
# ##' @param pkgname (Required). Package name
# ##' @return NULL
# ##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
# .onAttach <- function(libname, pkgname) {
#   options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
#   options(enrichR.live = TRUE)
#   packageStartupMessage("Welcome to enrichR\nChecking connection ... ", appendLF = TRUE)
#   options(modEnrichR.use = TRUE)
#   options(enrichR.sites.base.address = "https://maayanlab.cloud/")
#   options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"))
#   if (getOption("modEnrichR.use")) {
#     listEnrichrSites()
#   } else {
#     getEnrichr(url=paste0(getOption("enrichR.base.address"), "datasetStatistics"))
#     packageStartupMessage("Enrichr ... ", appendLF = FALSE)
#     if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")
#   }
# }


##' Helper function
##'
##' Helper function for GET
##' @title Helper function for GET
##' @param url (Required). URL address requested
##' @param ... (Optional). Additional parameters to pass to GET
##' @return same as GET
##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
##' @author I-Hsuan Lin \email{i-hsuan.lin@manchester.ac.uk}
##' @importFrom httr GET
##' @importFrom httr status_code
##' @importFrom httr http_status
getEnrichr <- function(url, ...) {
  options(enrichR.live = FALSE)
  tryCatch({
    x <- GET(url = url, ...)
    code <- status_code(x)
    if(code != 200) {
      # Error with status code
      message(http_status(code)$message)
    } else {
      # OK/success
      options(enrichR.live = TRUE)
      invisible(x)
    }
  },
  # Warning message
  warning = function(warn) {
    message(warn); message("") # force newline
  },
  # Error without status code
  error = function(err) {
    message(err); message("") # force newline
  },
  finally = function() {
    invisible(x)
  })
}

##' List modEnrichr Websites
##'
##' List Enrichr Websites
##' @title List Enrichr Websites
##' @return print Enrichr Website status
##' @author Alexander Blume
##' @param ... (Optional  Additional parameters)
listEnrichrSites <- function(...) {
  for (site in getOption("enrichR.sites")) {
    getEnrichr(url = paste0(getOption("enrichR.sites.base.address"), site, "/", "datasetStatistics"))
    packageStartupMessage(paste0(site, " ... "), appendLF = FALSE)
    if (paste0(getOption("enrichR.sites.base.address"), site, "/")  == getOption("enrichR.base.address")) {
      if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")
    } else
      if (getOption("enrichR.live")) packageStartupMessage("Connection is available!")

  }
}




# FROM HERE THE MODIFICATIONS WERE DONE


#' Run EnrichR function on multiple Geneset databases on either gene list or modules of intertest
#'
#' @param gene_list Here input either Genelist in SYMBOL or Ensembl ID format or input module object output from module function
#' @param gene_df Data frame including variable/module groups (column 1: group) and gene name (column2: gene). Can be used instead of gene_list
#' @param ID gene id is either SYMBOL, ENSEMBL, ENTREZ
#' @param dbs Enter the geneset database you want to run Enrichr on the default is MSigDB Hallmark 2020 to check available database run enrichR::listEnrichrDb()
#' @author Basilin Benson

#' @return list with enrichr result for each geneset database
#'
#' @export
#'
#' @examples
#' #Get gene names for enrichment
#' gene_list <- list(HRV1 = names(example.gene.list[[1]]),
#'                   HRV2 = names(example.gene.list[[2]]))
#' BIGenrichr(gene_list, ID = "ENSEMBL")
#'
#' # Use gene_df
#' gene_df <- data.frame(gs_name = c(rep("HRV1", 100), rep("HRV2",100)),
#'                       gene = c(names(example.gene.list[[1]]),
#'                       names(example.gene.list[[2]])))
#' BIGenrichr(gene_df=gene_df, ID="ENSEMBL")

BIGenrichr <- function(gene_list = NULL, gene_df = NULL,
                       ID = "SYMBOL",
                       dbs = c("MSigDB_Hallmark_2020")) {

  #### Intro messages ####
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  options(enrichR.quiet = TRUE)
  packageStartupMessage("Welcome to enrichR\nChecking connection ... ", appendLF = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"))
  if (getOption("modEnrichR.use")) {
    listEnrichrSites()
  } else {
    getEnrichr(url=paste0(getOption("enrichR.base.address"), "datasetStatistics"))
    packageStartupMessage("Enrichr ... ", appendLF = FALSE)
    if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")
  }




  #### Format data ####
  #Set variables to Null
  Adjusted.P.value <- Combined.Score <- FDR <- Genes <- Odds.Ratio <- Overlap <- P.value <- Term <- entrezgene_id <- genes <- group <- group_in_pathway <- `k/K` <- pathway <- pval <- size_group <- size_pathway <- ensembl_gene_id <- hgnc_symbol <- NULL
  #Setting the website for Enrichr to run
  enrichR::setEnrichrSite("Enrichr")

  #Convert to list if using gene_df
  if(!is.null(gene_df)){
    gene_list_format <- list()
    col1 <- colnames(gene_df)[1]
    col2 <- colnames(gene_df)[2]
    for(g in unique(gene_df[,1])){
      gene_list_format[[g]] <- gene_df %>%
        dplyr::filter(get(col1) == g) %>%
        dplyr::pull(get(col2)) %>% unique()
    }
  } else if(!is.null(gene_list)){
    gene_list_format <- gene_list

  } else{
    stop("Please provide either gene_list or gene_df.")
  }

  #Convert ENSEMBL and ENTREZ to HGNC SYMBOL
  if (ID == "SYMBOL") {
    gene_list_format2 <- gene_list_format
  }else if(ID == "ENSEMBL"){
    #Download Ensembl gene list to get HGNC symbols
    ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                                   dataset="hsapiens_gene_ensembl", mirror = "uswest")
    all_genes <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'),
                                mart = ensembl)
    #Convert within gene_list
    gene_list_format2 <- list()
    for(g in names(gene_list_format)){
      gene_list_format2[[g]] <- all_genes %>%
        dplyr::filter(ensembl_gene_id %in% gene_list_format[[g]]) %>%
        dplyr::pull(hgnc_symbol) %>% unique()
    }
  } else if(ID == "ENTREZ"){
    #Download Ensembl gene list to get HGNC symbols
    ensembl <- biomaRt::useEnsembl(biomart="ensembl",
                                   dataset="hsapiens_gene_ensembl")
    all_genes <- biomaRt::getBM(attributes=c('entrezgene_id', 'hgnc_symbol'),
                                mart = ensembl)
    #Convert within gene_list
    gene_list_format2 <- list()
    for(g in names(gene_list_format)){
      gene_list_format2[[g]] <- all_genes %>%
        dplyr::filter(entrezgene_id %in% gene_list_format[[g]]) %>%
        dplyr::pull(hgnc_symbol) %>% unique()
    }
  }else{
    stop("ID must be SYMBOL or ENSEMBL")
  }

  #Loop through the databases provided by user
  all.result.format.ls <- list()
  for (d in dbs){
    #Loop through the modules and genes within the modules
    all.result.format <- data.frame()
    for (g in names(gene_list_format2)) {
      #Running Enrichr
      result <- enrichR::enrichr(gene_list_format2[[g]], d)
      #Format results
      result.format <- result[[1]] %>%
        dplyr::rename(pathway=Term, pval=P.value, FDR=Adjusted.P.value) %>%
        tidyr::separate(Overlap, sep="/", into=c("group_in_pathway","size_pathway")) %>%
        dplyr::mutate_at(dplyr::vars("size_pathway","group_in_pathway"),
                         as.numeric) %>%
        #Calculate k/K
        dplyr::mutate("k/K"=group_in_pathway/size_pathway,
                      group=g,
                      size_group = length(gene_list_format2[[g]])) %>%
        #Split genes into vector
        dplyr::mutate(genes = strsplit(Genes, split=";")) %>%
        dplyr::select(group, size_group,
                      pathway, size_pathway, group_in_pathway, `k/K`,
                      pval, FDR, Odds.Ratio, Combined.Score, genes)
      all.result.format <- dplyr::bind_rows(result.format, all.result.format)
    }
    all.result.format.ls[[d]] <- all.result.format %>%
      dplyr::arrange(group, FDR)
  }
  return(all.result.format.ls)
}
