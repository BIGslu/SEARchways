#' @title BIGsea
#' @description Run gene set enrichment analysis (GSEA) on multiple list using fgsea. More information about GSEA and geneset databases check https://www.gsea-msigdb.org/gsea/msigdb
#' @param gene_list Named list object with named numeric vectors of gene symbols and logFC. List names must match dat$targetes column names if using rand="label"
#' @param gene_df Data frame including variable/module groups (column 1: group), gene name (column2: gene), and log fold change (column 3: logFC). Can be used instead of gene_list. group levels must match dat$targetes column names if using rand="label"
#' @param dat Voom object with gene expression and sample metadata. Used in label randomization only
#' @param ID Character string for type of ID used in gene_list. One of SYMBOL, ENTREZ, ENSEMBL. Default is "SYMBOL"
#' @param rand Character string specifying the type of randomization to create the null distribution. One of "multi" (default, randomize samples in groups and genes in gene sets), "label" (randomize genes in gene sets. Recommended for small sample groups), "simple" (randomize samples in group)
#' @param rand_var Character string specifying the variable to randomize in the label method
#' @param rand_est Data frame with random variable estimates. From prior run of BIGsea in the estimates slot
#' @param nperm Numeric permutations for P-value calculations. Default is 1000
#' @param species Character string denoting species of interest. Default is "human"
#' @param category Character string denoting Broad gene set database
#' @param subcategory Character string denoting Broad gene set sub-database. See https://www.gsea-msigdb.org/gsea/msigdb/
#' @param pw Character vector of pathway names to include. Must still provide category/subcategory. Format must exact match database such as HALLMARK_INTERFERON_GAMMA_RESPONSE
#' @param db If not using Broad databases, a data frame with gene ontology including gene set name (column 1: gs_name) and gene ID (column2: gene_symbol, entrez_gene, or ensembl_gene as matches your gene_list names)
#' @param minGeneSetSize Maximum overlap between a gene set and your list of query genes for hypergeometric enrichment to be calculated. Default is 10.
#' @param maxGeneSetSize Maximum size of a reference gene set for hypergeometric enrichment to be calculated. Default is 1e10
#' @param processors Numeric total processors to use. Default is 1
#' @param ... Additional parameters for kimma::kmFit needed if rand = "label"
#'
#' @return Data frame of enrichments including pathway, significance, and leading edge genes
#' @export
#'
#' @examples
#' BIGsea(example.gene.list, category="H", ID="ENSEMBL") #no result
#' BIGsea(example.gene.list, category="C2", subcategory="CP", ID="ENSEMBL")
#'
#' #Use gene_df. No overlap message
#' gene_df <- data.frame(gs_name = c("HRV1", rep("HRV2",100)),
#'                       gene = c("notReal",
#'                                names(example.gene.list[[2]])),
#'                      logFC = c(example.gene.list[[1]][1],
#'                                example.gene.list[[2]]))
#' BIGsea(gene_df=gene_df, category="C2", subcategory="CP", ID="ENSEMBL",
#'        rand="simple")
#'
#' #Use custom data base
#' db <- data.frame(module = c(rep("module1",10), rep("module2",10)),
#'                  symbol = sample(names(example.gene.list[[1]]), 20))
#' BIGsea(example.gene.list, ID="ENSEMBL", db=db)
#'
#' #Use label randomization
#' gene_df <- data.frame(gs_name = rep("virus",100),
#'                       gene = c(names(example.gene.list[[2]])),
#'                      logFC = c(example.gene.list[[2]]))
#' example.voom <- kimma::example.voom
#' test <- BIGsea(gene_df = gene_df, dat=example.voom, ID="ENSEMBL",
#'        category="C2", subcategory="CP",
#'        rand="label", rand_var="virus",
#'        model="~virus+median_cv_coverage",
#'        run_lm=TRUE, use_weights=TRUE,
#'        nperm=2, pw=c("REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION"),
#'        processors=2)
#'
#' #Use pre-calculated random estimates
#' BIGsea(gene_df = gene_df, dat=example.voom, ID="ENSEMBL",
#'        category="C2", subcategory="CP",
#'        rand="label", rand_var="virus", rand_est=test[['estimates']],
#'        model="~virus+median_cv_coverage",
#'        run_lm=TRUE, use_weights=TRUE,
#'        nperm=2, pw=c("REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION"))

BIGsea <- function(gene_list = NULL, gene_df = NULL,
                   dat = NULL,
                   rand="multi", nperm=1000,
                   rand_var=NULL, rand_est=NULL,
                   species="human", ID="SYMBOL",
                   category = NULL, subcategory = NULL, pw = NULL, db = NULL,
                   minGeneSetSize = 10, maxGeneSetSize = 1e10,
                   processors = 1, ...){
  gs_exact_source <- db_join <- pathway_GOID <- ensembl_gene <- entrez_gene <- gene_symbol <- gs_name <- gs_collection <- gs_subcollection <- padj <- pathway <- col1 <- db_species <- NULL

  #Blank list to hold results
  all.results <- list()

  #### Database ####
  #Load gene ontology
  if(!is.null(category)){
    #Check that category exists in msigdb
    all_cat <- msigdbr::msigdbr_collections() %>%
      dplyr::pull(gs_collection) %>% unique()
    if(!category %in% all_cat){
      stop("Category does not exist. Use msigdbr::msigdbr_collections() to see options.") }

    #Recode species
    if(species == "human"){
      species <- "Homo sapiens"
      db_species <- "HS"
    }
    if(species == "mouse"){
      species <- "Mus musculus"
      db_species <- "MM"
    }
    db.format <- msigdbr::msigdbr(species, db_species, collection=category)
    #Subset subcategory if selected
    if(!is.null(subcategory)){
      #Check that subcategory exists in msigdb
      all_subcat <- msigdbr::msigdbr_collections() %>%
        dplyr::pull(gs_subcollection) %>% unique()
      if(!subcategory %in% all_subcat){
        stop("Subcategory does not exist. Use msigdbr::msigdbr_collections() to see options.") }

      db.format <- db.format %>%
        dplyr::filter(grepl(paste0("^",subcategory), gs_subcollection))
    }
  } else if(!is.null(db)){
    db.format <- db %>%
      dplyr::mutate(gs_collection = "custom", gs_subcollection=NA)
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

  #Filter select pathways
  if(!is.null(pw)){
    db.format <- db.format %>% dplyr::filter(gs_name %in% pw)
    if(nrow(db.format)==0){
      stop("No pathways present in data base. Please check spelling in pw.")
    }
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
    gene_list_format <- NULL
    stop("Please provide either gene_list or gene_df.")
  }

  #### Loop ####
  #Loop through each list in the gene_list_format object
  for(g in names(gene_list_format)){
    message(paste("Running",g))
    if(rand %in% c("simple","multi")){
      #Extract 1 gene list
      genes.temp <- gene_list_format[[g]]
      #Check that genes exist in db
      gene_list_overlap <- base::intersect(names(genes.temp), unlist(db.ls))

      if(length(gene_list_overlap)>0){
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
        if(rand == "multi"){
          fg.result <- fgsea::fgseaMultilevel(
            pathways = db.ls, stats = genes.temp,
            nPermSimple = nperm,
            minSize = minGeneSetSize, maxSize = maxGeneSetSize,
            sampleSize = 101,
            #eps=0,
            scoreType = scoreType, nproc = processors,
            ...) %>%
            as.data.frame() %>%
            dplyr::mutate(method="multi", .before=0)
        } else if(rand == "simple"){
          fg.result <- fgsea::fgseaSimple(
            pathways = db.ls, stats = genes.temp,
            nperm = nperm,
            minSize = minGeneSetSize, maxSize = maxGeneSetSize,
            #eps=0,
            scoreType = scoreType, nproc = processors,
            ...) %>%
            as.data.frame() %>%
            dplyr::mutate(method="simple", .before=0)
        }
      } else {
        fg.result <- tibble::tibble(
          group=g, gs_collection=category, gs_subcollection=subcategory,
          pathway="No overlap of query genes and specified database.") }
      est.result <- NULL
    } else if(rand == "label"){
      #Check gene list/df name matches variable in targets data
      if(!g %in% colnames(dat$targets)){
        stop("All group names in gene_df or gene_list must be columns in dat$targets")
        }

      #Extract 1 gene list
      genes.temp <- gene_list_format[[g]]

      #### error about names
      #Check that genes exist in db
      gene_list_overlap <- base::intersect(rownames(dat$E), unlist(db.ls))

      if(length(gene_list_overlap)>0){
        #### FGSEA ####
        #Run GSEA with fgsea
        label.result <- SEARchways::fgseaLabel2(
          estimates = genes.temp,
          pathways = db.ls,
          dat = dat, label = g,
          rand_var = rand_var, rand_est = rand_est,
          nperm = nperm,
          minSize = minGeneSetSize, maxSize = maxGeneSetSize,
          #eps=0,
          nproc = processors, ...)

        fg.result <- label.result[["pvals"]] %>%
          as.data.frame() %>%
          dplyr::mutate(method="label", .before=0)
        est.result <- label.result[["estimates"]]
      } else {
        all.results[[g]] <- tibble::tibble(
          group=g, gs_collection=category, gs_subcollection=subcategory,
          pathway="No overlap of query genes and specified database.") }
    } else{
      stop("rand must be set to one of multi, simple, or label") }

    fg.result <- fg.result %>%
      dplyr::rename(FDR=padj) %>%
      dplyr::mutate(group=g, gs_collection=category, gs_subcollection=subcategory, .before=1)

    # add GO term reference ID to results
    if(!is.null(category)){
      if(category == "C5"){
        db_join <- db.format %>%
          dplyr::select(c("gs_name", "gs_exact_source")) %>%
          dplyr::distinct()
        fg.result <- fg.result %>%
          dplyr::left_join(db_join, by = c("pathway" = "gs_name")) %>%
          dplyr::rename(pathway_GOID = gs_exact_source) %>%
          dplyr::relocate(pathway_GOID, .after = pathway)
      }}

    #### Save ####
    all.results[[g]] <- fg.result
  }

  #### Format output ####
  #Unlist results into 1 df
  all.results.df <- dplyr::bind_rows(all.results)

  if(rand %in% c("simple","multi")){
    return(all.results.df)
  } else if (rand=="label"){
    all.results.ls <- list()
    all.results.ls[["gsea"]] <- all.results.df
    all.results.ls[["estimates"]] <- est.result
    return(all.results.ls)
  }
}

