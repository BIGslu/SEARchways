#' clusterProfiler of significant genes in Broad gene sets
#'
#' @param gene_list Named list object with gene symbols
#' @param gene_df Data frame including variable/module groups (column 1: group) and gene name (column2: gene). Can be used instead of gene_list
#' @param ID Character string for type of ID used in gene_list. One of SYMBOL, ENTREZ, ENSEMBL. Default is SYMBOL
#' @param species Character string denoting species of interest. Default is "human"
#' @param collection Character string denoting Broad gene set database
#' @param subcollection Character string denoting Broad gene set sub-database \cr
#' \tabular{rrrrr}{
#'  \strong{collection} \tab    \strong{subcollection}\cr
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
#' @return Data frame of enrichments including pathway, significance, and genes in pathway
#' @export
#'
#' @examples
#' #Get gene names for enrichment
#' gene_list <- list(HRV1 = names(example.gene.list[[1]]),
#'                   HRV2 = names(example.gene.list[[2]]))
#' BIGprofiler(gene_list, ID="ENSEMBL", collection="H", species = "human")
#'
#' # Use gene_df
#' gene_df <- data.frame(gs_name = c(rep("HRV1", 100), rep("HRV2",100)),
#'                       gene = c(names(example.gene.list[[1]]),
#'                                names(example.gene.list[[2]])))
#' BIGprofiler(gene_df=gene_df, ID="ENSEMBL", collection="H", species = "human")
#'
#' #Use custom data base
#' db <- data.frame(module = c(rep("module1",10), rep("module2",10)),
#'                  symbol = sample(gene_list[[1]], 20))
#'
#' BIGprofiler(gene_list, ID="ENSEMBL", db=db, species = "human")

BIGprofiler <- function(gene_list = NULL, gene_df = NULL, ID = "SYMBOL",
                        species = "human",
                        collection = NULL, subcollection = NULL,
                        db = NULL){
  db_join <- pathway_GOID <- gs_exact_source <- BgRatio <- Description <- FDR <- GeneRatio <- ensembl_gene <- entrez_gene <- geneID <- gene_symbol <- genes <- group <- group_in_cat.subcat <- group_in_pathway <- gs_collection <- gs_subcollection <- gs_name <- `k/K` <- p.adjust <- pathway <- pval <- pvalue <- qvalue <- size_cat.subcat <- size_group <- size_pathway <- db_species <- NULL

  #Recode species
  if(species %in% c("human", "Homo sapiens", "HS")){
    db_species <- "HS"
    species <- "human"
  }
  if(species %in% c("mouse", "Mus musculus", "MM")){
    db_species <- "MM"
    species = "mouse"
  }

  ##### Database #####
  #Load gene ontology
  if(!is.null(collection)){
    db.format <- msigdbr::msigdbr(species = species, collection = collection, db_species = db_species)
    #Subset subcollection if selected
    if(!is.null(subcollection)){
      db.format <- db.format %>%
        dplyr::filter(grepl(paste0("^",subcollection), gs_subcollection))
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
    stop("Please provide gene set information as Broad collection/subcollection or in a data frame as db.")
  }

  #Get gene ID
  if(ID == "SYMBOL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, gene_symbol, gs_collection, gs_subcollection)
  } else if(ID == "ENSEMBL"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, ensembl_gene, gs_collection, gs_subcollection)
  } else if(ID == "ENTREZ"){
    db.format2 <- db.format %>%
      dplyr::select(gs_name, entrez_gene, gs_collection, gs_subcollection)
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
    #run enrichment on gene list
    enrich.result <- clusterProfiler::enricher(gene=gene_list_format[[g]],
                                               TERM2GENE=db.format2[,1:2])

    #handle no enrichment results
    if(is.null(enrich.result)){
      if(!is.null(db)){
        result.clean <- data.frame(
          group=g,
          gs_collection="custom",
          pathway="No enriched terms")
      } else{
        result.clean <- data.frame(
          group=g,
          gs_collection=collection,
          gs_subcollection=subcollection,
          pathway="No enriched terms")
      }

    } else{
      #Format collection labels
      db.species.clean <- db.format2 %>%
        dplyr::distinct(gs_collection, gs_subcollection, gs_name) %>%
        dplyr::rename(pathway=gs_name)

      #Format results
      #Format gene column to vector
      result.clean <- enrich.result@result %>%
        tibble::remove_rownames() %>%
        dplyr::rename(pathway=Description, FDR=p.adjust, pval=pvalue) %>%
        dplyr::mutate(genes = strsplit(geneID, split="/")) %>%
        #Extract values from ratios
        tidyr::separate(BgRatio, into=c("size_pathway","size_cat.subcat"), sep="/") %>%
        tidyr::separate(GeneRatio, into=c("group_in_pathway",
                                          "group_in_cat.subcat"),
                        sep="/") %>%
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
                      gs_collection, gs_subcollection, size_cat.subcat,
                      group_in_cat.subcat,
                      pathway, size_pathway, group_in_pathway, `k/K`,
                      pval, FDR, qvalue, genes) %>%
        dplyr::arrange(FDR)

      # add GO term reference ID to results
      if(!is.null(collection)){
        if(collection == "C5"){
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
    }}

  ##### Save results #####
  #combine list of df results
  results.all.df <- dplyr::bind_rows(all.results)
  return(results.all.df)
}
