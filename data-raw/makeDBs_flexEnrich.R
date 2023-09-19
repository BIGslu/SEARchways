require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(tidyverse)

# human
default_human_bg_db <- msigdbr::msigdbr("human")
symbol.human.db.full <- unique(default_human_bg_db$gene_symbol)
ensembl.human.db.full <- unique(default_human_bg_db$ensembl_gene)
entrez.human.db.full <- unique(as.character(default_human_bg_db$entrez_gene))

## protein coding
### symbol
human_symbol_pc_filter <- hgnc::import_hgnc_dataset(file = hgnc::latest_archive_url())
human_symbol_pc_filter <- human_symbol_pc_filter[which(human_symbol_pc_filter$symbol %in% symbol.human.db.full & human_symbol_pc_filter$locus_group == "protein-coding gene"),]
symbol.human.db.pc <- unique(symbol.human.db.full[which(symbol.human.db.full %in% human_symbol_pc_filter$symbol)])

### ensembl
ensembl_mart_human <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
human_ensembl_pc_filter <- biomaRt::getBM(mart = ensembl_mart_human, attributes = c("ensembl_gene_id", "transcript_biotype"))
human_ensembl_pc_filter <- human_ensembl_pc_filter[which(human_ensembl_pc_filter$ensembl_gene_id %in% ensembl.human.db.full & human_ensembl_pc_filter$transcript_biotype == "protein_coding"),]
ensembl.human.db.pc <- ensembl.human.db.full[which(ensembl.human.db.full %in% human_ensembl_pc_filter$ensembl_gene_id)]

### entrez
human_entrez_pc_filter <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez.human.db.full,
                                                columns = c("ENTREZID", "GENETYPE"),
                                                keytype = "ENTREZID")
human_entrez_pc_filter<- human_entrez_pc_filter[which(human_entrez_pc_filter$ENTREZID %in% entrez.human.db.full & human_entrez_pc_filter$GENETYPE == "protein-coding"),]
entrez.human.db.pc <- entrez.human.db.full[which(entrez.human.db.full %in% human_entrez_pc_filter$ENTREZID)]

# mouse
default_mouse_bg_db <- msigdbr::msigdbr("mouse")
symbol.mouse.db.full <- unique(default_mouse_bg_db$gene_symbol)
ensembl.mouse.db.full <- unique(default_mouse_bg_db$ensembl_gene)
entrez.mouse.db.full <- unique(as.character(default_mouse_bg_db$entrez_gene))

## protein coding
### ensembl
ensembl_mart_mouse <- biomaRt::useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

mouse_ensembl_pc_filter <- biomaRt::getBM(mart = ensembl_mart_mouse, attributes = c("ensembl_gene_id", "transcript_biotype"))

mouse_ensembl_pc_filter <- mouse_ensembl_pc_filter[which(mouse_ensembl_pc_filter$ensembl_gene_id %in% ensembl.mouse.db.full & mouse_ensembl_pc_filter$transcript_biotype == "protein_coding"),]
ensembl.mouse.db.pc <- unique(ensembl.mouse.db.full[which(ensembl.mouse.db.full %in% mouse_ensembl_pc_filter$ensembl_gene_id)])
### entrez
mouse_entrez_pc_filter <- AnnotationDbi::select(org.Mm.eg.db, keys = entrez.mouse.db.full,
                                                columns = c("ENTREZID", "GENETYPE"),
                                                keytype = "ENTREZID")
mouse_entrez_pc_filter<- mouse_entrez_pc_filter[which(as.character(mouse_entrez_pc_filter$ENTREZID) %in% as.character(entrez.mouse.db.full) & mouse_entrez_pc_filter$GENETYPE == "protein-coding"),]
entrez.mouse.db.pc <- unique(entrez.mouse.db.full[which(as.character(entrez.mouse.db.full) %in% as.character(mouse_entrez_pc_filter$ENTREZID))])

### symbol
mouse_symbol_pc_filter <- AnnotationDbi::select(org.Mm.eg.db, keys = symbol.mouse.db.full,
                                                columns = c("SYMBOL", "GENETYPE"),
                                                keytype = "SYMBOL")
mouse_symbol_pc_filter<- mouse_symbol_pc_filter[which(as.character(mouse_symbol_pc_filter$SYMBOL) %in% as.character(symbol.mouse.db.full) & mouse_symbol_pc_filter$GENETYPE == "protein-coding"),]
symbol.mouse.db.pc <- unique(symbol.mouse.db.full[which(as.character(symbol.mouse.db.full) %in% as.character(mouse_symbol_pc_filter$SYMBOL))])

# save
## human
### full
write.csv(symbol.human.db.full, file = "data-raw/symbol.human.db.full.csv", row.names = F)
write.csv(entrez.human.db.full, file = "data-raw/entrez.human.db.full.csv", row.names = F)
write.csv(ensembl.human.db.full, file = "data-raw/ensembl.human.db.full.csv", row.names = F)
### protein coding
write.csv(symbol.human.db.pc, file = "data-raw/symbol.human.db.pc.csv", row.names = F)
write.csv(entrez.human.db.pc, file = "data-raw/entrez.human.db.pc.csv", row.names = F)
write.csv(ensembl.human.db.pc, file = "data-raw/ensembl.human.db.pc.csv", row.names = F)

## mouse
### protein coding
write.csv(symbol.mouse.db.pc, file = "data-raw/symbol.mouse.db.pc.csv", row.names = F)
write.csv(entrez.mouse.db.pc, file = "data-raw/entrez.mouse.db.pc.csv", row.names = F)
write.csv(ensembl.mouse.db.pc, file = "data-raw/ensembl.mouse.db.pc.csv", row.names = F)


