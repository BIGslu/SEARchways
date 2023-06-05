#### Example gene fold change list ####
library(kimma)
library(tidyverse)

# Example expression data from kimma
example.voom <- example.voom

# Calculate mean fold change
FC <- as.data.frame(example.voom$E) %>%
  rownames_to_column() %>%
  pivot_longer(-rowname, names_to = "libID") %>%
  inner_join(example.voom$targets) %>%
  select(rowname, ptID, virus, value) %>%
  pivot_wider(names_from = virus) %>%
  mutate(FC = HRV-none) %>%
  group_by(rowname) %>%
  summarise(meanFC=mean(FC, na.rm = TRUE), .groups="drop")

# select 2 sets of 100 genes
# Format as list
set.seed(1879)
HRV1 <- FC %>%
  slice_sample(n=100) %>%
  pull(meanFC, name=rowname)

HRV2 <- FC %>%
  slice_sample(n=100) %>%
  pull(meanFC, name=rowname)

example.gene.list <- list(
  "HRV1" = HRV1,
  "HRV2" = HRV2
)

usethis::use_data(example.gene.list, overwrite = TRUE)

### Make flexEnrich backgrounds ###
# human
## full
symbol.human.db.full <- as.matrix(read.csv("data-raw/symbol.human.db.full.csv", header = T))[,1]
entrez.human.db.full <- as.matrix(read.csv("data-raw/entrez.human.db.full.csv", header = T))[,1]
ensembl.human.db.full <- as.matrix(read.csv("data-raw/ensembl.human.db.full.csv", header = T))[,1]
usethis::use_data(symbol.human.db.full, overwrite = TRUE)
usethis::use_data(entrez.human.db.full, overwrite = TRUE)
usethis::use_data(ensembl.human.db.full, overwrite = TRUE)

## protein coding
symbol.human.db.pc <- as.matrix(read.csv("data-raw/symbol.human.db.pc.csv", header = T))[,1]
entrez.human.db.pc <- as.matrix(read.csv("data-raw/entrez.human.db.pc.csv", header = T))[,1]
ensembl.human.db.pc <- as.matrix(read.csv("data-raw/ensembl.human.db.pc.csv", header = T))[,1]
usethis::use_data(symbol.human.db.pc, overwrite = TRUE)
usethis::use_data(entrez.human.db.pc, overwrite = TRUE)
usethis::use_data(ensembl.human.db.pc, overwrite = TRUE)

# mouse
## protein coding
symbol.mouse.db.pc <- as.matrix(read.csv("data-raw/symbol.mouse.db.pc.csv", header = T))[,1]
entrez.mouse.db.pc <- as.matrix(read.csv("data-raw/entrez.mouse.db.pc.csv", header = T))[,1]
ensembl.mouse.db.pc <- as.matrix(read.csv("data-raw/ensembl.mouse.db.pc.csv", header = T))[,1]
usethis::use_data(symbol.mouse.db.pc, overwrite = TRUE)
usethis::use_data(entrez.mouse.db.pc, overwrite = TRUE)
usethis::use_data(ensembl.mouse.db.pc, overwrite = TRUE)


