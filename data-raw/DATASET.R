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
