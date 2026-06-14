library(maaslin3)
library(tidyverse)
library(phyloseq)


sample_data(ps.relabund.pups) <- metadata %>% 
  mutate(PrePostFVT = recode(PrePostFVT, "Post-Final" = "Post")) %>%
  sample_data()

ps.relabund.pups.otutable <- ps.relabund.pups %>% 
  #subset_samples(PrePostFVT != "Post-Final") %>% 
  tax_glom(taxrank = "genus") %>% 
  subset_samples(Group == 1) %>% 
  otu_table() %>% as.data.frame() %>% rownames_to_column("tax_id")
ps.relabund.pups.otutable <- left_join(ps.relabund.pups.otutable, emu_tax %>% select(tax_id, species))
rownames(ps.relabund.pups.otutable) <- paste(ps.relabund.pups.otutable$species, ps.relabund.pups.otutable$tax_id, sep = "-")
ps.relabund.pups.otutable <- ps.relabund.pups.otutable %>% select(-c("tax_id", "species"))


metadata.nofinal <- metadata %>% 
  mutate(PrePostFVT = recode(PrePostFVT, "Post-Final" = "Post"))

set.seed(123)
fit_out <- maaslin3(input_data = ps.relabund.pups.otutable,
                    input_metadata = metadata.nofinal,
                    output = 'maaslin3-out-16s-group1',
                    formula = '~  PrePostFVT + (1 | AnimalID)',
                    #formula = '~ PrePostFVT + (1 | AnimalID)',
                    normalization = 'NONE',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 8,
                    save_models = FALSE)




maaslin.results.g1 <- read_tsv("maaslin3-out-16s-group1/all_results.tsv")
maaslin.results.g2 <- read_tsv("maaslin3-out-16s-group2/all_results.tsv")

maaslin.results.g1.signif <- maaslin.results.g1 %>% 
  filter(model == "abundance") %>% 
  filter(qval_individual <= 0.1) %>%
  #filter(grepl(":", name)) %>%
  filter(name == "PrePostFVTPost") %>% 
  separate_wider_delim(feature, "-", names = c("feature", "tax_id")) %>% 
  mutate(fold_change = 2^coef) %>% 
  mutate(type = case_when(coef > 0 ~ "up",
                          coef < 0 ~ "down",
                          TRUE ~ "ns"))


maaslin.results.g2.signif <- maaslin.results.g2 %>% 
  filter(model == "abundance") %>% 
  filter(qval_individual <= 0.1) %>%
  #filter(grepl(":", name)) %>%
  filter(name == "PrePostFVTPost") %>% 
  separate_wider_delim(feature, "-", names = c("feature", "tax_id")) %>% 
  mutate(fold_change = 2^coef)



maaslin.results.g1.signif %>% 
  filter(tax_id %in% setdiff(maaslin.results.g1.signif$tax_id, maaslin.results.g2.signif$tax_id))

#no taxa in group1 who are DA. only one was, but it was also DA in group2...