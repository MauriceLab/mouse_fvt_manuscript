library(maaslin3)
library(tidyverse)
library(phyloseq)
library(ggrepel)
library(ggtext)

ps.relabund.pups.otutable <- ps.relabund.pups %>% 
  subset_samples(Group %in% c(2)) %>% # repeat with group set to 1, 2
  otu_table() %>% as.data.frame() %>% rownames_to_column("tax_id")
ps.relabund.pups.otutable <- left_join(ps.relabund.pups.otutable, emu_tax %>% select(tax_id, species))
rownames(ps.relabund.pups.otutable) <- paste(ps.relabund.pups.otutable$species, ps.relabund.pups.otutable$tax_id, sep = "-")
ps.relabund.pups.otutable <- ps.relabund.pups.otutable %>% select(-c("tax_id", "species"))

fit_out <- maaslin3(input_data = ps.relabund.pups.otutable,
                    input_metadata = metadata,
                    output = 'maaslin3-out-16s-group2',
                    formula = '~ PrePostFVT + (1 | AnimalID)',
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

maaslin.results.g1 <- read_tsv("analysis-exp2/maaslin3-out-16s-group1/all_results.tsv")
maaslin.results.g2 <- read_tsv("analysis-exp2/maaslin3-out-16s-group2/all_results.tsv")

maaslin.results.g2.signif <- maaslin.results.g2 %>% 
  filter(model == "abundance") %>% 
  filter(qval_individual <= 0.1) %>%
  filter(name == "PrePostFVTPost") %>% 
  separate_wider_delim(feature, "-", names = c("feature", "tax_id")) %>% 
  mutate(fold_change = 2^coef) %>% 
  pull(tax_id)

'%!in%' <- function(x,y)!('%in%'(x,y))
maaslin.results.g1.signif <- maaslin.results.g1 %>% 
  filter(model == "abundance") %>% 
  filter(name == "PrePostFVTPost") %>% 
  separate_wider_delim(feature, "-", names = c("feature", "tax_id")) %>% 
  mutate(signif = ifelse((tax_id %!in% maaslin.results.g2.signif) & qval_individual <= 0.1, T, F)) %>% 
  mutate(type = case_when(coef > 0 & signif ~ "up",
                          coef < 0 & signif ~ "down",
                          TRUE ~ "ns"))


cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

set.seed(123)
fig.3c1 <- maaslin.results.g1.signif %>% 
  ggplot(aes(x = coef, y = -log10(qval_individual), label = feature, colour = type)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(-4, 4, 2), limits = c(-4, 4)) +
  scale_y_continuous(breaks = seq(0, 5, 2.5), limits = c(0, 5)) +
  xlab("log<sub>2</sub>(FoldChange)") +
  ylab("-log<sub>10</sub>(P-value)") +
  scale_colour_manual(values = cols) +
  geom_label_repel(data = maaslin.results.g1.signif %>%
                    filter(type %in% c("up", "down")), aes(label = feature), size = 1.75, fill = "white", colour = "black", min.segment.length = 0) +
  theme_bw() +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        #panel.border = element_rect(colour = "black", fill = NA, size= 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        aspect.ratio = 1, 
        legend.position = "none") +
  ggtitle("Bacteriome - Experiment #2")


fig.3c1

fig.s3a <- ps.relabund.pups %>% 
  psmelt() %>% 
  filter(OTU %in% (maaslin.results.g1.signif %>% filter(signif))$tax_id) %>% 
  filter(Group %in% c(1, 2)) %>% 
  ggplot(aes(x = PrePostFVT, y = Abundance, group = interaction(Dam, Group, PrePostFVT), color = Dam)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ species, ncol = 2, scales="free_y") +
  scale_color_manual(values = c("#00A087FF", "#006622", "#E64B35FF"), name = "Dam", labels = c("F2\n(Live FVT)", "F5\n(Live FVT)", "F8\n(TI-FVT)")) +
  scale_y_continuous(labels=scales::percent) +
  xlab("Timepoint (pre/during/post FVT)") +
  ylab("Relative abundance (%)") +
  ggtitle("Bacteriome - Experiment #2")


fig.s3a
