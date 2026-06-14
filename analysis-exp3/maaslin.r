library(maaslin3)
library(phyloseq)
library(tidyverse)
library(ggrepel)
library(ggtext)

### votu level ###
phage_metadata <- phage_metadata %>% 
  mutate(PrePostFVT = recode(PrePostFVT, "Post-Final" = "Post"))

sample_data(ps.phage.relabund) <- phage_metadata %>% 
  mutate(PrePostFVT = recode(PrePostFVT, "Post-Final" = "Post")) %>%
  sample_data()

ps.phage.relabund.otutable <- ps.phage.relabund %>% 
  subset_samples(Group %in% c(3)) %>% # repeat for groups 1, 2, 3
  otu_table() %>% 
  as.data.frame()

set.seed(123)
fit_out <- maaslin3(input_data = ps.phage.relabund.otutable,
                    input_metadata = phage_metadata,
                    output = 'maaslin3-out-vlp-group3',
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


# do volcano plot for G1 DA vOTUs

maaslin.results.g1 <- read_tsv("analysis-exp3/maaslin3-out-vlp-group1/all_results.tsv")
maaslin.results.g2 <- read_tsv("analysis-exp3/maaslin3-out-vlp-group2/all_results.tsv")
maaslin.results.g3 <- read_tsv("analysis-exp3/maaslin3-out-vlp-group3/all_results.tsv")


maaslin.results.g2.signif <- maaslin.results.g2 %>% 
  filter(model == "abundance") %>% 
  filter(qval_individual <= 0.1) %>%
  filter(name == "PrePostFVTPost") %>% 
  mutate(fold_change = 2^coef) %>% 
  pull(feature)

maaslin.results.g3.signif <- maaslin.results.g3 %>% 
  filter(model == "abundance") %>% 
  filter(qval_individual <= 0.1) %>%
  filter(name == "PrePostFVTPost") %>% 
  mutate(fold_change = 2^coef) %>% 
  pull(feature)

'%!in%' <- function(x,y)!('%in%'(x,y))
maaslin.results.g1.signif <- maaslin.results.g1 %>% 
  filter(model == "abundance") %>% 
  filter(name == "PrePostFVTPost") %>% 
  mutate(signif = ifelse((feature %!in% maaslin.results.g2.signif) & (feature %!in% maaslin.results.g3.signif) & qval_individual <= 0.1, T, F)) %>% 
  mutate(type = case_when(coef > 0 & signif ~ "up",
                          coef < 0 & signif ~ "down",
                          TRUE ~ "ns"))




cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

set.seed(123)
fig.3d <- maaslin.results.g1.signif %>% 
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
  ggtitle("Virome - Experiment #3")

fig.3d


fig.s3c <- ps.phage.relabund %>% 
  subset_samples(SampleType == "Pup") %>% 
  psmelt() %>% 
  mutate(PrePostFVT = recode(PrePostFVT, "Post-Final" = "Post")) %>% 
  #filter(OTU %in% maaslin.results.g2.signif$tax_id) %>% 
  filter(OTU %in% (maaslin.results.g1.signif %>% filter(signif))$feature) %>% 
  filter(Group %in% c(1,2,3)) %>% 
  ggplot(aes(x = PrePostFVT, y = Abundance, group = interaction(Group, PrePostFVT), color = Group)) +
  #geom_point() +
  geom_boxplot() +
  facet_wrap(~ OTU, scales="free", nrow = 2) +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), name = "Group", labels = c("Live FVT", "TI-FVT", "rPBS")) +
  theme_bw() +
  scale_y_continuous(labels=scales::percent) +
  xlab("Timepoint (pre/during/post FVT)") +
  ylab("Relative abundance (%)") +
  ggtitle("Virome - Experiment #3")

fig.s3c


