library(tidyverse)
library(phyloseq)
library(reshape2)
library(vegan)
library(ggrepel)
library(ggtext)

phage_cov <- read_tsv("analysis-exp3/phage-coverage.tsv")
phage_cov_filt <- phage_cov %>% 
  filter(meandepth > 1, coverage > 75) %>% 
  filter(endpos >= 1000)

phage_cov.summary <- phage_cov_filt %>%
  group_by(sample, contig) %>%
  summarize(phage_cov = meandepth)

phage_cov.totals <- phage_cov.summary %>% 
  group_by(sample) %>% 
  summarize(sample_total_cov = sum(phage_cov)) # get total coverage per sample

phage_cov.joined <- inner_join(phage_cov.summary, phage_cov.totals, by = c("sample" = "sample"))

phage_cov.joined$relabund_cov <- (phage_cov.joined$phage_cov / phage_cov.joined$sample_total_cov)


fvt.dates <- data.frame(
  "Dam" = c("F14", "F12", "F2"),
  "FVT.date" = c("FEB19", "FEB19", "FEB11")
) %>% 
  mutate(FVT.date = lubridate::mdy(paste0(FVT.date, "/2025")))


sample.depth <- read_csv("analysis-exp3/vivo8-vlp-depth.csv") %>% 
  group_by(sample) %>% 
  summarize(ReadDepth = sum(sum_len))
colnames(sample.depth) <- c("SampleName", "ReadDepth")

phage_metadata <- read.csv("analysis-exp3/vlp-metadata.csv")

phage_metadata$Group <- factor(phage_metadata$Group)
phage_metadata$PrePostFVT <- factor(phage_metadata$PrePostFVT, levels = c("Pre", "During", "Post"))

phage_metadata <- phage_metadata %>% 
  mutate(Date.formatted = if_else(grepl("DEC", Date, fixed = T), 
                                  lubridate::mdy(paste0(Date, "/2024")), 
                                  lubridate::mdy(paste0(Date, "/2025")))) %>% 
  left_join(fvt.dates, by = c("Dam" = "Dam")) %>% 
  mutate(RelativeDay = as.numeric(Date.formatted - FVT.date)) %>% 
  left_join(sample.depth)

rownames(phage_metadata) <- phage_metadata$SampleName

#phage_metadata$DateFormatted <- as.Date(paste0(phage_metadata$Date, "-2024"), format = "%b%d-%Y")
#phage_metadata$SampleName <- factor(phage_metadata$SampleName, levels = unique(phage_metadata$SampleName[order(phage_metadata$DateFormatted)]))

phage_cov.joined <- left_join(phage_cov.joined, phage_metadata, by = c("sample" = "SampleName"))



votu.table <- as.data.frame(pivot_wider(phage_cov.joined[, c("sample", "contig", "phage_cov")], id_cols = "contig", names_from = "sample", values_from = "phage_cov"))
rownames(votu.table) <- votu.table$contig
votu.table <- votu.table[, -1]
votu.table <- floor(votu.table)
votu.table[is.na(votu.table)] <- 0

PHAGE_OTU <- otu_table(votu.table, taxa_are_rows = T)

PHAGE_METADATA <- sample_data(phage_metadata)

ps.phage <- phyloseq(PHAGE_OTU, PHAGE_METADATA)
ps.phage.relabund <- transform_sample_counts(ps.phage, function(x) x / sum(x) )


ps.phage.relabund.pups <- ps.phage.relabund %>% 
  subset_samples(SampleType == "Pup")

set.seed(123)
ps.phage.relabund.pups.prepost <- ps.phage.relabund.pups %>% subset_samples(PrePostFVT %in% c("Pre", "Post"))
ord.nmds.bray <- ordinate(ps.phage.relabund.pups.prepost, method="NMDS", distance="bray")
fig.s2f <- plot_ordination(ps.phage.relabund.pups, ord.nmds.bray, color="Group", title="") +
  stat_ellipse(geom="polygon", type = "t", alpha = 0.1, aes(fill = Group)) +
  theme_bw() +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "TI-FVT", "PBS")) +
  scale_fill_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "TI-FVT", "PBS")) +
  facet_wrap(~PrePostFVT) +
  theme(aspect.ratio = 1, legend.position = "bottom") +
  guides(fill = "none") +
  ggtitle("Virome beta diversity")

fig.s2f


# PERMANOVA: did pups shift from baseline?
dist.vlp.bray <- distance(ps.phage.relabund.pups.prepost %>% subset_samples(Group == 1), method = "bray")
dist.vlp.bray.metadata <- phage_metadata %>% filter(SampleName %in% (dist.vlp.bray %>% as.matrix() %>% colnames()))
adonis2(dist.vlp.bray ~ PrePostFVT, 
        data = dist.vlp.bray.metadata, 
        strata = dist.vlp.bray.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, PrePostFVT p = 0.22, R2 = 0.071

# PERMANOVA: are live FVT pups distinct from HK-FVT/rPBS after treatment?
dist.vlp.bray.final <- distance(ps.phage.relabund.pups.prepost %>% subset_samples(PrePostFVT == "Post"), method = "bray")
dist.vlp.bray.final.metadata <- phage_metadata %>% filter(SampleName %in% (dist.vlp.bray.final %>% as.matrix() %>% colnames()))
adonis2(dist.vlp.bray.final ~ Group + PostnatalDay, 
        data = dist.vlp.bray.final.metadata, 
        strata = dist.vlp.bray.final.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, Group p = 0.13, PostnatalDay p = 0.13



phage.richness <- estimate_richness(ps.phage, measures=c("Shannon", "Observed")) %>% 
  rownames_to_column("sample") %>% 
  mutate(sample = gsub("\\.", '-', sample)) %>% 
  left_join(phage_metadata, by = c("sample" = "SampleName")) #%>% 
  #filter(SampleType == "Pup")

# modelling of pup richness
library(mgcv)
library(emmeans)

richness.data <- 
  phage.richness %>% 
  filter(SampleType == "Pup") %>% 
  mutate(AnimalID = as.factor(AnimalID))

model.richness.gam <- gam(Observed ~ Group + log10(ReadDepth) + s(RelativeDay, by=Group) + s(AnimalID, bs = "re"), 
                          data = richness.data, family = gaussian(), method = "REML")

model.shannon.gam <- gam(Shannon ~ Group +  log10(ReadDepth) + s(RelativeDay, by=Group) + s(AnimalID, bs = "re"), 
                         data = richness.data, family = gaussian(), method = "REML")


summary(model.richness.gam)
draw(model.richness.gam)

library(gratia)

fig.s2d2 <- model.richness.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(metric = "richness") %>% 
  rbind(
    model.shannon.gam %>% 
      smooth_estimates() %>% 
      add_confint() %>% 
      mutate(metric = "shannon")
  ) %>% 
  filter(.smooth %in% c("s(RelativeDay):Group1", "s(RelativeDay):Group2", "s(RelativeDay):Group3")) %>% 
  ggplot(aes(x = RelativeDay, y = .estimate, color = .smooth)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth), linewidth = 0.1, alpha = 0.1) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "TI-FVT", "PBS")) +
  scale_fill_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "TI-FVT", "PBS")) +
  theme_bw() + 
  labs(title = "Virome alpha diversity",
       x = "Days relative to first FVT dose",
       y = "Partial effect") +
  theme(legend.title = element_blank(),
        legend.key.spacing.y = unit(0.15, 'cm'),
        legend.position = "bottom",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold")
  ) +
  scale_x_continuous(breaks = seq(-14, 45, by = 7)) +
  facet_wrap(~metric, scale = "free", ncol = 1, labeller = labeller(metric = c("richness" = "Richness", "shannon" = "Shannon")))


fig.s2d2



### FIND NEW vOTUS IN 'DURING' VS 'PRE' TIMEPOINT
new_votus_by_group <- ps.phage.relabund.pups %>% 
  psmelt() %>% 
  filter(PrePostFVT %in% c("Pre", "During")) %>% 
  filter(Abundance > 0) %>% 
  distinct(AnimalID, Group, PrePostFVT, OTU) %>% 
  group_by(AnimalID, Group) %>% 
  mutate(new_votu = OTU %in% OTU[PrePostFVT == "During"] & !(OTU %in% OTU[PrePostFVT == "Pre"])) %>% 
  ungroup() %>% 
  filter(new_votu, PrePostFVT == "During") %>% 
  count(Group, OTU, sort = TRUE, name = "n_animals") %>% 
  arrange(Group, desc(n_animals))
  

new_votus_by_group %>% 
  filter(n_animals > 2) %>% 
  ggplot(aes(x = OTU, y = n_animals)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Group, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ps.phage.relabund.pups %>% 
  psmelt() %>% 
  filter(Group == 1) %>% 
  filter(OTU %in% (new_votus_by_group %>% filter(n_animals > 2, Group == 1) %>% pull(OTU))) %>% 
  ggplot(aes(x = PostnatalDay, y = Abundance, group = AnimalID, color=PrePostFVT)) +
  geom_point() +
  geom_line()+
  facet_wrap(~OTU, scales = "free")



### look at votus from the FVT in the pups
fvt.votus <- ps.phage.relabund %>% 
  psmelt() %>% 
  filter(SampleName == "VIVO8_VLP_DONOR2") %>% 
  filter(Abundance > 0) %>% 
  pull(OTU)


ps.phage.relabund.pups %>% 
  psmelt() %>% 
  filter(Group == 1) %>% 
  filter(OTU %in% fvt.votus) %>% 
  #filter(Abundance > 0) %>% 
  filter(OTU %in% c("NODE_912_length_8852_cov_15.026486", "NODE_994_length_4515_cov_8.534753")) %>% 
  ggplot(aes(x = PostnatalDay, y = Abundance, group = AnimalID, color=PrePostFVT)) +
  geom_point() +
  geom_line()+
  facet_wrap(~OTU, scales = "free")



### FVT engraftment analysis ###
fvt.votus <- ps.phage.relabund %>% 
  psmelt() %>% 
  filter(SampleName == "VIVO8_VLP_DONOR2", Abundance > 0) %>% 
  pull(OTU)

ps.phage.relabund %>% 
  psmelt() %>% 
  filter(Group %in% c(1), Abundance > 0, SampleType == "Pup") %>% 
  filter(OTU %in% fvt.votus) %>% 
  group_by(OTU, PrePostFVT, Group) %>% 
  summarise(n = n()) %>% view()



### DISTANCES FROM PUPS TO FVT (H DONOR INOCULUM) ###
ps.phage.relabund.pups.fvt <- ps.phage.relabund %>% 
  subset_samples(SampleType %in% c("Inoculum", "Pup")) %>% 
  subset_samples(SampleName != "VIVO8_VLP_DONOR1")

pups.fvt.dist <- ps.phage.relabund.pups.fvt %>% 
  distance(method = "bray")

pups.fvt.dist.metadata <- data.frame(sample_data(ps.phage.relabund.pups.fvt))

pups.fvt.dist %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName") %>% 
  melt() %>% 
  filter(value > 0) %>% 
  filter(variable == "VIVO8_VLP_DONOR2") %>% 
  left_join(pups.fvt.dist.metadata) %>% 
  ggplot(aes(x = RelativeDay, y = value)) +
  geom_point() +
  facet_wrap(~Group) +
  ylab("Distance to FVT")


#### distance from first sample (pre-FVT) ####
dist.pups.bray <-  as.matrix(distance(ps.phage.relabund %>% subset_samples(SampleType == "Pup"), method = "bray")) %>% 
  melt() %>% 
  filter(Var1 != Var2)

colnames(dist.pups.bray) <- c("sample1", "sample2", "dist")

dist.pups.bray.metadata <- dist.pups.bray %>% 
  left_join(phage_metadata, by = c("sample1" = "SampleName")) %>% 
  left_join(phage_metadata, by = c("sample2" = "SampleName"))

dist.pups.bray.metadata.samepup <- dist.pups.bray.metadata %>% 
  filter(AnimalID.x == AnimalID.y)

dist.pups.bray.metadata.samepup.baseline <- dist.pups.bray.metadata.samepup %>% 
  group_by(AnimalID.x) %>% 
  mutate(first_timepoint = min(PostnatalDay.x)) %>% 
  filter(PostnatalDay.x == first_timepoint) %>% 
  filter(PostnatalDay.y != first_timepoint) %>% 
  mutate(timepoint_diff = PostnatalDay.y - first_timepoint)


dist.pups.bray.metadata.samepup.baseline %>% 
  ggplot(aes(x = RelativeDay.y, y = dist, group = AnimalID.x, color = Dam.x)) +
  #ggplot(aes(x = PrePostFVT.y, y = dist, group = interaction(PrePostFVT.y,Dam.x, Group.x), color = PrePostFVT.y)) +
  geom_line(color = "darkgray") +
  geom_point() +
  #geom_boxplot() +
  facet_wrap(~Group.x*Dam.x, ncol=3, labeller = as_labeller(c("1" = "Live FVT", "2" = "hkFVT", "3" = "rPBS", "F2"="F2", "F14"="F14", "F12"="F12"))) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Postnatal Day", y = "Distance to baseline") +
  theme(aspect.ratio = 1)

dist.pups.bray.metadata.samepup.baseline %>% 
  slice_max(PostnatalDay.y, n = 1, with_ties = FALSE) %>%
  filter(PostnatalDay.y > 70) %>% 
  mutate(Group.x = fct_recode(Group.x, "Live FVT" = "1", "hkFVT" = "2", "rPBS" = "3")) %>% 
  ggplot(aes(x = Group.x, y = dist)) +
  geom_boxplot() +
  geom_jitter(aes(color = Dam.x), position = position_dodge2(width = 0.5)) +
  theme_bw() +
  labs(x = "Group", y = "Distance to baseline", subtitle = "VIVO8 - Virome") +
  theme(aspect.ratio = 1)







