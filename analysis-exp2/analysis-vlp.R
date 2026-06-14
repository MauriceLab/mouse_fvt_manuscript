library(tidyverse)
library(phyloseq)
library(reshape2)
library(vegan)
library(ggrepel)
library(ggtext)

phage_cov <- read_tsv("analysis-exp2/phage-coverage.tsv")
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
  "Dam" = c("F2", "F5", "F8", "F11"),
  "FVT.date" = c("MAY28", "JUN11", "JUN19", "MAY28")
) %>% 
  mutate(FVT.date = lubridate::mdy(paste0(FVT.date, "/2024")))

sample.depth <- read_csv("analysis-exp2/vivo7-vlp-depth.csv") %>% 
  group_by(SampleName) %>% 
  summarize(ReadDepth = sum(sum_len))
colnames(sample.depth) <- c("SampleName", "ReadDepth")

phage_metadata <- read.csv("analysis-exp2/vlp-metadata.csv")

phage_metadata$Group <- factor(phage_metadata$Group)
phage_metadata$PrePostFVT <- factor(phage_metadata$PrePostFVT, levels = c("Pre", "During", "Post"))
phage_metadata$DateFormatted <- as.Date(paste0(phage_metadata$Date, "-2024"), format = "%b%d-%Y")
phage_metadata$SampleName <- factor(phage_metadata$SampleName, levels = unique(phage_metadata$SampleName[order(phage_metadata$DateFormatted)]))

phage_metadata <- phage_metadata %>% 
  mutate(Date.formatted = lubridate::mdy(paste0(Date, "/2024"))) %>% 
  left_join(fvt.dates, by = c("Dam" = "Dam")) %>% 
  mutate(RelativeDay = as.numeric(Date.formatted - FVT.date)) %>% 
  left_join(sample.depth)

rownames(phage_metadata) <- phage_metadata$SampleName

phage_cov.joined <- left_join(phage_cov.joined, phage_metadata, by = c("sample" = "SampleName"))

#=========

votu.table <- as.data.frame(pivot_wider(phage_cov.joined[, c("sample", "contig", "phage_cov")], id_cols = "contig", names_from = "sample", values_from = "phage_cov"))
rownames(votu.table) <- votu.table$contig
votu.table <- votu.table[, -1]
votu.table <- floor(votu.table)
votu.table[is.na(votu.table)] <- 0

PHAGE_OTU <- otu_table(votu.table, taxa_are_rows = T)

PHAGE_METADATA <- sample_data(phage_metadata)

ps.phage <- phyloseq(PHAGE_OTU, PHAGE_METADATA)
ps.phage.relabund <- transform_sample_counts(ps.phage, function(x) x / sum(x) )


ps.phage.relabund.pups.fvtgroups <- ps.phage.relabund %>% 
  subset_samples(SampleType == "Pup") %>% 
  subset_samples(Group %in% c(1, 2)) %>% 
  subset_samples(PrePostFVT %in% c("Pre", "Post")) #%>% 
  #subset_samples(Dam != "F5") # %>% 
  #subset_samples(Dam %in% c("F2","F8","F11"))

set.seed(123)
ord.nmds.bray <- ordinate(ps.phage.relabund.pups.fvtgroups, method="NMDS", distance="bray")
fig.s2c <- plot_ordination(ps.phage.relabund.pups.fvtgroups, ord.nmds.bray, color="Dam", title="") +
  stat_ellipse(geom="polygon", type = "t", alpha = 0.1, aes(fill = Dam)) +
  theme_bw() +
  scale_color_manual(values = c("#00A087FF", "#006622", "#E64B35FF"), name = "Dam", labels = c("F2\n(Live FVT)", "F5\n(Live FVT)", "F8\n(TI-FVT)")) +
  scale_fill_manual(values = c("#00A087FF", "#006622", "#E64B35FF"), name = "Dam", labels = c("F2\n(Live FVT)", "F5\n(Live FVT)", "F8\n(TI-FVT)")) +
  facet_wrap(~PrePostFVT) +
  guides(fill = "none") +
  theme(legend.position = "bottom") +
  ggtitle("Virome beta diversity")

fig.s2c

# PERMANOVA: did pups shift from baseline?
dist.vlp.bray <- distance(ps.phage.relabund.pups.fvtgroups %>% subset_samples(Group == 1), method = "bray")
dist.vlp.bray.metadata <- phage_metadata %>% filter(SampleName %in% (dist.vlp.bray %>% as.matrix() %>% colnames()))
adonis2(dist.vlp.bray ~ PrePostFVT, 
        data = dist.vlp.bray.metadata, 
        strata = dist.vlp.bray.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, PrePostFVT p = 0.029, R2 = 0.043

# PERMANOVA: are live FVT pups distinct from HK-FVT after treatment?
dist.vlp.bray.final <- distance(ps.phage.relabund.pups.fvtgroups %>% subset_samples(PrePostFVT == "Post"), method = "bray")
dist.vlp.bray.final.metadata <- phage_metadata %>% filter(SampleName %in% (dist.vlp.bray.final %>% as.matrix() %>% colnames()))
adonis2(dist.vlp.bray.final ~ Group + PostnatalDay, 
        data = dist.vlp.bray.final.metadata, 
        strata = dist.vlp.bray.final.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, Group p = 0.64, PostnatalDay p = 0.64




### alpha diversity
phage.richness <- estimate_richness(ps.phage, measures=c("Shannon", "Observed")) %>% 
  rownames_to_column("sample") %>% 
  mutate(sample = gsub("\\.", '-', sample)) %>% 
  left_join(phage_metadata, by = c("sample" = "SampleName"))


# modelling of pup richness
library(mgcv)
library(emmeans)
library(gratia)

richness.data <- 
  phage.richness %>% 
  filter(SampleType == "Pup") %>% 
  mutate(AnimalID = as.factor(AnimalID), Dam = as.factor(Dam))

model.richness.gam <- gam(Observed ~ Group + log10(ReadDepth) + s(RelativeDay, by=Group) + s(AnimalID, bs = "re") + s(Dam, bs = "re"), 
                          data = richness.data, family = gaussian(), method = "REML")

model.shannon.gam <- gam(Shannon ~ Group +  log10(ReadDepth) + s(RelativeDay, by=Group) + s(AnimalID, bs = "re") + s(Dam, bs = "re"), 
                         data = richness.data, family = gaussian(), method = "REML")


summary(model.richness.gam)
draw(model.richness.gam)

fig.s2a2 <- model.richness.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(metric = "richness") %>% 
  rbind(
    model.shannon.gam %>% 
      smooth_estimates() %>% 
      add_confint() %>% 
      mutate(metric = "shannon")
  ) %>% 
  #filter(.smooth %in% c("s(PostnatalDay):Group1", "s(PostnatalDay):Group2", "s(PostnatalDay):Group4")) %>% 
  filter(.smooth %in% c("s(RelativeDay):Group1", "s(RelativeDay):Group2")) %>% 
  ggplot(aes(x = RelativeDay, y = .estimate, color = .smooth)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth), linewidth = 0.1, alpha = 0.1) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF"), labels = c("Live FVT", "TI-FVT")) +
  scale_fill_manual(values = c("#00A087FF", "#E64B35FF"), labels = c("Live FVT", "TI-FVT")) +
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

fig.s2a2

ds.richness <- data_slice(model.richness.gam, 
                          RelativeDay = seq(-15, 45, length.out = 100), 
                          Group = unique(richness.data$Group))
fv.richness <- fitted_values(model.richness.gam, data = ds.richness, scale = "response") %>% mutate(metric = "richness")

ds.shannon <- data_slice(model.shannon.gam, 
                         RelativeDay = seq(-15, 45, length.out = 100), 
                         Group = unique(richness.data$Group))
fv.shannon <- fitted_values(model.shannon.gam, data = ds.shannon, scale = "response") %>% mutate(metric = "shannon")



#### FVT engraftment analysis ####
fvt.votus <- ps.phage.relabund %>% 
  psmelt() %>% 
  filter(SampleName == "VIVO7-DONOR2", Abundance > 0) %>% 
  pull(OTU)

ps.phage.relabund %>% 
  psmelt() %>% 
  filter(Group %in% c(1, 2), Abundance > 0) %>% 
  filter(OTU %in% fvt.votus) %>% 
  group_by(OTU, PrePostFVT, Group) %>% 
  summarise(n = n())



### DISTANCES FROM PUPS TO FVT (H DONOR INOCULUM) ###
ps.phage.relabund.pups.fvt <- ps.phage.relabund %>% 
  subset_samples(SampleType %in% c("Inoculum", "Pup")) %>% 
  subset_samples(SampleName != "VIVO7-DONOR1")

pups.fvt.dist <- ps.phage.relabund.pups.fvt %>% 
  distance(method = "bray")

pups.fvt.dist.metadata <- data.frame(sample_data(ps.phage.relabund.pups.fvt))

pups.fvt.dist %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName") %>% 
  melt() %>% 
  filter(value > 0) %>% 
  filter(variable == "VIVO7-DONOR2") %>% 
  left_join(pups.fvt.dist.metadata) %>% 
  ggplot(aes(x = RelativeDay, y = value, colour = Dam)) +
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
  facet_wrap(~Group.x*Dam.x, ncol=4, labeller = as_labeller(c("1" = "Live FVT", "2" = "hkFVT", "4" = "H control/PBS", "F2"="F2", "F5"="F5", "F8"="F8", "F11"="F11"))) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Postnatal Day", y = "Distance to baseline") +
  theme(aspect.ratio = 1)



