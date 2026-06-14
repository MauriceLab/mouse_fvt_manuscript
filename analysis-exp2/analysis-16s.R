library(tidyverse)
library(phyloseq)
library(vegan)
library(reshape2)
library(microshades)
library(cowplot)
library(grid)
library(gridExtra)
library(ggtext)
library(ggrepel)

emu_tax <- read_tsv("analysis-exp2/emu-combined-taxonomy-tax_id.tsv") %>% filter(tax_id != "unassigned")
emu_tax.mat <- as.matrix(emu_tax[, -1])
rownames(emu_tax.mat) <- emu_tax$tax_id
emu_tax.mat <- emu_tax.mat[, rev(seq_len(ncol(emu_tax.mat)))]


emu_counts <- read_tsv("analysis-exp2/emu-combined-abundance-tax_id-counts.tsv") %>% filter(tax_id != "unassigned")
emu_counts.mat <- as.matrix(emu_counts[, -1])
rownames(emu_counts.mat) <- emu_counts$tax_id
colnames(emu_counts.mat) <- gsub(".fq", "", colnames(emu_counts.mat))
emu_counts.mat <- floor(emu_counts.mat)

### Make Phyloseq object
TAX <- tax_table(emu_tax.mat)
emu_counts.mat.zeros <- emu_counts.mat
emu_counts.mat.zeros[is.na(emu_counts.mat.zeros)] <- 0
OTU <- otu_table(emu_counts.mat.zeros, taxa_are_rows = T,)
TREE <- read_tree("emu_tax_tree_species.nwk")


fvt.dates <- data.frame(
  "Dam" = c("F2", "F5", "F8", "F11"),
  "FVT.date" = c("MAY28", "JUN11", "JUN19", "MAY28")
  ) %>% 
  mutate(FVT.date = lubridate::mdy(paste0(FVT.date, "/2024")))

sample.depth <- colSums(emu_counts.mat, na.rm = T) %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName")
colnames(sample.depth) <- c("SampleName", "ReadDepth")

metadata <- read.csv("analysis-exp2/16s-metadata.csv")
metadata$Group <- factor(metadata$Group)
metadata$DonorID <- factor(metadata$DonorID)
metadata$AnimalID <- factor(metadata$AnimalID)
metadata$SampleType <- factor(metadata$SampleType)
metadata$Dam <- factor(metadata$Dam)
metadata$PrePostFVT <- factor(metadata$PrePostFVT, levels = c("Pre", "During", "Post"))

metadata <- metadata %>% 
  mutate(Date.formatted = lubridate::mdy(paste0(Date, "/2024"))) %>% 
  left_join(fvt.dates, by = c("Dam" = "Dam")) %>% 
  mutate(RelativeDay = as.numeric(Date.formatted - FVT.date)) %>% 
  left_join(sample.depth)

rownames(metadata) <- metadata$SampleName

METADATA <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, TREE, METADATA)
ps.relabund <- transform_sample_counts(ps, function(x) x/sum(x))

ps.relabund.pups <- ps.relabund %>% subset_samples(SampleType == "Pup")


#pcoa
#pcoa by timepoint
ord.pcoa.uni <- ordinate(ps.relabund.pups %>% 
                           subset_samples(PrePostFVT %in% c("Pre", "Post")) %>% 
                           subset_samples(Group %in% c(1, 2)), method="PCoA", distance="wunifrac")
fig.s2b <- plot_ordination(ps.relabund.pups, ord.pcoa.uni, color="Dam", title="") +
  stat_ellipse(geom="polygon",type = "t", alpha = 0.1, aes(fill=Dam)) +
  theme_bw() +
  scale_color_manual(values = c("#00A087FF", "#006622", "#E64B35FF"), name = "Dam", labels = c("F2\n(Live FVT)", "F5\n(Live FVT)", "F8\n(TI-FVT)")) +
  scale_fill_manual(values = c("#00A087FF", "#006622", "#E64B35FF"), name = "Dam", labels = c("F2\n(Live FVT)", "F5\n(Live FVT)", "F8\n(TI-FVT)")) +
  facet_wrap(~PrePostFVT) +
  theme(legend.position = "bottom") +
  labs(col = "Dam") +
  guides(fill = "none", shape = "none") +
  xlab("Axis.1 [45.8%]") +
  ylab("Axis.2 [36.0%]") +
  ggtitle("Bacteriome beta diversity")

fig.s2b

metadata.pups <- sample_data(ps.relabund.pups)
class(metadata.pups) <- "data.frame"
dist.wuni <- distance(ps.relabund.pups, method = "wunifrac")
set.seed(123)
permanova <- adonis2(dist.wuni ~ PrePostFVT*Group*AnimalID*Dam*PostnatalDay, data = metadata.pups, by = "terms")
permanova

# PERMANOVA: did pups shift from baseline?
dist.wuni <- distance(ps.relabund.pups %>% subset_samples(Group %in% c(1, 2)) %>% subset_samples(Group == 1), method = "wunifrac")
dist.wuni.metadata <- metadata %>% filter(SampleName %in% (dist.wuni %>% as.matrix() %>% colnames()))
adonis2(dist.wuni ~ PrePostFVT, 
        data = dist.wuni.metadata, 
        strata = dist.wuni.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, PrePostFVT p = 0.302, R2 = 0.03



# PERMANOVA: are live FVT pups distinct from HK-FVT after treatment?
dist.wuni.final <- distance(ps.relabund.pups %>% subset_samples(Group %in% c(1, 2)) %>% subset_samples(PrePostFVT == "Post"), method = "wunifrac")
dist.wuni.final.metadata <- metadata %>% filter(SampleName %in% (dist.wuni.final %>% as.matrix() %>% colnames()))
adonis2(dist.wuni.final ~ Group + PostnatalDay, 
        data = dist.wuni.final.metadata, 
        strata = dist.wuni.final.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, Group p = 0.172, PostnatalDay p = 0.162

pairwiseAdonis::pairwise.adonis2(dist.wuni.final ~ Dam,
                                 data = dist.wuni.final.metadata,
                                 strata = "AnimalID")



# alpha diversity
library(mgcv)
library(emmeans)
library(gratia)

richness <- estimate_richness(ps, measures=c("Shannon", "Observed")) %>% 
  rownames_to_column("sample") %>% 
  mutate(sample = gsub("\\.", '-', sample)) %>% 
  left_join(metadata, by = c("sample" = "SampleName"))

richness.data <- 
  richness %>% 
  filter(SampleType == "Pup") %>% 
  filter(Group %in% c(1, 2)) %>% 
  mutate(AnimalID = as.factor(AnimalID), Dam = as.factor(Dam))

model.richness.gam <- gam(Observed ~ Group + s(RelativeDay, by = Group) + s(AnimalID, bs = "re") + s(Dam, bs = "re"), 
                          data = richness.data, family = nb(), method = "REML")
summary(model.richness.gam)
draw(model.richness.gam)

model.shannon.gam <- gam(Shannon ~ Group + s(RelativeDay, by = Group) + s(AnimalID, bs = "re") + s(Dam, bs = "re"), 
                          data = richness.data, family = gaussian(), method = "REML")


fig.s2a1 <- model.richness.gam %>% 
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
  labs(title = "Bacteriome alpha diversity",
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

fig.s2a1

ds.richness <- data_slice(model.richness.gam, 
                          RelativeDay = seq(-15, 45, length.out = 100), 
                          Group = unique(richness.data$Group))
fv.richness <- fitted_values(model.richness.gam, data = ds.richness, scale = "response") %>% mutate(metric = "richness")

ds.shannon <- data_slice(model.shannon.gam, 
                          RelativeDay = seq(-15, 45, length.out = 100), 
                          Group = unique(richness.data$Group))
fv.shannon <- fitted_values(model.shannon.gam, data = ds.shannon, scale = "response") %>% mutate(metric = "shannon")


fv.richness %>% rbind(fv.shannon) %>% 
  filter(Group %in% c(1, 2)) %>% 
  ggplot(aes(x = RelativeDay, y = .fitted, color = Group)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = Group), linewidth = 0.1, alpha = 0.07) +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF"), labels = c("Live FVT", "TI-FVT")) +
  scale_fill_manual(values = c("#00A087FF", "#E64B35FF"), labels = c("Live FVT", "TI-FVT")) +
  theme_bw() + 
  labs(title = "Bacteriome alpha diversity",
       x = "Days relative to first FVT dose",
       y = "Fitted value\n(total effect)") +
  theme(legend.title = element_blank(),
        legend.key.spacing.y = unit(0.15, 'cm'),
        legend.position = "bottom"
  ) +
  facet_wrap(~ metric, scale = "free_y", labeller = as_labeller(c(richness = "Richness", shannon = "Shannon")), nrow = 2) +
  scale_x_continuous(breaks = seq(-14, 45, by = 7))

# distance from previous timepoint (stability)

dist.pups.wuni <-  as.matrix(distance(ps.relabund.pups, method = "wunifrac")) %>% 
  melt() %>% 
  filter(Var1 != Var2)

colnames(dist.pups.wuni) <- c("sample1", "sample2", "dist")

dist.pups.wuni.metadata <- dist.pups.wuni %>% 
  left_join(metadata.pups, by = c("sample1" = "SampleName")) %>% 
  left_join(metadata.pups, by = c("sample2" = "SampleName"))

dist.pups.wuni.metadata.samepup <- dist.pups.wuni.metadata %>% 
  filter(AnimalID.x == AnimalID.y)


dist.pups.wuni.metadata.samepup.stability <- 
  aggregate(dist ~ AnimalID.x + Group.x + Dam.x + PostnatalDay.x + PostnatalDay.y, 
            data = dist.pups.wuni.metadata.samepup, 
            FUN = mean) %>% 
  mutate(timepoint_after = PostnatalDay.y > PostnatalDay.x) %>% 
  filter(timepoint_after) %>% 
  mutate(timepoint_diff = PostnatalDay.y - PostnatalDay.x) %>% 
  group_by(AnimalID.x, PostnatalDay.x) %>% 
  arrange(timepoint_diff) %>% 
  slice(1)


dist.pups.wuni.metadata.samepup.stability %>% 
  ggplot(aes(x = PostnatalDay.y, y = 1 - dist, group = AnimalID.x, color = Dam.x)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Group.x) +
  theme_bw() +
  labs(x = "Postnatal Day", y = "Stability (1 - distance)") +
  theme(aspect.ratio = 1)


# distance from first sample (pre-FVT)
dist.pups.wuni.metadata.samepup.baseline <- dist.pups.wuni.metadata.samepup %>% 
  group_by(AnimalID.x) %>% 
  mutate(first_timepoint = min(PostnatalDay.x)) %>% 
  filter(PostnatalDay.x == first_timepoint) %>% 
  filter(PostnatalDay.y != first_timepoint) %>% 
  mutate(timepoint_diff = PostnatalDay.y - first_timepoint)


dist.pups.wuni.metadata.samepup.baseline %>% 
  ggplot(aes(x = PostnatalDay.y, y = dist, group = AnimalID.x, color = Dam.x)) +
  #ggplot(aes(x = PrePostFVT.y, y = dist, group = interaction(PrePostFVT.y,Dam.x, Group.x), color = PrePostFVT.y)) +
  geom_line(color = "darkgrey") +
  geom_point() +
  #geom_boxplot() +
  facet_wrap(~Group.x*Dam.x, ncol=4, labeller = as_labeller(c("1" = "Live FVT", "2" = "hkFVT", "4" = "H control/PBS", "F2"="F2", "F5"="F5", "F8"="F8", "F11"="F11"))) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Postnatal Day", y = "Distance to baseline") +
  theme(aspect.ratio = 1)

# LAST TIMEPOINT ONLY
dist.pups.wuni.metadata.samepup.baseline %>% 
  slice_max(PostnatalDay.y, n = 1, with_ties = FALSE) %>%
  filter(PostnatalDay.y > 70) %>% 
  mutate(Group.x = fct_recode(Group.x, "Live FVT" = "1", "hkFVT" = "2", "H control" = "4")) %>% 
  ggplot(aes(x = Group.x, y = dist)) +
  geom_boxplot() +
  geom_jitter(aes(color = Dam.x), position = position_dodge2(width = 0.5)) +
  theme_bw() +
  labs(x = "Group", y = "Distance to baseline", subtitle = "VIVO7 - Bacteriome") +
  theme(aspect.ratio = 1)


# pup distance from dam over time
dist.pups.to.dams.wuni <-  as.matrix(distance(ps.relabund, method = "wunifrac")) %>% 
  melt() %>% 
  filter(Var1 != Var2)

colnames(dist.pups.to.dams.wuni) <- c("sample1", "sample2", "dist")

dist.pups.to.dams.wuni.metadata <- dist.pups.to.dams.wuni %>% 
  left_join(metadata, by = c("sample1" = "SampleName")) %>% 
  left_join(metadata, by = c("sample2" = "SampleName")) %>% 
  filter(SampleType.x == "Pup") %>%
  filter(SampleType.y == "Breeder") %>% 
  filter(as.character(Dam.x) == as.character(AnimalID.y)) %>% 
  filter(Date.y %in% c("APR30", "MAY16", "MAY27", "MAY2")) #dam postpartum timepoint (collected ~1wk after birth)


dist.pups.to.dams.wuni.metadata %>% 
  ggplot(aes(x = PostnatalDay.x, y = dist, group = AnimalID.x)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Group.x + Dam.x)



### COMPARE BASELINE BACTERIOME FOR TWO LIVE-FVT LITTERS ###
ord.pcoa.uni <- ordinate(ps.relabund.pups %>% 
                           subset_samples(PrePostFVT %in% c("Pre")) %>% 
                           subset_samples(Group %in% c(1)), method="PCoA", distance="wunifrac")
plot_ordination(ps.relabund.pups, ord.pcoa.uni, color="Dam", title="") +
  stat_ellipse(geom="polygon",type = "t", alpha = 0.1, aes(fill=Dam)) +
  theme_bw() +
  #scale_color_manual(values = c("#00A087FF", "#E64B35FF"), labels = c("Live FVT", "hkFVT")) +
  #scale_fill_manual(values = c("#00A087FF", "#E64B35FF"), labels = c("Live FVT", "hkFVT")) +
  ##scale_color_manual(values = c("#00A087FF", "#3C5488FF"), labels = c("Live FVT", "H control")) +
  ##scale_fill_manual(values = c("#00A087FF", "#3C5488FF"), labels = c("Live FVT", "H control")) +
  #scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "hkFVT", "Healthy\ncontrol")) +
  #scale_fill_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF")) +
  #  scale_shape_manual(values = c(16, 8, 1), labels = c("Pre-FVT", "During FVT", "Post-FVT")) +
  facet_wrap(~PrePostFVT) +
  theme(aspect.ratio = 1) +
  labs(col = "Treatment") +
  guides(fill = "none", shape = "none")



