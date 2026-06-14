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


emu_tax <- read_tsv("analysis-exp3/emu-combined-taxonomy-tax_id.tsv") %>% filter(tax_id != "unassigned")
emu_tax.mat <- as.matrix(emu_tax[, -1])
rownames(emu_tax.mat) <- emu_tax$tax_id
emu_tax.mat <- emu_tax.mat[, rev(seq_len(ncol(emu_tax.mat)))]

emu_counts <- read_tsv("analysis-exp3/emu-combined-abundance-tax_id-counts.tsv") %>% filter(tax_id != "unassigned")
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
  "Dam" = c("F14", "F12", "F2"),
  "FVT.date" = c("FEB19", "FEB19", "FEB11")
) %>% 
  mutate(FVT.date = lubridate::mdy(paste0(FVT.date, "/2025")))

sample.depth <- colSums(emu_counts.mat, na.rm = T) %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName")
colnames(sample.depth) <- c("SampleName", "ReadDepth")

metadata <- read.csv("analysis-exp3/16s-metadata.csv")
metadata$Group <- factor(metadata$Group)
metadata$DonorID <- factor(metadata$DonorID)
metadata$AnimalID <- factor(metadata$AnimalID)
metadata$SampleType <- factor(metadata$SampleType)
metadata$Dam <- factor(metadata$Dam)
metadata$PrePostFVT <- factor(metadata$PrePostFVT, levels = c("Pre", "During", "Post", "Post-Final"))

metadata <- metadata %>% 
  mutate(Date.formatted = if_else(grepl("DEC", Date, fixed = T), 
                                  lubridate::mdy(paste0(Date, "/2024")), 
                                  lubridate::mdy(paste0(Date, "/2025")))) %>% 
  left_join(fvt.dates, by = c("Dam" = "Dam")) %>% 
  mutate(RelativeDay = as.numeric(Date.formatted - FVT.date)) %>% 
  left_join(sample.depth)
  

rownames(metadata) <- metadata$SampleName

METADATA <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, TREE, METADATA)
ps.relabund <- transform_sample_counts(ps, function(x) x/sum(x))



ps.relabund.pups <- ps.relabund %>% subset_samples(SampleType == "Pup")
ps.relabund.dams <- ps.relabund %>% subset_samples(SampleType == "Breeder")




# beta diversity
ord.pcoa.uni <- ordinate(ps.relabund.pups %>% 
                           subset_samples(PrePostFVT %in% c("Pre", "Post")), method="PCoA", distance="wunifrac")
fig.s2e <- plot_ordination(ps.relabund.pups, ord.pcoa.uni, color="Group", title="") +
  stat_ellipse(geom="polygon",type = "t", alpha = 0.1, aes(fill=Group)) +
  theme_bw() +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "TI-FVT", "PBS")) +
  scale_fill_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), labels = c("Live FVT", "TI-FVT", "PBS")) +
  facet_wrap(~PrePostFVT) +
  theme(aspect.ratio = 1, legend.position = "bottom") +
  labs(col = "Group") +
  guides(fill = "none", shape = "none") +
  xlab("Axis.1 [46.8%]") +
  ylab("Axis.2 [21.0%]") +
  ggtitle("Bacteriome beta diversity")


fig.s2e

#####


# PERMANOVA: did pups shift from baseline?
dist.wuni <- distance(ps.relabund.pups %>% subset_samples(Group %in% c(1, 2, 3)) %>% subset_samples(Group == 1), method = "wunifrac")
dist.wuni.metadata <- metadata %>% filter(SampleName %in% (dist.wuni %>% as.matrix() %>% colnames()))
adonis2(dist.wuni ~ PrePostFVT, 
        data = dist.wuni.metadata, 
        strata = dist.wuni.metadata$AnimalID,
        permutations = 999,
        by = "terms")
# not significant, PrePostFVT p = 0.646, R2 = 0.088

# PERMANOVA: are live FVT pups distinct from HK-FVT/rPBS after treatment?
dist.wuni.final <- distance(ps.relabund.pups %>% subset_samples(Group %in% c(1, 2, 3)) %>% subset_samples(PrePostFVT == "Post"), method = "wunifrac")
dist.wuni.final.metadata <- metadata %>% filter(SampleName %in% (dist.wuni.final %>% as.matrix() %>% colnames()))
adonis2(dist.wuni.final ~ Group + PostnatalDay, 
        data = dist.wuni.final.metadata, 
        permutations = 999,
        by = "terms")
# not significant, Group p = 0.208, PostnatalDay p = 0.247


richness <- estimate_richness(ps, measures=c("Shannon", "Observed", "Chao1")) %>% 
  rownames_to_column("sample") %>% 
  mutate(sample = gsub("\\.", '-', sample)) %>% 
  left_join(metadata, by = c("sample" = "SampleName")) #%>% 

# modelling of pup richness
library(mgcv)
library(emmeans)

richness.data <- 
  richness %>% 
  filter(SampleType == "Pup") %>% 
  mutate(AnimalID = as.factor(AnimalID))

model.richness.gam <- gam(Observed ~ Group + log10(ReadDepth) + s(RelativeDay, by=Group) + s(AnimalID, bs = "re"), 
                          data = richness.data, family = gaussian(), method = "REML")

model.shannon.gam <- gam(Shannon ~ Group + log10(ReadDepth) + s(RelativeDay, by=Group) + s(AnimalID, bs = "re"), 
                          data = richness.data, family = gaussian(), method = "REML")


summary(model.richness.gam)
draw(model.richness.gam)

library(gratia)
fig.s2d1 <- model.richness.gam %>% 
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

fig.s2d1

# distance from previous timepoint (stability)
dist.pups.wuni <-  as.matrix(distance(ps.relabund.pups, method = "wunifrac")) %>% 
  melt() %>% 
  filter(Var1 != Var2)

colnames(dist.pups.wuni) <- c("sample1", "sample2", "dist")

dist.pups.wuni.metadata <- dist.pups.wuni %>% 
  left_join(metadata, by = c("sample1" = "SampleName")) %>% 
  left_join(metadata, by = c("sample2" = "SampleName"))

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
  facet_wrap(~Group.x, scales = "free") +
  #scale_y_continuous(limits = c(0, 1)) +
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
  #ggplot(aes(x = PostnatalDay.y, y = dist, group = AnimalID.x, color = Dam.x)) +
  ggplot(aes(x = RelativeDay.y, y = dist, group = AnimalID.x)) +
  geom_line(color = "black") +
  geom_point() +
  #geom_smooth(se = F, method="lm",
  #            formula=  y ~ splines::bs(x, degree = 2)) + 
  facet_wrap(~Group.x, labeller = as_labeller(c("1" = "Live FVT", "2" = "hkFVT", "3" = "PBS"))) +
  #scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  #scale_color_viridis_d() +
  labs(x = "Postnatal Day", y = "Distance to baseline") +
  theme(aspect.ratio = 1)

# LAST TIMEPOINT ONLY
dist.pups.wuni.metadata.samepup.baseline %>% 
  slice_max(PostnatalDay.y, n = 1, with_ties = FALSE) %>%
  filter(PostnatalDay.y > 70) %>% 
  mutate(Group.x = fct_recode(Group.x, "Live FVT" = "1", "hkFVT" = "2", "rPBS" = "3")) %>% 
  ggplot(aes(x = Group.x, y = dist)) +
  geom_boxplot() +
  geom_jitter(aes(color = Dam.x), position = position_dodge2(width = 0.5)) +
  theme_bw() +
  labs(x = "Group", y = "Distance to baseline", subtitle = "VIVO8 - Bacteriome") +
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
  filter(Date.y %in% c("JAN17", "JAN23", "JAN28")) #dam postpartum timepoint (collected ~1wk after birth)


dist.pups.to.dams.wuni.metadata %>% 
  ggplot(aes(x = PostnatalDay.x, y = dist, group = AnimalID.x)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Group.x)



### CHECK GROUP 1 PUPS FOR FVT-TARGETING TAXA AT BASELINE ###


votus.pups.pre <- ps.relabund %>% 
  subset_samples(Group == 1) %>% 
  subset_samples(PrePostFVT == "Pre") %>% 
  #sample_data() %>% view()
  psmelt() %>% 
  filter(Abundance > 0) %>% 
  pull(OTU) %>% 
  unique()

emu_tax %>% filter(tax_id %in% votus.pups.pre) %>% pull(genus) %>% unique() %>% 
  dput()

