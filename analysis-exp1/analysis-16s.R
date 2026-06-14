library(tidyverse)
library(reshape2)
library(phyloseq)
library(vegan)
library(mgcv)
library(emmeans)
library(gratia)

#### BEGIN 16S ANALYSIS ####
emu_tax <- read_tsv("analysis-exp1/emu-combined-taxonomy-tax_id.tsv") %>% filter(tax_id != "unassigned")
emu_tax.mat <- as.matrix(emu_tax[, -1])
rownames(emu_tax.mat) <- emu_tax$tax_id

emu_counts <- read_tsv("analysis-exp1/emu-combined-abundance-tax_id-counts.tsv") %>% filter(tax_id != "unassigned")
emu_counts.mat <- as.matrix(emu_counts[, -1])
rownames(emu_counts.mat) <- emu_counts$tax_id
colnames(emu_counts.mat) <- gsub(".fq", "", colnames(emu_counts.mat))
emu_counts.mat <- floor(emu_counts.mat)

depth <- colSums(emu_counts.mat, na.rm = T) %>% 
  melt() %>% 
  rownames_to_column("SampleName")
colnames(depth) <- c("SampleName", "Depth")

#### PHYLOSEQ SETUP ####
TAX <- tax_table(emu_tax.mat)
emu_counts.mat.zeros <- emu_counts.mat
emu_counts.mat.zeros[is.na(emu_counts.mat.zeros)] <- 0
OTU <- otu_table(emu_counts.mat.zeros, taxa_are_rows = T,)
TREE <- read_tree("emu_tax_tree_species.nwk")


metadata <- read.csv("analysis-exp1/16s-metadata.csv")
metadata$Batch <- factor(metadata$Batch)
metadata$Group <- factor(metadata$Group)
metadata$DonorID <- factor(metadata$DonorID)
metadata$SampleShortName <-  factor(metadata$SampleShortName, 
                                    levels = c("D1-INOCULUM", "D2-INOCULUM", "D3-INOCULUM", "D4-INOCULUM", "F11-GF", "F12-GF", "F14-GF", "F15-GF", "F17-GF", "F2-GF", "F3-GF", "F8-GF", "F9-GF", "F11-FMT+14", "F12-FMT+14", "F14-FMT+14", "F15-FMT+14", "F17-FMT+14", "F2-FMT+14", "F3-FMT+14", "F8-FMT+14", "F9-FMT+14", "F11-PRENATAL", "F12-PRENATAL", "F14-PRENATAL", "F15-PRENATAL", "F17-PRENATAL", "F2-PRENATAL", "F3-PRENATAL", "F8-PRENATAL", "F9-PRENATAL", "F11-POSTNATAL", "F12-POSTNATAL", "F14-POSTNATAL", "F15-POSTNATAL", "F17-POSTNATAL", "F2-POSTNATAL", "F3-POSTNATAL", "F8-POSTNATAL", "F9-POSTNATAL", "F11-WEAN", "F12-WEAN", "F14-WEAN", "F15-WEAN", "F17-WEAN", "F2-WEAN", "F3-WEAN", "F8-WEAN", "F9-WEAN", "POOL-F15-P20", "POOL-F8-P20", "POOL-F11-P21", "POOL-F14-P21", "POOL-F9-P22", "F23-P21", "F24-P21", "F25-P21", "F26-P21", "F27-P21", "F28-P21", "M19-P21", "M20-P21", "M21-P21", "M22-P21", "F45-P22", "F46-P22", "F48-P22", "F49-P22", "F50-P22", "F30-P23", "F31-P23", "M29-P23", "F26-P24", "F27-P24", "F28-P24", "F41-P26", "M39-P26", "M42-P26", "M43-P26", "M44-P26", "F23-P27", "F24-P27", "F25-P27", "F26-P27", "F28-P27", "F33-P27", "F34-P27", "F38-P27", "M22-P27", "M32-P27", "M36-P27", "M37-P27", "M40-P28", "F48-P30", "F49-P30", "F23-P32", "F24-P32", "F25-P32", "M19-P32", "M20-P32", "M21-P32", "M22-P32", "F30-P33", "F31-P33", "F38-P33", "M29-P33", "M36-P33", "M37-P33", "M42-P33", "M43-P33", "M44-P33", "F45-P34", "F46-P34", "F38-P35", "F41-P35", "F48-P35", "F49-P35", "F50-P35", "M36-P35", "M37-P35", "M39-P35", "M40-P35", "M42-P35", "M43-P35", "M44-P35", "F33-P37", "F34-P37", "M32-P37", "F23-P38", "F24-P38", "F25-P38", "F26-P38", "F27-P38", "F28-P38", "F41-P38", "M19-P38", "M20-P38", "M21-P38", "M22-P38", "M39-P38", "M40-P38", "F30-P39", "F31-P39", "F38-P39", "M29-P39", "M36-P39", "M37-P39", "M42-P40", "M43-P40", "M44-P40", "F45-P41", "F46-P41", "F48-P42", "F49-P42", "F33-P43", "F34-P43", "M32-P43"))
                                    #ordered by timepoint

metadata <- metadata %>% left_join(depth)
rownames(metadata) <- metadata$SampleName

METADATA <- sample_data(metadata)

ps <- phyloseq(OTU, TAX, TREE, METADATA)
ps.relabund <- transform_sample_counts(ps, function(x) x/sum(x))


ps.relabund.pups <- ps.relabund %>% subset_samples(SampleType == "Pup")
metadata.pups <- data.frame(sample_data(ps.relabund.pups))


#### ALPHA DIVERSITY ####
ps.richness <- ps %>% 
  subset_samples(SampleType == "Pup") %>% 
  estimate_richness(measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "SampleName") %>%
  mutate(SampleName = str_replace_all(SampleName, "\\.", "-")) %>% 
  left_join(metadata, by = c("SampleName" = "SampleName"))

ps.richness %>% 
  melt() %>% 
  ggplot(aes(x = PostnatalDay, y = value, color = Group, group = interaction(Group, PostnatalDay))) +
    #geom_point() +
    geom_smooth(method = "lm", se = F, formula = y ~ splines::bs(x, 3), aes(group = Group)) +
    geom_boxplot() +
    facet_wrap(. ~ variable, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    theme_bw()


ps.richness$AnimalID <- as.factor(ps.richness$AnimalID)
ps.richness$Group <- as.factor(ps.richness$Group)
ps.richness$DonorID <- as.factor(ps.richness$DonorID)
ps.richness$logDepth <- log10(ps.richness$Depth)
model.richness.gam <- gam(Observed ~ s(PostnatalDay) + DonorID + Group + logDepth + s(PostnatalDay, by = Group) + s(AnimalID, bs = "re"),
                          data = ps.richness,
                          family = nb(),
                          method = "REML")
summary(model.richness.gam)
richness.emms <- emmeans(model.richness.gam, pairwise ~ Group | PostnatalDay, 
                         at = list(PostnatalDay = seq(20, 45, by = 2)),
                         type = "response")
as.data.frame(richness.emms$contrasts) %>%
  filter(p.value < 0.05) %>% 
  view()

model.shannon.gam <- gam(Shannon ~ s(PostnatalDay) + DonorID + Group + logDepth + s(PostnatalDay, by = Group) + s(AnimalID, bs = "re"),
                         data = ps.richness,
                         family = gaussian(),
                         method = "REML")
summary(model.shannon.gam)
shannon.emms <- emmeans(model.shannon.gam, pairwise ~ Group | PostnatalDay, 
                         at = list(PostnatalDay = seq(20, 45, by = 2)),
                         type = "response")
as.data.frame(shannon.emms$contrasts) %>%
  filter(p.value < 0.05) %>% 
  view()


ds.richness <- data_slice(model.richness.gam, 
                 PostnatalDay = seq(20, 45, length.out = 100), 
                 Group = unique(ps.richness$Group))
fv.richness <- fitted_values(model.richness.gam, data = ds.richness, scale = "response") %>% mutate(metric = "richness")


ds.shannon <- data_slice(model.shannon.gam, 
                          PostnatalDay = seq(20, 45, length.out = 100), 
                          Group = unique(ps.richness$Group))
fv.shannon <- fitted_values(model.shannon.gam, data = ds.richness, scale = "response") %>% mutate(metric = "shannon")

fig.2a <- 
  fv.richness %>% rbind(fv.shannon) %>% 
  ggplot(aes(x = PostnatalDay, y = .fitted, color = Group)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = Group), linewidth = 0.1, alpha = 0.07) +
  scale_colour_manual(values = c("#DCB0F2FF", "#9EB9F3FF", "#D3B484FF"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
  scale_fill_manual(values = c("#DCB0F2FF", "#9EB9F3FF", "#D3B484FF"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
  theme_bw() + 
  labs(title = "Alpha diversity",
       x = "Age (days)",
       y = "Fitted value\n(total effect)") +
  theme(legend.key.spacing.y = unit(0.15, 'cm'),
        legend.position = "bottom"
  ) +
  facet_wrap(~ metric, scale = "free_y", labeller = as_labeller(c(richness = "Richness", shannon = "Shannon"))) +
  scale_x_continuous(breaks = seq(0, 42, by = 6)) +
  guides(colour = guide_legend(nrow = 2))

fig.2a


#### BETA DIVERSITY PLOTS ####
ord.pcoa.wuni <- ordinate(ps.relabund.pups, 
                         method="PCoA", distance="wunifrac")

ord.pcoa.wuni.axis1.variance <- round(ord.pcoa.wuni$values$Relative_eig[1]*100, 1)
ord.pcoa.wuni.axis2.variance <- round(ord.pcoa.wuni$values$Relative_eig[2]*100, 1)

pcoa.wuni.age <- 
  ord.pcoa.wuni$vectors %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName") %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = PostnatalDay)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "lightblue", high = "navyblue") + 
  xlab(paste("Axis.1 [", ord.pcoa.wuni.axis1.variance, "%]", sep = "")) +
  ylab(paste("Axis.2 [", ord.pcoa.wuni.axis2.variance, "%]", sep = "")) +
  ggtitle("Age (days)") +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank())

pcoa.wuni.age

pcoa.wuni.group <- 
  ord.pcoa.wuni$vectors %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName") %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point() +
  theme_bw() +
  #paletteer::scale_colour_paletteer_d("ggthemes::Classic_Green_Orange_12", labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
  scale_colour_manual(values = c("#c57bea", "#779dee", "#d0ae7c"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
  scale_fill_manual(values = c("#c57bea", "#779dee", "#d0ae7c"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
  xlab(paste("Axis.1 [", ord.pcoa.wuni.axis1.variance, "%]", sep = "")) +
  ylab(paste("Axis.2 [", ord.pcoa.wuni.axis2.variance, "%]", sep = "")) +
  ggtitle("Group") +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank()) +
  guides(colour = guide_legend(nrow = 2))

pcoa.wuni.group

pcoa.wuni.donor <- 
  ord.pcoa.wuni$vectors %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleName") %>% 
  left_join(metadata) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = DonorID)) +
  geom_point() +
  theme_bw() +
  paletteer::scale_colour_paletteer_d("rcartocolor::Pastel", labels = c("Donor 1 (H)", "Donor 2 (H)", "Donor 3 (S)", "Donor 4 (S)")) +
  xlab(paste("Axis.1 [", ord.pcoa.wuni.axis1.variance, "%]", sep = "")) +
  ylab(paste("Axis.2 [", ord.pcoa.wuni.axis2.variance, "%]", sep = "")) +
  ggtitle("Donor") +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank()) +
  guides(colour = guide_legend(nrow = 2))

pcoa.wuni.donor

#### PERMANOVA ####
dist.wuni <- distance(ps.relabund.pups, method = "wunifrac")
perm <- how(plots = Plots(strata = metadata.pups$AnimalID), nperm = 999)
set.seed(123)
permanova.dist.wuni <- adonis2(dist.wuni ~ log10(Depth) + PostnatalDay + DonorID + Group, 
                               data = metadata.pups,
                               permutations = perm,
                               by = "terms")
permanova.dist.wuni



