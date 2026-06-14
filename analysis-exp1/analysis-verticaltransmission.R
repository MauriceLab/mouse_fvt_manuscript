library(tidyverse)
library(reshape2)
library(phyloseq)


#### BACTERIOME OVERLAP WITH DAMS ####
ps.vertical.relabund <- ps.relabund %>% 
  subset_samples((SampleType == "Breeder" & SampleShortName %in% c("F11-PRENATAL", "F12-PRENATAL", "F14-PRENATAL", "F15-PRENATAL", "F17-PRENATAL", "F2-PRENATAL", "F3-PRENATAL", "F8-PRENATAL", "F9-PRENATAL")) | 
                   SampleType == "Pup")
metadata.vertical <- sample_data(ps.relabund)

otutable.f2.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F2-PRENATAL") | (Dam == "VIVO6-F2"))) %>% as.data.frame()
otutable.f3.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F3-PRENATAL") | (Dam == "VIVO6-F3"))) %>% as.data.frame()
otutable.f8.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F8-PRENATAL") | (Dam == "VIVO6-F8"))) %>% as.data.frame()
otutable.f9.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F9-PRENATAL") | (Dam == "VIVO6-F9"))) %>% as.data.frame()
otutable.f11.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F11-PRENATAL") | (Dam == "VIVO6-F11"))) %>% as.data.frame()
otutable.f12.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F12-PRENATAL") | (Dam == "VIVO6-F12"))) %>% as.data.frame()
otutable.f14.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F14-PRENATAL") | (Dam == "VIVO6-F14"))) %>% as.data.frame()
otutable.f15.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F15-PRENATAL") | (Dam == "VIVO6-F15"))) %>% as.data.frame()
otutable.f17.pups <- otu_table(ps.vertical.relabund %>% subset_samples((SampleShortName == "F17-PRENATAL") | (Dam == "VIVO6-F17"))) %>% as.data.frame()

overlap.percent.f2.pups <- otutable.f2.pups %>% 
  filter(`VIVO6-FEC-16S-MAY5_F2` > 0) %>%
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f3.pups <- otutable.f3.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY5_F3` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f8.pups <- otutable.f8.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY5_F8` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f9.pups <- otutable.f9.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY3_F9` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f11.pups <- otutable.f11.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY5_F11` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f12.pups <- otutable.f12.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY5_F12` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f14.pups <- otutable.f14.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY5_F14` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f15.pups <- otutable.f15.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY12_F15` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")

overlap.percent.f17.pups <- otutable.f17.pups %>% ### CHANGE
  filter(`VIVO6-FEC-16S-MAY3_F17` > 0) %>% ### CHANGE
  pivot_longer(cols = everything(), names_to = 'sample', values_to = 'abundance') %>% 
  filter(abundance > 0) %>% 
  left_join(metadata.vertical, by = c("sample" = "SampleName")) %>% 
  group_by(sample) %>% 
  summarize(relabund.common = sum(abundance), relabund.other = 1-sum(abundance)) %>% 
  left_join(metadata.vertical[, c("SampleName", "DonorID", "PostnatalDay", "SampleType", "Dam", "Group")], by = c("sample" = "SampleName")) %>% 
  filter(SampleType != "Breeder")


overlap.percent.all.pups <- rbind(overlap.percent.f2.pups, overlap.percent.f3.pups, overlap.percent.f8.pups, overlap.percent.f9.pups, overlap.percent.f11.pups, overlap.percent.f12.pups, overlap.percent.f14.pups, overlap.percent.f15.pups, overlap.percent.f17.pups) %>% 
  pivot_longer(cols = c(relabund.common, "relabund.other"), names_to = "names", values_to = "values") %>% 
  group_by(DonorID, names, PostnatalDay, Dam, Group) %>% 
  summarise(abundance = mean(values), sd = sd(values))
overlap.percent.all.pups[overlap.percent.all.pups == "relabund.common"] <- "Shared with dam"
overlap.percent.all.pups[overlap.percent.all.pups == "relabund.other"] <- "Pup only"
overlap.percent.all.pups$PostnatalDay <- as.factor(overlap.percent.all.pups$PostnatalDay)
overlap.percent.all.pups["sd"][overlap.percent.all.pups["names"] == "Pup only"] <- NA

# fig 1d
fig.1c <- rbind(overlap.percent.f2.pups, overlap.percent.f3.pups, overlap.percent.f8.pups, overlap.percent.f9.pups, overlap.percent.f11.pups, overlap.percent.f12.pups, overlap.percent.f14.pups, overlap.percent.f15.pups, overlap.percent.f17.pups) %>% 
  pivot_longer(cols = c(relabund.common, "relabund.other"), names_to = "names", values_to = "values") %>% 
  group_by(DonorID, names, PostnatalDay, Dam, Group) %>% 
  filter(names == "relabund.common") %>% 
  ungroup() %>% 
  left_join(metadata %>% select(SampleName, AnimalID, SampleShortName), 
            by = c("sample" = "SampleName")) %>% 
  mutate(Group = as.factor(Group)) %>% 
  group_by(sample) %>% 
  mutate(
    x_jitter = PostnatalDay + runif(n(), -0.15, 0.15),
    y_jitter = values + runif(n(), -0.015, 0.015)
  ) %>% 
  ungroup() %>% 
  filter(!grepl("POOL", sample)) %>% 
  ggplot(aes(x = x_jitter, y = y_jitter, color = Group, group = AnimalID)) +
    geom_point(alpha = 1) +
    geom_line(alpha = 0.7) +
    #scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), name = "Group", labels = c("H donor + CON diet",  "S donor + MAL diet", "H donor + MAL diet")) +
    scale_colour_manual(values = c("#c57bea", "#779dee", "#d0ae7c"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(labels = seq(from = 21, to = 49, by = 7), breaks = seq(from = 21, to = 49, by = 7), limits = c(21, 45)) +
    xlab("Age of pup (days)") +
    ylab("Rel. abundance of\nvertically-transmitted taxa") +
    theme(aspect.ratio = 1, legend.position = "none",
          axis.title.y = element_text(size = 10))

fig.1c
