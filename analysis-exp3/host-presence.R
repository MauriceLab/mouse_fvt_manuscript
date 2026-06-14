library(tidyverse)
library(phyloseq)

### IPHOP IMPORT
iphop_out <- read_csv("Host_prediction_to_genome_m90.csv")

iphop_out.dedup <- iphop_out %>% group_by(Virus) %>% arrange(desc(`Confidence score`)) %>% slice(1)

iphop_out.dedup <- iphop_out.dedup[, c("Virus", "Host taxonomy", "Host genome")]

iphop_out.split <- iphop_out.dedup %>% 
  separate(col = `Host taxonomy`, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

iphop_out.split$domain <- gsub("d__", "", iphop_out.split$domain)
iphop_out.split$phylum <- gsub("p__", "", iphop_out.split$phylum)
iphop_out.split$class <- gsub("c__", "", iphop_out.split$class)
iphop_out.split$order <- gsub("o__", "", iphop_out.split$order)
iphop_out.split$family <- gsub("f__", "", iphop_out.split$family)
iphop_out.split$genus <- gsub("g__", "", iphop_out.split$genus)
iphop_out.split$species <- gsub("s__", "", iphop_out.split$species)

iphop_out.split$phylum <- sapply(strsplit(iphop_out.split$phylum, "_"), `[`, 1)
iphop_out.split$class <- sapply(strsplit(iphop_out.split$class, "_"), `[`, 1)
iphop_out.split$order <- sapply(strsplit(iphop_out.split$order, "_"), `[`, 1)
iphop_out.split$genus <- sapply(strsplit(iphop_out.split$genus, "_"), `[`, 1)
iphop_out.split$species <- sapply(strsplit(iphop_out.split$species, "_"), `[`, 1)

fvt.host.targets.species <- ps.phage.relabund %>% 
  subset_samples(SampleType == "Inoculum") %>% 
  subset_samples(SampleName == "VIVO8_VLP_DONOR2") %>% 
  psmelt() %>% 
  select(OTU, Abundance) %>% 
  filter(Abundance > 0) %>% 
  left_join(iphop_out.split %>% select(Virus, family, genus, species), by = c("OTU" = "Virus")) %>% 
  filter(!is.na(species)) %>% 
  group_by(species, genus) %>% 
  summarise(sum.relabund = sum(Abundance))

fvt.host.targets.genus <- fvt.host.targets.species %>% 
  group_by(genus) %>% 
  summarise(sum.relabund = sum(sum.relabund))

# stats on percent vOTUs assigned host, their cumulative relative abundance
ps.phage.relabund %>% 
  subset_samples(SampleType == "Inoculum") %>% 
  subset_samples(SampleName == "VIVO8_VLP_DONOR2") %>% 
  psmelt() %>% 
  select(OTU, Abundance) %>% 
  filter(Abundance > 0) %>% 
  left_join(iphop_out.split %>% select(Virus, family, genus, species), by = c("OTU" = "Virus")) %>% 
  select(species, Abundance) %>% 
  ungroup() %>% 
  summarize(
    percent_assigned = sum(!is.na(species)) / n(),
    n_assigned = sum(!is.na(species)),
    n_na = sum(is.na(species)),
    n_total = n(),
    sum.relabund.assigned = sum(Abundance[!is.na(species)])
  )

bact.hosts.present <- ps.relabund %>% 
  subset_samples(SampleType == "Pup") %>% 
  subset_samples(PrePostFVT == "Pre") %>% 
  psmelt() %>% 
  select(SampleName, Abundance, AnimalID, Dam, Group, PostnatalDay, genus, species) %>% 
  filter(Abundance > 0) %>% 
  group_by(SampleName, AnimalID, Group, Dam, PostnatalDay) %>% 
  summarize(species_list = list(species), genus_list = list(unique(genus)))


host.presence.results <- bact.hosts.present %>% 
  rowwise() %>% 
  mutate(
    species_targets_present = 
      list(intersect(species_list, fvt.host.targets.species$species)),
    prop_species_targets_present = 
      length(intersect(species_list, fvt.host.targets.species$species)) / 
      length(fvt.host.targets.species$species),
    genus_targets_present = 
      list(intersect(genus_list, fvt.host.targets.genus$genus)),
    prop_genus_targets_present = 
      length(intersect(genus_list, fvt.host.targets.genus$genus)) / 
      length(fvt.host.targets.genus$genus)
    
  ) %>%
  ungroup()

host.presence.results %>% 
  filter(Group == 1) %>% 
  group_by(AnimalID) %>% 
  arrange(desc(PostnatalDay)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  summarize(mean_species = mean(prop_species_targets_present),
            mean_genus = mean(prop_genus_targets_present))


host.presence.results %>% 
  filter(Group == 1) %>% 
  group_by(AnimalID) %>% 
  arrange(desc(PostnatalDay)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  pull(genus_targets_present)
# Bacteroides, Bifidobacterium