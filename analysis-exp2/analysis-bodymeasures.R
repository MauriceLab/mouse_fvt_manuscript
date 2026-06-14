library(tidyverse)
library(reshape2)

# import pups measures from exp 1, 2, 3
vivo6.pups.measures.g1 <- readRDS("analysis-exp2/vivo6-pups-measures-g1.rds")

vivo8.pups.measures <- readRDS("analysis-exp2/vivo8-pups-measures.rds") %>% 
  mutate(experiment = "VIVO8")

pups.measures <- read_csv("analysis-exp2/vivo7-pups-measures.csv") %>% 
  mutate(experiment = "VIVO7")
pups.measures$sex <- as.factor(pups.measures$sex)
pups.measures$animal <- as.factor(pups.measures$animal)
pups.measures$dam <- as.factor(pups.measures$dam)
pups.measures$cage <- as.factor(pups.measures$cage)
pups.measures$group <- as.factor(pups.measures$group)

pups.measures.merged <- pups.measures %>% 
  plyr::rbind.fill(vivo8.pups.measures) %>% 
  filter(group %in% c(1,2,3))

# modelling of pup body measures
library(mgcv)
library(emmeans)

model.weight.gam <- gam(weight ~ s(postnatal_day, by=group) + group + s(dam, bs = "re") + s(dam, postnatal_day, bs="re"), data = pups.measures, method = "REML")
summary(model.weight.gam)
plot(model.weight.gam, pages=1, se=TRUE)
anova(model.weight.gam)
emmeans(model.weight.gam, pairwise ~ group | postnatal_day)


model.tail.gam <- gam(length_tail ~ s(postnatal_day, by=group) + group + s(dam, bs = "re") + s(dam, postnatal_day, bs="re"), data = pups.measures, method = "REML")
summary(model.tail.gam)
plot(model.tail.gam, pages=1, se=TRUE)
anova(model.tail.gam)
emmeans(model.tail.gam, pairwise ~ group | postnatal_day)


model.length.gam <- gam(length_body ~ s(postnatal_day, by=group) + group + s(dam, bs = "re") + s(dam, postnatal_day, bs="re"), data = pups.measures, method = "REML")
summary(model.length.gam)
plot(model.length.gam, pages=1, se=TRUE)
anova(model.length.gam)
emmeans(model.length.gam, pairwise ~ group | postnatal_day)


# figure 3b
vivo6.pups.measures.g1.long <- 
  vivo6.pups.measures.g1 %>% 
  select(postnatal.day, weight, length.tail, length.body) %>% 
  rename(length_tail = length.tail, length_body = length.body) %>%
  pivot_longer(c(weight, length_tail, length_body), names_to = "metric", values_to = "measure") %>% 
  filter(!is.na(measure))
 
colnames(vivo6.pups.measures.g1.long) <- c("day", "metric", "measure")

fig.3b <- pups.measures.merged %>% 
  select(day, weight, length_body, experiment, group) %>% 
  pivot_longer(c(weight, length_body), names_to = "metric", values_to = "measure") %>% 
  filter(!is.na(measure)) %>% 
  ggplot(aes(x = day, y = measure, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(data = vivo6.pups.measures.g1.long %>% filter(metric %in% c("length_body", "weight")), method = "lm", 
              formula = y ~ splines::ns(x, 4), 
              aes(x = day , y = measure), 
              linetype = "dashed", se = F, color = "black", linewidth = 0.75) +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 4), linewidth = 0.75) +
  xlab("Age (days)") +
  ylab("Body measure") +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), name = "Group", labels = c("Live FVT", "TI-FVT", "rPBS")) +
  theme_bw() +
  facet_grid(metric~experiment, scale = "free_y", 
             labeller = as_labeller(c(length_body = "Body length (mm)",
                                      weight = "Body weight (g)",
                                      VIVO7 = "Experiment #2",
                                      VIVO8 = "Experiment #3") )) + 
theme(legend.position = "bottom") +
guides(col = guide_legend(nrow = 1))

fig.3b

#### SUPPLEMENTARY FIGURE 1 ####

fig.s1b <- pups.measures.merged %>% 
  select(day, length_tail, experiment, group) %>% 
  pivot_longer(c(length_tail), names_to = "metric", values_to = "measure") %>% 
  filter(!is.na(measure)) %>% 
  ggplot(aes(x = day, y = measure, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(data = vivo6.pups.measures.g1.long %>% 
                filter(metric == "length_tail"), method = "lm", 
              formula = y ~ splines::ns(x, 4), 
              aes(x = day , y = measure), 
              linetype = "dashed", se = F, color = "black", linewidth = 0.75) +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 4), linewidth = 0.75) +
  xlab("Age (days)") +
  ylab("Tail length (mm)") +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), name = "Group", labels = c("Live FVT", "TI-FVT", "rPBS")) +
  theme_bw() +
  facet_wrap(~experiment, scale = "free_y", 
             labeller = as_labeller(c(VIVO7 = "Experiment #2",
                                      VIVO8 = "Experiment #3") )) + 
  theme(aspect.ratio=1, legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2))

fig.s1b
