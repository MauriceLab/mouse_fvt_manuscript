library(tidyverse)
library(reshape2)

# pups weights

vivo6.pups.measures.g1 <- readRDS("~/Desktop/VIVO7/analysis/vivo6-pups-measures-g1.rds")

pups.measures <- read_csv("vivo8-pups-measures-R.csv")

pups.measures$sex <- as.factor(pups.measures$sex)
pups.measures$animal <- as.factor(pups.measures$animal)
pups.measures$dam <- as.factor(pups.measures$dam)
pups.measures$cage <- as.factor(pups.measures$cage)
pups.measures$group <- as.factor(pups.measures$group)


saveRDS(pups.measures, "vivo8-pups-measures.rds")

graph.pups.weight <- pups.measures %>% 
  ggplot(aes(x = day, y = weight, group = group, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(data = vivo6.pups.measures.g1, method = "lm", formula = y ~ splines::ns(x, 4), aes(x = postnatal.day , y = as.numeric(weight)), linetype = "dashed", se = F, color = "black") +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 4)) +
  xlab("Postnatal day") +
  ylab("Body weight (g)") +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF"), name = "Group", labels = c("S donor + FVT (MAL diet)", "S donor + HK-FVT (MAL diet)", "S donor + rPBS (MAL diet)")) +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = "none")#, axis.title.x=element_blank())

graph.pups.weight

graph.pups.tail <- pups.measures %>% 
  ggplot(aes(x = day, y = length_tail, group = group, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(data = vivo6.pups.measures.g1, method = "lm", formula = y ~ splines::ns(x, 4), aes(x = postnatal.day , y = as.numeric(length.tail)), linetype = "dashed", se = F, color = "black") +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 2)) +
  xlab("Postnatal day") +
  ylab("Tail length (mm)") +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF")) +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = "none")#, axis.title.x=element_blank())

graph.pups.tail

graph.pups.body <- pups.measures %>% 
  ggplot(aes(x = day, y = length_body, group = group, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 2)) +
  geom_smooth(data = vivo6.pups.measures.g1, method = "lm", formula = y ~ splines::ns(x, 4), aes(x = postnatal.day , y = as.numeric(length.body)), linetype = "dashed", se = F, color = "black") +
  xlab("Postnatal day") +
  ylab("Body length (mm)") +
  scale_color_manual(values = c("#00A087FF", "#E64B35FF", "#3C5488FF")) +
  theme_bw() +
  theme(aspect.ratio=1, legend.position = "none")#, axis.title.x=element_blank())

graph.pups.body

ggpubr::ggarrange(graph.pups.weight, graph.pups.tail, graph.pups.body, nrow = 1, ncol = 3, common.legend = TRUE, legend = "right")



# modelling of pup body measures
library(mgcv)
library(emmeans)

model.weight.gam <- gam(weight ~ s(day, by=group) + group + s(dam, bs = "re") + s(dam, day, bs="re"), data = pups.measures, method = "REML")
summary(model.weight.gam)
plot(model.weight.gam, pages=1, se=TRUE)
anova(model.weight.gam)
emmeans(model.weight.gam, pairwise ~ group | day)


model.tail.gam <- gam(length_tail ~ s(day, by=group) + group + s(dam, bs = "re") + s(dam, day, bs="re"), data = pups.measures, method = "REML")
summary(model.tail.gam)
plot(model.tail.gam, pages=1, se=TRUE)
anova(model.tail.gam)
emmeans(model.tail.gam, pairwise ~ group | day)


model.length.gam <- gam(length_body ~ s(day, by=group) + group + s(dam, bs = "re") + s(dam, day, bs="re"), data = pups.measures, method = "REML")
summary(model.length.gam)
plot(model.length.gam, pages=1, se=TRUE)
anova(model.length.gam)
emmeans(model.length.gam, pairwise ~ group | day)


