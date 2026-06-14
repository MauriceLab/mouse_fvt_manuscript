library(tidyverse)
library(reshape2)

#### DIET OVERVIEW ####

diet <- data.frame(
  macronutrient = c("Fat", "Protein", "Carbohydrates"),
  CON = c(0.15, 0.2, 0.65),
  MAL = c(0.05, 0.07, 0.88)
)

diet.plot <- diet %>% 
  melt() %>% 
  ggplot(aes(x = variable, y = 100*value, fill = macronutrient)) +
  geom_bar(stat = "identity", color = "black") +
  ylab("% of calories") +
  scale_fill_manual(values = c("#E64B35FF", "#00A087FF", "#3C5488FF")) +
  theme_bw() +
  theme(legend.position = "right", axis.title.x=element_blank(), 
        legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


#### PLOT PUP BODY MEASURES ####
pups.measures <- read_csv("analysis-exp1/vivo6-pups-measures.csv")
pups.metadata <- read_csv("analysis-exp1/metadata-pups.csv")

pups.measures <- left_join(pups.measures, pups.metadata, by = c("litter" = "mouse"))

pups.measures$group <- as.factor(pups.measures$group)
pups.measures$postnatal.day <- as.numeric(pups.measures$date - pups.measures$litter.dob)

# pups.measures.g1 <- pups.measures %>% filter(group == "1")
# saveRDS(pups.measures.g1, file = "~/Desktop/vivo6-pups-measures-g1.rds")

graph.pups.weight <- pups.measures %>% 
  ggplot(aes(x = postnatal.day, y = as.double(weight), group = group, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 4)) +
  xlab("Age (days)") +
  ylab("Body weight (g)") +
  scale_colour_manual(values = c("#c57bea", "#779dee", "#d0ae7c"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet"), name = "Group") +
  theme_bw() +
  scale_x_continuous(breaks = c(7, 14, 21, 28, 35, 42, 49)) +
  theme(aspect.ratio=1, legend.position = "none")#, axis.title.x=element_blank())

graph.pups.weight


graph.pups.tail <- pups.measures %>% mutate(experiment = "VIVO6") %>% 
  ggplot(aes(x = postnatal.day, y = as.double(length.tail), group = group, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 2)) +
  xlab("Age (days)") +
  ylab("Tail length (mm)") +
  scale_colour_manual(values = c("#c57bea", "#779dee", "#d0ae7c"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet"), name = "Group") +
  theme_bw() +
  scale_x_continuous(breaks = c(7, 14, 21, 28, 35, 42, 49)) +
  theme(aspect.ratio=1, legend.position = "none")#, axis.title.x=element_blank())


fig.s1a <- graph.pups.tail +
  theme(aspect.ratio=1, legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2)) +
  facet_wrap(~experiment, scale = "free_y", 
             labeller = as_labeller(c(VIVO6 = "Experiment #1") ))
  
fig.s1a

graph.pups.body <- pups.measures %>% 
  ggplot(aes(x = postnatal.day, y = as.double(length.body), group = group, color = group)) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 4)) +
  xlab("Age (days)") +
  ylab("Body length (mm)") +
  scale_colour_manual(values = c("#c57bea", "#779dee", "#d0ae7c"), labels = c("H donor +\nCON diet", "S donor +\nMAL diet","H donor +\nMAL diet")) +
  theme_bw() +
  scale_x_continuous(breaks = c(7, 14, 21, 28, 35, 42, 49)) +
  theme(aspect.ratio=1, legend.position = "none")#, legend.position = "none")#, axis.title.x=element_blank())


shared_legend <- get_legend(
  graph.pups.weight + theme(legend.position = "bottom")
)

fig.1bc <- 
  plot_grid(
    plot_grid(
      graph.pups.weight, 
      graph.pups.body, 
      ncol = 2
      ),
    shared_legend,
    nrow = 2,
    rel_heights = c(1, 0.2))

fig.1bc


# modelling of pup body measures
library(mgcv)
library(emmeans)
library(gratia)

model.weight.gam <- gam(as.double(weight) ~ s(postnatal.day, by=group) + group + s(litter, bs = "re") + s(litter, postnatal.day, bs = "re"), 
                        data = pups.measures, method = "REML")
summary(model.weight.gam)
plot(model.weight.gam, pages=1, se=TRUE)

emms_weight <- emmeans(model.weight.gam, pairwise ~ group | postnatal.day, 
                       at = list(postnatal.day = seq(6, 45, by = 2)), 
                       type = "response")
as.data.frame(emms_weight$contrasts) %>% 
  #filter(p.value < 0.05) %>% 
  view()


library(gratia)
pups.measures %>% 
  group_by(group, postnatal.day) %>% 
  summarise(mean.weight = mean(weight),
            mean.tail = mean(length.tail),
            mean.body = mean(length.body)) %>% 
  filter(26 < postnatal.day) %>% 
  filter(postnatal.day < 34) %>% view()


model.tail.gam <- gam(as.double(length.tail) ~ s(postnatal.day, by=group) + group + s(litter, bs = "re"), 
                      data = pups.measures, method = "REML")
summary(model.tail.gam)
plot(model.tail.gam, pages=1, se=TRUE)

emms_tail <- emmeans(model.tail.gam, pairwise ~ group | postnatal.day, 
                       at = list(postnatal.day = seq(6, 45, by = 2)), 
                       type = "response")
as.data.frame(emms_tail$contrasts) %>% 
  #filter(p.value < 0.05) %>% 
  view()

model.length.gam <- gam(as.double(length.body) ~ s(postnatal.day, by=group) + group + s(litter, bs = "re") + s(litter, postnatal.day, bs="re"), data = pups.measures, method = "REML")
summary(model.length.gam)
plot(model.length.gam, pages=1, se=TRUE)
anova(model.length.gam)
emmeans(model.length.gam, pairwise ~ group | postnatal.day)

emms_length <- emmeans(model.length.gam, pairwise ~ group | postnatal.day, 
                     at = list(postnatal.day = seq(6, 45, by = 2)), 
                     type = "response")
as.data.frame(emms_length$contrasts) %>% 
  #filter(p.value < 0.05) %>% 
  view()
