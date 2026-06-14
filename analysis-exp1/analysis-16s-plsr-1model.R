library(pls)
library(phyloseq)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

#### COMPARE EFFECTS OF DIET - TRAIN ON GROUP 1, TEST ON GROUP 3 ####
donor1.group1.pups.metadata <- metadata %>% 
  filter(Group == 1, SampleType == "Pup")


#ordinate with group 1 & group 3 (donor 1 only), then train PLSR
ord.pcoa.uni <- ordinate(ps.relabund.pups %>% subset_samples(Group %in% c(1, 3)) %>% subset_samples(DonorID == 1), method="PCoA", distance="wunifrac")

ord.pcoa.uni.pcs <- ord.pcoa.uni$vectors %>% as.data.frame()
ord.pcoa.uni.pcs$SampleName <- rownames(ord.pcoa.uni.pcs)

ord.pcoa.uni.pcs <- left_join(ord.pcoa.uni.pcs, 
                              metadata %>% select(SampleName, PostnatalDay))

ord.pcoa.uni.pcs <- ord.pcoa.uni.pcs %>% filter(SampleName %in% donor1.group1.pups.metadata$SampleName)

rownames(ord.pcoa.uni.pcs) <- ord.pcoa.uni.pcs$SampleName
ord.pcoa.uni.pcs <- ord.pcoa.uni.pcs %>% subset(select = -c(SampleName))

set.seed(123)

plsr.age <- plsr(PostnatalDay~., data = ord.pcoa.uni.pcs, validation = "CV")
summary(plsr.age)

validationplot(plsr.age, xlim = c (0, 15), ylim=c(0, 15), type = "b")
plot(cumsum(explvar(plsr.age)))
plot(drop(R2(plsr.age, estimate = "train", intercept = FALSE)$val))

best.dims <- 5 #which.min(RMSEP(plsr.age)$val[1,,])-1

plsr.age <- plsr(PostnatalDay~., data = ord.pcoa.uni.pcs, validation = "CV", ncomp = best.dims)


#so let's use this model on our training data
plsr.train <- predict(plsr.age, ord.pcoa.uni.pcs[, 1:ncol(ord.pcoa.uni.pcs)-1], ncomp = best.dims)

plsr.train.df <- cbind(plsr.train, ord.pcoa.uni.pcs %>% select(PostnatalDay))
colnames(plsr.train.df) <- c("PostnatalDay_pred", "PostnatalDay")
plsr.train.df$dataset <- "Group 1 (train)"

mylims.train <- range(with(plsr.train.df, c(PostnatalDay, PostnatalDay_pred)))
plsr.train.df %>% 
  ggplot(aes(x = PostnatalDay, y = PostnatalDay_pred)) +
  geom_point() +
  coord_cartesian(xlim = mylims.train, ylim = mylims.train) + 
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 2)) +
  theme(aspect.ratio=1)


#test group 3 (donor 1 only) using model
test.ordination.plsr.g3 <- ord.pcoa.uni$vectors %>% as.data.frame()

test.ordination.plsr.g3$SampleName <- rownames(test.ordination.plsr.g3)

test.ordination.plsr.g3 <- test.ordination.plsr.g3 %>% 
  filter(SampleName %in% (metadata %>% filter(Group == 3, DonorID == 1, SampleType == "Pup"))$SampleName)

test.ordination.plsr.g3 <- test.ordination.plsr.g3 %>% select(-SampleName)

plsr.test.g3 <- predict(plsr.age, test.ordination.plsr.g3, ncomp = best.dims)

plsr.test.metadata.g3 <- ps.relabund.pups %>% subset_samples(Group == 3) %>% subset_samples(DonorID == 1) %>% sample_data()
class(plsr.test.metadata.g3) <- "data.frame"

plsr.test.df.g3 <- cbind(plsr.test.g3, plsr.test.metadata.g3 %>% select(PostnatalDay))
colnames(plsr.test.df.g3) <- c("PostnatalDay_pred", "PostnatalDay")
plsr.test.df.g3$dataset <- "Group 3 (test)"

plsr.test.df.g3 %>% 
  ggplot(aes(x = PostnatalDay, y = PostnatalDay_pred)) +
  geom_point() +
  theme(aspect.ratio=1) +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 2))

 
# #test group 2 using model
# test.ordination.plsr.g2 <- ord.pcoa.uni$vectors %>% as.data.frame()
# 
# test.ordination.plsr.g2$SampleName <- rownames(test.ordination.plsr.g2)
# 
# test.ordination.plsr.g2 <- test.ordination.plsr.g2 %>% filter(SampleName %in% (metadata %>% filter(Group == 2))$SampleName)
# 
# test.ordination.plsr.g2 <- test.ordination.plsr.g2 %>% select(-SampleName)
# 
# plsr.test.g2 <- predict(plsr.age, test.ordination.plsr.g2, ncomp = best.dims)
# 
# plsr.test.metadata.g2 <- ps.relabund.pups %>% subset_samples(Group == 2) %>% sample_data()
# class(plsr.test.metadata.g2) <- "data.frame"
# 
# plsr.test.df.g2 <- cbind(plsr.test.g2, plsr.test.metadata.g2 %>% select(PostnatalDay))
# colnames(plsr.test.df.g2) <- c("PostnatalDay_pred", "PostnatalDay")
# plsr.test.df.g2$dataset <- "Group 2 (test)"
# 
# plsr.test.df.g2 %>% 
#   ggplot(aes(x = PostnatalDay, y = PostnatalDay_pred)) +
#   geom_point() +
#   theme(aspect.ratio=1) +
#   geom_smooth(method = "lm", se = F, formula = y ~ splines::bs(x, 2))


model1.plsr.to.plot <- rbind(plsr.train.df, plsr.test.df.g3)
mylims <- range(with(model1.plsr.to.plot, c(PostnatalDay, PostnatalDay_pred)))

model1.test.plot <- model1.plsr.to.plot
model1.test.plot$SampleName <- rownames(model1.test.plot)

#add donor ID info
model1.test.plot <- model1.test.plot %>%
  left_join(metadata %>% select(SampleName, DonorID))

model1.test.plot

fig.2c <- model1.test.plot %>% 
  ggplot(aes(x = PostnatalDay, y = PostnatalDay_pred, color = interaction(dataset, DonorID), group = interaction(dataset, DonorID))) +
  # geom_rect(aes(xmin = age.kmeans.ranges[1, 2]-10, xmax = age.kmeans.ranges[1, 3]+1, ymin = -Inf, ymax = Inf),
  #           fill = "#f9f0e2", alpha = 0.1, color = NA) +
  # geom_rect(aes(xmin = age.kmeans.ranges[2, 2]-1, xmax = age.kmeans.ranges[2, 3]+1, ymin = -Inf, ymax = Inf),
  #           fill = "#f0d9b8", alpha = 0.1, color = NA) +
  # geom_rect(aes(xmin = age.kmeans.ranges[3, 2]-1, xmax = age.kmeans.ranges[3, 3]+10, ymin = -Inf, ymax = Inf),
  #           fill = "#e7c28e", alpha = 0.1, color = NA) +
  geom_point(alpha = 0.5, position = "jitter") +
  xlab("Actual age (days)") +
  ylab("Predicted age (days)") +
  #ggtitle("Bacteriome - PLSR analysis", "model #1, dist = wuni, ncomp = 3") +
  coord_cartesian(xlim = mylims, ylim = mylims) +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 2), linewidth = 1.3) +
  scale_color_manual(values = c("#c57bea", "#d0ae7c"), name = "Group", labels = c("H donor + CON diet (donor 1)", "H donor + MAL diet (donor 1)")) +
  theme_bw() +
  theme(aspect.ratio=1)

fig.2c

### PERFORMANCE METRICS ###
test_rmse <- sqrt(mean((plsr.test.df.g3$PostnatalDay_pred - plsr.test.df.g3$PostnatalDay)^2)) # 11.29 days
test_r2 <- cor(plsr.test.df.g3$PostnatalDay_pred, plsr.test.df.g3$PostnatalDay)^2 # R2 = 0.054, r = 0.232
train_rmse <- sqrt(mean((plsr.train.df$PostnatalDay_pred - plsr.train.df$PostnatalDay)^2)) # 2.97 days
train_r2 <- cor(plsr.train.df$PostnatalDay_pred, plsr.train.df$PostnatalDay)^2 # R2 = 0.803, r = 0.896

### GAMM
library(mgcv)
library(emmeans)

model1.test.plot$dataset <- as.factor(model1.test.plot$dataset)

model.plsr.gam <- gam(PostnatalDay_pred ~ s(PostnatalDay, by=dataset) + dataset, data = model1.test.plot, method = "REML")
summary(model.plsr.gam)
plot(model.plsr.gam, pages=1, se=TRUE)
anova(model.plsr.gam)
emms_microbiome_age <- emmeans(model.plsr.gam, pairwise ~ dataset | PostnatalDay, 
                       at = list(PostnatalDay = seq(20, 45, by = 2)), 
                       type = "response")
as.data.frame(emms_microbiome_age$contrasts) %>% 
  #filter(p.value < 0.05) %>% 
  view()

### COEFFICIENTS START, MODEL 1 (GROUPS 1 AND 3, DONOR 1 ONLY) ###

coefficients <- as.data.frame(coef(plsr.age))
coefficients$axis <- rownames(coefficients)

colnames(coefficients) <- c("plsr.coeff", "orig.axis")

coefficients %>% 
  mutate(orig.axis = fct_reorder(orig.axis, plsr.coeff)) %>%
  ggplot(aes(x = orig.axis, y = plsr.coeff)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

coefficients.mat <- coefficients %>% 
  reshape(idvar = "plsr.coeff", timevar = "orig.axis", direction = "wide") %>% 
  t()

loadings <- plot_ordination(ps.relabund.pups %>% subset_samples(Group %in% c(1, 3)) %>% subset_samples(DonorID == 1), ord.pcoa.uni, axes = 1:28, justDF = T, type = "biplot") %>%
  filter(id.type == "Taxa") %>% 
  select(coefficients$orig.axis, species)

loadings.trf <- t(sapply(1:nrow(loadings), function(i) loadings[i,1:ncol(loadings)-1] * coefficients.mat)) %>% 
  apply(2, function(x) sapply(x, `[`, 1))
loadings.trf <- cbind(loadings.trf, loadings$species)

loadings.trf <- loadings.trf[apply(loadings.trf[, 1:ncol(loadings.trf)-1], 1, function(x) all(x != "NaN")), ]

loadings.trf <- as.data.frame(loadings.trf)
#loadings.trf[, 1:ncol(loadings.trf)-1] <- loadings.trf[, 1:ncol(loadings.trf)-1] %>% apply(2, as.double) %>% as.data.frame()
names(loadings.trf)[ncol(loadings.trf)] <- "species"

loadings.trf <- cbind(as.data.frame(loadings.trf[, 1:ncol(loadings.trf)-1] %>% apply(2, as.numeric)), loadings.trf$species)
names(loadings.trf)[ncol(loadings.trf)] <- "species"

loadings.trf.mag <- as.data.frame(loadings.trf[, 1:ncol(loadings.trf)-1] %>% apply(1, function(x) norm(x, type = "2")))
loadings.trf.mag <- cbind(loadings.trf.mag, loadings.trf$species)
colnames(loadings.trf.mag) <- c("magnitude", "species")

### COEFFICIENTS END ###

### HEATMAP OF IMPORTANT TAXA FOR HEALTHY DEVELOPMENT ###
library(pheatmap)

taxtable.genus <- taxtable %>% 
  group_by(genus) %>% 
  summarise(magnitude = sum(magnitude))

ps.relabund.summary.genus <- ps.relabund.pups %>% 
  subset_samples(Group %in% c(1, 3)) %>% 
  subset_samples(DonorID == 1) %>% 
  psmelt() %>%                    
  filter(Abundance > 0) %>% 
  filter(genus %in% taxtable.genus$genus) %>% 
  left_join(taxtable.genus) %>% 
  group_by(Sample, genus) %>% 
  summarise(Group = Group, Abundance = sum(Abundance), PostnatalDay = PostnatalDay, magnitude = magnitude) %>% 
  distinct() %>% 
  ungroup() %>% 
  group_by(Group, genus) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  filter(!is.na(genus))

heatmap.g1 <- ps.relabund.summary.genus %>% 
  filter(Group == 1) %>% 
  group_by(PostnatalDay, genus) %>% 
  summarise(mean = mean(Abundance), PostnatalDay = PostnatalDay) %>% 
  distinct() %>% 
  #spread(key = PostnatalDay, value = mean) %>%
  pivot_wider(values_from = mean, names_from = PostnatalDay) %>% 
  as.data.frame()
rownames(heatmap.g1) <- heatmap.g1[, 1]
heatmap.g1 <- heatmap.g1[, -1]

heatmap.g1[is.na(heatmap.g1)] <- 0

#### calculate R2 for PHFs ####
important.taxa.r2 <- ps.relabund.summary.genus %>% 
  filter(Group == 1) %>% 
  select(genus, Abundance, PostnatalDay) %>% 
  group_by(genus) %>% 
  summarize(cor = sign(cor(Abundance, PostnatalDay)),
            r2 = cor(Abundance, PostnatalDay)^2,
            sign_r2 = r2*cor)

important.taxa <- ps.relabund.summary.genus %>% 
  left_join(important.taxa.r2, by = c("genus" = "genus")) #%>% 
  #top_n(20, magnitude)

max_abs_R2 <- max(abs(important.taxa$sign_r2), na.rm = TRUE)
symm_limits <- c(-max_abs_R2, max_abs_R2)

neg_colors <- viridisLite::mako(100,begin = 0.3)
pos_colors <- viridisLite::magma(100,begin = 0.6)
all_colors <- c(neg_colors, rev(pos_colors))

col_r2 = colorRamp2(
  seq(symm_limits[1], symm_limits[2], length.out = length(all_colors)), 
  all_colors
)

lgd_r2 = Legend(
  title = bquote(bold("sign * R"^2)),
  col_fun = col_r2, 
  at = c(symm_limits[1], 0, symm_limits[2]), 
  labels = c(round(symm_limits[1], 2), "0", round(symm_limits[2], 2))
)

heatmap_colors = viridisLite::viridis(100, option = "viridis")

lgd_abundance = Legend(
  title = "Scaled relative\nabundance", 
  col_fun = colorRamp2(seq(0, 1, length.out = 100), heatmap_colors), 
  at = c(0, 1), 
  labels = c("min", "max"),
  direction = "vertical"
)

heatmap_scaled <- t(apply(heatmap.g1, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
row_data <- important.taxa[match(rownames(heatmap_scaled), important.taxa$genus), ]

row_anno = rowAnnotation(
  "sign * R2" = anno_simple(row_data$sign_r2, col = col_r2, width = unit(5, "mm"), border = T),
  
  "spacer" = anno_empty(border = FALSE, width = unit(0.5, "mm")),
  
  "Model importance" = anno_barplot(
    row_data$magnitude,
    baseline = 0,
    axis_param = list(side = "bottom"),
    gp = gpar(fill = "#D1D1D1", col = "white"),
    width = unit(3, "cm")
  ),
  annotation_label = list(
    "sign * R2" = expression("sign * R"^2)
  ),
  show_annotation_name = c("sign * R2" = TRUE, "spacer" = FALSE, "Model importance" = TRUE)
)


CombinedPlot <- Heatmap(
  heatmap_scaled,
  col = heatmap_colors,
  show_heatmap_legend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  right_annotation = row_anno,
  column_names_side = "bottom",
  column_title_side = "bottom",
  column_names_rot = 45,
  row_names_side = "left",
  row_names_gp = gpar(
    fontface = "italic"
  ),
  column_title = "Age (days)"
)

plot.heatmap <- ComplexHeatmap::draw(CombinedPlot, 
                                   annotation_legend_list = packLegend(lgd_abundance, lgd_r2), 
                                   annotation_legend_side = "right")

plot.heatmap

pdf("~/Desktop/Figure2D.pdf", height = 4, width = 7)
draw(plot.heatmap)
dev.off()

