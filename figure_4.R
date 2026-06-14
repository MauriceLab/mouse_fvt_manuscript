library(tidyverse)
library(ggVennDiagram)
library(cowplot)

# Stacked bar charts (abundance)
abundance_data <- data.frame(
  Experiment = c("Experiment #2", "Experiment #2", "Experiment #3", "Experiment #3"),
  Category = c("Host assigned", "Host unassigned", "Host assigned", "Host unassigned"),
  RelAbund = c(95, 5, 11, 89),
  Count = c(58, 29, 47, 36)
)

abundance_data$Category <- factor(abundance_data$Category, levels = c("Host unassigned", "Host assigned"))
abundance_data$Label <- paste0("(n = ", abundance_data$Count, ")")

panel_a <- ggplot(abundance_data, aes(x = Experiment, y = RelAbund/100, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = Label, color = Category), 
            position = position_stack(vjust = 0.5), 
            size = 3, 
            show.legend = FALSE) +
  scale_fill_manual(values = c("Host unassigned" = "#E0E0E0", "Host assigned" = "#2C7FB8")) +
  scale_color_manual(values = c("Host unassigned" = "black", "Host assigned" = "white")) +
  #theme_classic() +
  theme_bw() +
  labs(y = "Relative abundance in FVT (%)", fill = NULL) +
  theme(
    text = element_text(size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top"
  ) +
  scale_y_continuous(labels = scales::percent)

panel_a

# Venn diagrams

# Experiment #2 Data
exp2_fvt_predicted <- c("Bacteroides", "Bifidobacterium", "Clostridium", "Collinsella", 
                        "Companilactobacillus", "Faecalibacterium", "Haemophilus", "Pauljensenia", 
                        "Phocaeicola", "Prevotella", "Streptococcus")
exp2_pup_baseline <- c("Bacteroides", "Mediterraneibacter", "Enterococcus", "Erysipelatoclostridium", 
                       "Eubacterium", "Escherichia", "Klebsiella", "Bifidobacterium", 
                       "Anaerostipes", "Proteus", "Gordonibacter", "Anaeromassilibacillus", 
                       "Clostridium", "Neglecta", "Lactobacillus", "Hungatella", "Lachnoclostridium", 
                       "Duncaniella", "Laceyella", "Shigella", "Clostridioides", 
                       "Flavonifractor", "Anaerotruncus")

exp2_venn_list <- list(
  "FVT" = exp2_fvt_predicted,
  "Pups" = exp2_pup_baseline
)

# Experiment #3 Data
exp3_fvt_predicted <- c("Bifidobacterium", "Prevotella", "Collinsella", "HGM10836", 
                        "CAG-279", "Bacteroides", "Mycobacterium")
exp3_pup_baseline <- c("Streptococcus", "Anaerostipes", "Hungatella", "Romboutsia", 
                       "Staphylococcus", "Faecalicoccus", "Gordonibacter", "Enterococcus", 
                       "Lactococcus", "Terrisporobacter", "Paeniclostridium", "Erysipelatoclostridium", 
                       "Ruthenibacterium", "Clostridium", "Enterobacter", "Anaerotignum", 
                       "Bifidobacterium", "Turicibacter", "Butyricimonas", "Lachnoclostridium", 
                       "Pseudoflavonifractor", "Parabacteroides", "Enterocloster", 
                       "Blautia", "Bacteroides", "Flintibacter", "Flavonifractor", "Lacrimispora", 
                       "Mediterraneibacter", "Phocaeicola", "Laceyella", "Faecalitalea", 
                       "Escherichia", "Slackia", "Eggerthella")

exp3_venn_list <- list(
  "FVT" = exp3_fvt_predicted,
  "Pups" = exp3_pup_baseline
)

panel_b1 <- ggVennDiagram(exp2_venn_list, label = "count", label_alpha = 0, edge_size = 0.5, set_size = 4) +
  scale_fill_gradient(low = "white", high = "#92C5DE") + 
  theme(legend.position = "none", 
        plot.title = element_text(size=12),
        plot.margin = margin(t = 15, r = 5, b = 5, l = -100, unit = "pt")) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  #scale_y_continuous(expand = expansion(mult = 0.2)) +
  ggtitle("Experiment #2 host availability (genus-level)")

panel_b1


panel_b2 <- ggVennDiagram(exp3_venn_list, label = "count", label_alpha = 0, edge_size = 0.5, set_size = 4) +
  scale_fill_gradient(low = "white", high = "#92C5DE") + 
  theme(legend.position = "none", 
        plot.title = element_text(size=12),
        plot.margin = margin(t = 15, r = 5, b = 5, l = -100, unit = "pt")) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  #scale_y_continuous(expand = expansion(mult = 0.2)) +
  ggtitle("Experiment #3 host availability (genus-level)")


panel_b2


# Lifecycle

exp2.lifecycle <- read_tsv("vivo7-bacphlip.tsv") %>% 
  mutate(lifecycle = ifelse(temperate > 0.5, "temperate", "non-temperate")) %>% 
  filter(virus %in% c("NODE_4_length_53798_cov_5011.541131", "NODE_1_length_98866_cov_307.607594", 
                      "NODE_49_length_7768_cov_112.421107", "NODE_12_length_33149_cov_93.647217", 
                      "NODE_15_length_25124_cov_76.089912", "NODE_19_length_18409_cov_46.580691", 
                      "NODE_32_length_10340_cov_43.589499", "NODE_344_length_2139_cov_42.039347", 
                      "NODE_10_length_40318_cov_25.031642", "NODE_26_length_11514_cov_23.622393", 
                      "NODE_2_length_65832_cov_23.005367", "NODE_711_length_1467_cov_22.435552", 
                      "NODE_92_length_5255_cov_17.769231", "NODE_23_length_16638_cov_20.790870", 
                      "NODE_64_length_16604_cov_16.296695", "NODE_16_length_23040_cov_13.240505", 
                      "NODE_20_length_18213_cov_12.170448|provirus_1_13875", "NODE_14_length_25566_cov_11.862765", 
                      "NODE_38_length_46787_cov_228.578148", "NODE_148_length_3654_cov_11.196443", 
                      "NODE_16_length_27554_cov_174.849813", "NODE_214_length_2934_cov_7.283084", 
                      "NODE_22_length_17082_cov_8.175838", "NODE_749_length_1427_cov_7.816327", 
                      "NODE_726_length_1456_cov_7.029979", "NODE_72_length_6345_cov_7.043561", 
                      "NODE_129_length_4038_cov_6.656289", "NODE_46_length_8021_cov_6.572307", 
                      "NODE_524_length_1713_cov_6.638721", "NODE_53_length_7445_cov_6.344655", 
                      "NODE_95_length_5065_cov_7.959681", "NODE_40_length_8342_cov_5.657053", 
                      "NODE_450_length_1835_cov_5.998876", "NODE_272_length_2532_cov_5.148567", 
                      "NODE_1263_length_1121_cov_4.600375", "NODE_163_length_3442_cov_4.336286", 
                      "NODE_180_length_3287_cov_4.732983", "NODE_36_length_9075_cov_4.509534", 
                      "NODE_50_length_7713_cov_4.500522", "NODE_89_length_5439_cov_4.326152", 
                      "NODE_94_length_5224_cov_4.401625", "NODE_187_length_3184_cov_3.568872", 
                      "NODE_160_length_3481_cov_3.781086", "NODE_188_length_3170_cov_3.798716", 
                      "NODE_210_length_2952_cov_3.932689", "NODE_230_length_2867_cov_4.246444", 
                      "NODE_256_length_2665_cov_3.679310", "NODE_385_length_1994_cov_3.674059", 
                      "NODE_633_length_1551_cov_3.543449", "NODE_736_length_1448_cov_4.040919", 
                      "NODE_764_length_1408_cov_3.855137", "NODE_80_length_3898_cov_3.087952", 
                      "NODE_1017_length_1228_cov_3.350384", "NODE_1146_length_1173_cov_3.042039", 
                      "NODE_1423_length_1061_cov_3.210736", "NODE_152_length_3616_cov_3.334457", 
                      "NODE_236_length_2821_cov_3.253073", "NODE_237_length_2809_cov_3.084241", 
                      "NODE_251_length_2687_cov_3.181611", "NODE_254_length_2668_cov_3.052047", 
                      "NODE_297_length_2403_cov_3.211670", "NODE_362_length_2084_cov_3.660917", 
                      "NODE_367_length_2068_cov_3.520616", "NODE_935_length_1270_cov_3.364609", 
                      "NODE_990_length_1241_cov_3.139966", "NODE_161_length_5158_cov_6.566921", 
                      "NODE_1461_length_1047_cov_2.663306", "NODE_1558_length_1010_cov_2.902618", 
                      "NODE_268_length_2575_cov_2.801190", "NODE_295_length_2417_cov_3.035563", 
                      "NODE_343_length_2140_cov_2.683933", "NODE_478_length_1783_cov_2.836806", 
                      "NODE_527_length_1708_cov_2.523291", "NODE_944_length_1267_cov_3.127063", 
                      "NODE_115_length_6212_cov_4.749229", "NODE_1501_length_1032_cov_1.957011", 
                      "NODE_205_length_4285_cov_4.502128", "NODE_184_length_3248_cov_2.385218", 
                      "NODE_269_length_2575_cov_2.501587", "NODE_312_length_2324_cov_2.236668", 
                      "NODE_408_length_1927_cov_2.191239", "NODE_724_length_1458_cov_2.410549", 
                      "NODE_791_length_1385_cov_2.360150", "NODE_178_length_4823_cov_5.003775", 
                      "NODE_332_length_2960_cov_3.930120", "NODE_459_length_2232_cov_5.176849", 
                      "NODE_104_length_6596_cov_5.204862")) %>% 
  select(virus, lifecycle) %>% 
  summarize(n_temperate = sum(lifecycle == "temperate"), n_lytic = sum(lifecycle == "non-temperate"), 
            percent_temperate = n_temperate/n(), percent_lytic = n_lytic/n(), 
            experiment = "Experiment #2")

exp3.lifecycle <- read_tsv("vivo8-bacphlip.tsv") %>% 
  mutate(lifecycle = ifelse(temperate > 0.5, "temperate", "non-temperate")) %>% 
  filter(virus %in% c("NODE_16_length_20004_cov_2375.001805", "NODE_24_length_16285_cov_42.554713", 
                      "NODE_479_length_2068_cov_51.126180", "NODE_912_length_8852_cov_15.026486", 
                      "NODE_41_length_11398_cov_19.446619", "NODE_10_length_29026_cov_12.843982", 
                      "NODE_5_length_36112_cov_12.717004", "NODE_18_length_19544_cov_11.024270", 
                      "NODE_64_length_8589_cov_10.134521", "NODE_38_length_12381_cov_9.720509", 
                      "NODE_81_length_6522_cov_9.763260", "NODE_110_length_5458_cov_9.223024", 
                      "NODE_28_length_15120_cov_8.141653", "NODE_1387_length_1039_cov_20.370935", 
                      "NODE_44_length_11166_cov_7.982990", "NODE_131_length_4785_cov_7.429598", 
                      "NODE_46_length_10779_cov_6.092503", "NODE_79_length_11271_cov_35.365549", 
                      "NODE_109_length_5517_cov_5.794764", "NODE_449_length_1931_cov_5.841151", 
                      "NODE_1172_length_1138_cov_33.322253", "NODE_626_length_1564_cov_5.337309", 
                      "NODE_113_length_5349_cov_5.386286", "NODE_17_length_19741_cov_4.852687", 
                      "NODE_229_length_3048_cov_5.127631", "NODE_79_length_6809_cov_5.143471", 
                      "NODE_95_length_5968_cov_4.964654", "NODE_176_length_3854_cov_4.569624", 
                      "NODE_630_length_1562_cov_5.102853", "NODE_103_length_5688_cov_4.401740", 
                      "NODE_1076_length_1174_cov_4.462913", "NODE_122_length_4946_cov_4.345328", 
                      "NODE_180_length_3805_cov_4.354133", "NODE_184_length_3654_cov_4.500695", 
                      "NODE_200_length_3419_cov_4.433413", "NODE_256_length_2773_cov_4.676233", 
                      "NODE_621_length_1567_cov_4.782407", "NODE_940_length_1266_cov_4.521057", 
                      "NODE_98_length_5873_cov_4.399966", "NODE_349_length_2569_cov_4.017900", 
                      "NODE_114_length_5287_cov_4.093272", "NODE_115_length_5265_cov_4.248752", 
                      "NODE_119_length_5193_cov_4.218568", "NODE_153_length_4167_cov_3.723492", 
                      "NODE_215_length_3236_cov_4.272870", "NODE_217_length_3222_cov_3.829807", 
                      "NODE_248_length_2813_cov_4.138144", "NODE_72_length_7180_cov_3.890667", 
                      "NODE_84_length_6403_cov_3.902804", "NODE_88_length_6325_cov_3.779426", 
                      "NODE_99_length_5851_cov_3.793651", "NODE_1077_length_1174_cov_3.741734", 
                      "NODE_1322_length_1047_cov_3.806452", "NODE_160_length_4008_cov_3.618771", 
                      "NODE_246_length_2849_cov_3.388332", "NODE_266_length_2659_cov_3.225422", 
                      "NODE_274_length_2619_cov_3.736349", "NODE_358_length_2209_cov_3.369545", 
                      "NODE_53_length_9847_cov_3.295752", "NODE_716_length_1459_cov_3.289174", 
                      "NODE_186_length_3606_cov_2.664038", "NODE_231_length_3028_cov_2.509250", 
                      "NODE_280_length_2575_cov_2.821825", "NODE_585_length_1621_cov_2.642401", 
                      "NODE_591_length_1610_cov_2.979421", "NODE_819_length_1362_cov_2.496557", 
                      "NODE_994_length_4515_cov_8.534753", "NODE_1059_length_1184_cov_2.370239", 
                      "NODE_1127_length_1147_cov_2.521978", "NODE_1174_length_1125_cov_1.920561", 
                      "NODE_1225_length_1099_cov_2.553640", "NODE_1321_length_1048_cov_2.348439", 
                      "NODE_1338_length_1041_cov_2.307302", "NODE_1382_length_1024_cov_2.011352", 
                      "NODE_1411_length_1012_cov_2.218391", "NODE_195_length_3486_cov_2.403964", 
                      "NODE_201_length_3412_cov_2.440870", "NODE_770_length_1403_cov_2.066766", 
                      "NODE_809_length_1367_cov_2.485518", "NODE_825_length_1360_cov_2.013027", 
                      "NODE_870_length_1319_cov_1.970728", "NODE_735_length_1436_cov_1.792904", 
                      "NODE_793_length_1381_cov_1.760935")) %>% 
  select(virus, lifecycle) %>% 
  summarize(n_temperate = sum(lifecycle == "temperate"), n_lytic = sum(lifecycle == "non-temperate"), 
            percent_temperate = n_temperate/n(), percent_lytic = n_lytic/n(), 
            experiment = "Experiment #3")


lifecycle_long <- rbind(exp2.lifecycle, exp3.lifecycle) %>% 
  pivot_longer(
    cols = -experiment,
    names_to = c(".value", "lifecycle"),
    names_sep = "_"
  ) %>% 
  mutate(lifecycle = gsub("lytic", "non-temperate", lifecycle)) %>% 
  mutate(
    lifecycle = str_to_title(lifecycle),
    lifecycle = factor(lifecycle, levels = c("Temperate", "Non-Temperate"))
  )

lifecycle_plot <- ggplot(lifecycle_long, aes(x = experiment, y = percent, fill = lifecycle)) +
  geom_col(color = "black", width = 0.6) +
  geom_text(aes(label = paste0("(n = ", n, ")"), color = lifecycle), 
            position = position_stack(vjust = 0.5),
            size = 4, show.legend = F) +
  # Pick colors that match your paper's aesthetic (e.g., grey and blue)
  scale_fill_manual(values = c("Temperate" = "grey80", "Non-Temperate" = "#2166AC")) +
  scale_color_manual(values = c("Temperate" = "black", "Non-Temperate" = "white")) +
  labs(y = "Proportion of FVT vOTUs", fill = "Predicted lifecycle") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 13, color = "black"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    axis.title.x = element_blank(),
  ) +
  scale_y_continuous(labels = scales::percent)

lifecycle_plot

save_plot("~/Desktop/FigureS4.pdf", lifecycle_plot, dpi = 600, base_height = 4, base_width = 6)

# Combine panels
fig.4 <- plot_grid(
  panel_a,
  plot_grid(
    panel_b1,
    panel_b2,
    ncol = 1
  ),
  ncol = 2,
  rel_widths = c(0.75, 0.8), 
  labels = c("A", "B")
)

save_plot("~/Desktop/Figure4.pdf", fig.4, dpi = 600, base_height = 4, base_width = 8)




