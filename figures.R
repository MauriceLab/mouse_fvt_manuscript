library(tidyverse)
library(cowplot)
library(magick)

#### FIGURE 1 ####
exp1.model.overview <- ggdraw() + draw_image(image_read_pdf("~/Desktop/Figure1A1.pdf", density = 600))

diet.plot

fig.1a <- plot_grid(
  exp1.model.overview,
  #NULL,
  diet.plot,
  ncol = 2,
  rel_widths = c(1, 0.5),
  labels = c("A", "")
  )

shared_legend <- get_legend(
  graph.pups.weight + 
    scale_color_manual(values = c("#c57bea", "#779dee", "#d0ae7c"),
                         labels = c("Group 1\n(H donor +\nCON diet)", "Group 2\n(S donor +\nMAL diet)","Group 3\n(H donor +\nMAL diet)"),
                       name = "Group") +
    theme(legend.position = "bottom")
)

fig.1bc <- 
  plot_grid(
    plot_grid(
      graph.pups.weight, 
      graph.pups.body, 
      fig.1c,
      ncol = 3,
      labels = c("B", "", "C")
    ),
    shared_legend,
    nrow = 2,
    rel_heights = c(1, 0.2))
fig.1bc

fig.1 <- plot_grid(
  fig.1a,
  fig.1bc,
  nrow = 2,
  rel_heights = c(0.9, 1)
)

fig.1

save_plot("~/Desktop/Figure1.pdf", fig.1, dpi = 600, base_height = 5, base_width = 8)


#### FIGURE 2 ####
fig.2a

fig.2b <- plot_grid(
  pcoa.wuni.age,
  pcoa.wuni.group +
    scale_color_manual(values = c("#c57bea", "#779dee", "#d0ae7c"),
                       labels = c("Grp. 1 (H donor +\nCON diet)", "Grp. 2 (S donor +\nMAL diet)","Grp. 3 (H donor +\nMAL diet)"),
                       name = "Group"),
  pcoa.wuni.donor,
  ncol = 3,
  nrow = 1
  )
fig.2b

fig.2c

fig.2d <- ggdraw() + draw_image(image_read_pdf("~/Desktop/Figure2D.pdf", density = 600))


fig.2 <- plot_grid(
  plot_grid(fig.2a +
              scale_color_manual(values = c("#c57bea", "#779dee", "#d0ae7c"),
                                 labels = c("Group 1 (H donor +\nCON diet)", "Group 2 (S donor +\nMAL diet)","Group 3 (H donor +\nMAL diet)"),
                                 name = "Group") +
              scale_fill_manual(values = c("#c57bea", "#779dee", "#d0ae7c"),
                                labels = c("Group 1 (H donor +\nCON diet)", "Group 2 (S donor +\nMAL diet)","Group 3 (H donor +\nMAL diet)"),
                                name = "Group"),
            fig.2b,
            ncol = 2,
            rel_widths = c(0.55, 1),
            labels = c("A", "B")),
  plot_grid(fig.2c +
              scale_color_manual(values = c("#c57bea", "#d0ae7c"),
                                 labels = c("Group 1\n(H donor + CON diet, donor 1)", "Group 3\n(H donor + MAL diet, donor 1)"),
                                 name = "Group"),
            fig.2d,
            ncol = 2,
            rel_widths = c(0.75, 1),
            labels = c("C", "D")),
  nrow = 2,
  rel_heights = c(1, 0.8)
)

fig.2

save_plot("~/Desktop/Figure2.pdf", fig.2, dpi = 600, base_height = 6, base_width = 12)



#### FIGURE 3 ####
exp2_3.model.overview <- ggdraw() + draw_image(image_read_pdf("~/Desktop/Figure3A.pdf", density = 600))

fig.3b

fig.3c <- plot_grid(
  fig.3c1,
  fig.3c2,
  nrow = 1
)

fig.3d

fig.3 <- plot_grid(
  plot_grid(
    #exp2_3.model.overview,
    NULL,
    fig.3b,
    nrow = 1,
    rel_widths = c(1, 0.9),
    labels = c("A", "B")
  ),
  plot_grid(
    fig.3c,
    fig.3d,
    nrow = 1,
    rel_widths = c(1, 0.5),
    labels = c("C", "D")
  ),
  nrow = 2
)

fig.3

save_plot("~/Desktop/Figure3.pdf", fig.3, dpi = 600, base_height = 7, base_width = 9)

#### SUPPL FIGURE 1 ####
fig.s1a

fig.s1b

fig.s1 <- plot_grid(
  fig.s1a +
    scale_color_manual(values = c("#c57bea", "#779dee", "#d0ae7c"),
                       labels = c("Group 1\n(H donor + CON diet)", "Group 2\n(S donor + MAL diet)", "Group 3\n(H donor + MAL diet)"),
                       name = "Group"),
  fig.s1b,
  nrow = 1,
  labels = c("A", "B"),
  rel_widths = c(0.8, 1))

fig.s1

save_plot("~/Desktop/SupplFigure1.pdf", fig.s1, dpi = 600, base_height = 3, base_width = 7.5)



#### SUPPL FIGURE 2 ####

fig.s2a1
fig.s2a2
fig.s2b
fig.s2c

fig.s2d1
fig.s2d2
fig.s2e
fig.s2f

fig.s2.exp2.title <- ggdraw() + 
  draw_label("Experiment #2", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

fig.s2.exp3.title <- ggdraw() + 
  draw_label("Experiment #3", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

fig.s2 <- plot_grid(
  fig.s2.exp2.title,
  plot_grid(
    plot_grid(
      fig.s2a1,
      fig.s2a2,
      nrow = 1
    ),
    fig.s2b + theme(aspect.ratio = 1),
    fig.s2c + theme(aspect.ratio = 1),
    nrow = 1,
    labels = c("A", "B", "C"),
    rel_widths = c(1, 0.67, 0.67)
  ),
  fig.s2.exp3.title,
  plot_grid(
    plot_grid(
      fig.s2d1,
      fig.s2d2,
      nrow = 1
    ),
    fig.s2e + theme(aspect.ratio = 1),
    fig.s2f + theme(aspect.ratio = 1),
    nrow = 1,
    labels = c("D", "E", "F"),
    rel_widths = c(1, 0.67, 0.67)
  ),
  nrow = 4,
  rel_heights = c(0.1, 1, 0.1, 1)
)

fig.s2

save_plot("~/Desktop/SupplFigure2.pdf", fig.s2, dpi = 600, base_height = 7.5, base_width = 13)

#### SUPPL FIGURE 3 ####

fig.s3 <- plot_grid(
  plot_grid(
    fig.s3a,
    fig.s3b,
    nrow = 1,
    rel_widths = c(1, 0.7),
    labels = c("A", "B")
  ),
  plot_grid(
    fig.s3c, 
    NULL,
    nrow = 1,
    rel_widths = c(0.7, 1),
    labels = c("C", "")
  ),
  rel_heights = c(1, 0.7),
  nrow = 2
)

fig.s3

save_plot("~/Desktop/SupplFigure3.pdf", fig.s3, dpi = 600, base_height = 8, base_width = 10)

#### SUPPL FIGURE 4 ####

