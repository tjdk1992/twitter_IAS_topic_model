#-----------------------------------------------------------------------------#
# Script Name: 19-NIS-topic-dissimilarity.R                                   #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script calculate the Bray-Curtis dissimilarity among taxon           #
# according to topic distribution over NIS.                                   #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Packages
pacman::p_load(tidyverse,  # for data manipulation
               hrbrthemes, # for nice visualization
               magrittr, 　# for data manipulation
               pals,
               patchwork,  # to use color palette
               reshape2,
               rstatix,
               vegan,
               glue,
               ggtext
               )

# Data
NIS_topic <- read_csv("data/NIS-popular-topic.csv")

# Color palette
pal_orig <- pals::cols25(7)[c(5, 4, 7, 2, 1, 6, 3)]

# Levels of taxonomic group
taxon_arrange <- c("Mammal", "Bird", "Reptile", "Amphibian", "Fish", "Invertebrate", "Plant")

# Topicの類似度(pairwise) -----------------------------------------------------

# 種レベルのトピック分布
df_dissim <- NIS_topic %>% 
  dplyr::select(TP01:TP25) %>% 
  as.data.frame()
rownames(df_dissim) <- NIS_topic$name_sp

# Calculate Bray-Curtis dissimilarlity
df_bray <- df_dissim %>% 
  vegdist(method = "bray", upper = TRUE) %>% 
  as.matrix() %>% 
  melt(varnames = c("var1", "var2")) %>% 
  as_tibble() %>% 
  mutate(target_rm = if_else(var1 == var2, "T", "F")) %>% 
  filter(target_rm != "T") %>% 
  dplyr::select(-target_rm)

# Remove combination between different taxonomimc group
df_bray_pairwise <- df_bray %>% 
  left_join(
    NIS_topic %>% 
      transmute(var1 = name_sp, group_biol_1 = group_biol), by = "var1") %>% 
  left_join(
    NIS_topic %>% 
      transmute(var2 = name_sp, group_biol_2 = group_biol), by = "var2") %>% 
  mutate(check = if_else(group_biol_1 == group_biol_2, "T", "F")) %>% 
  filter(check == "T") %>% 
  transmute(var1, var2, value, group_biol = group_biol_1)

# Remove doubled duplicate
df_bray_sp <- data.frame()
for (i in 1:nrow(df_bray_pairwise)) {
  t_df_bray_sp <- df_bray_pairwise[i, ]
  t_df_bray_sp$id_distinct <- paste(
    str_sort(c(t_df_bray_sp$var1, t_df_bray_sp$var2)), collapse = "_")
  df_bray_sp <- rbind(df_bray_sp, t_df_bray_sp)
}

df_bray_sp %<>% 
  distinct(id_distinct, .keep_all = TRUE)

# count n in each group
n_ias_group <- df_bray_pairwise %>% 
  group_by(group_biol) %>% 
  summarise(n = n())

label_n_ias_group <- n_ias_group %>% 
  mutate(group_biol = str_to_title(group_biol),
                            group_biol = factor(group_biol, levels = taxon_arrange)) %>% 
  arrange(group_biol) %>% 
  transmute(n_group_biol = str_c(group_biol, "\n(", n, ")")) %>% 
  pull()

# Plot
df_bray_sp %>% 
  left_join(n_ias_group, by = "group_biol") %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange)) %>% 
  ggplot() +
  geom_boxplot(aes(x = group_biol, y = value, fill = group_biol), 
               show.legend = FALSE) +
  coord_flip() +
  # geom_jitter(aes(x = group_biol, y = value), width = 0.3, alpha = 0.3) +
  scale_x_discrete(labels = label_n_ias_group) +
  labs(x = "", 
       y = "Bray-Curtis dissimilarity") +
  scale_fill_manual(values = pal_orig) + # cols25かalphabet2のどちらかが良さそう。
  theme_ipsum(base_size = 10,
              axis_title_size = 10,
              strip_text_size = 10,
              axis_text_size = 8,
              axis_title_just = "center",
              base_family = "Helvetica",
              plot_margin = margin(5, 5, 5, 5)) + 
  theme(legend.position = "right",
        legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.key.size = unit(5, 'mm'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Save the visualized result
ggsave("fig/boxplot_sp-pairwise-topic-dissimilarity.png",
       units = "mm", width = 84, height = 60)
ggsave("fig/boxplot_sp-pairwise-topic-dissimilarity.eps",
       units = "mm", width = 84, height = 60, device = cairo_ps)

# Distribution of dissimilarlity
df_bray_sp %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(. ~ group_biol)

# Multiple comparison by Kruskal-Wallis method
# (res.kruskal <- kruskal_test(df_bray_sp, value ~ group_biol))
(res_anova <- anova_test(df_bray_sp, value ~ group_biol))
summary(res_anova)
# Multiple pairwise-comparisons
# df_bray_sp %>% 
#   dunn_test(value ~ group_biol, p.adjust.method = "bonferroni") %>% 
#   as.data.frame()
res_anova <- aov(value ~ group_biol, data = df_bray_sp) %>% 
  tukey_hsd() %>% 
  as.data.frame()
# All combination patterns in Tukey HSD
res_anova_t <- res_anova %>% 
  mutate(estimate = -1*estimate,
         conf.low = -1*conf.low,
         conf.high = -1*conf.high) %>% 
  rename(group1 = group2, group2 = group1)

res_anova <- res_anova %>% 
  bind_rows(res_anova_t) %>% 
  mutate(comb = str_c(group1, "-", group2)) %>% 
  as_tibble()

vis_tukey <- list()
l_biol <- c("mammal", "bird", "reptile", "amphibian", 
            "fish", "invertebrate", "plant")
for (i in 1:length(l_biol)) {
  # taxonomic group
  biol <- l_biol[i]
  # df of the taxonomic group
  target_res_anova <- res_anova %>% 
    filter(group1 == biol) %>% 
    mutate(signif = if_else(p.adj.signif == "ns", "ns", "signif"),
           signif = factor(signif, 
                           levels = c("signif", "ns")))
  # check whether there are any significant combination
  n_signif <- target_res_anova %>% 
    group_by(signif) %>% 
    summarise(n = n()) %>% 
    nrow()
  # Visualization
  if (n_signif == 1) {
    vis_tukey[[i]] <- target_res_anova %>% 
      ggplot() +
      geom_segment(aes(x = conf.low, xend = conf.high, 
                       y = comb, yend = comb), colour = "grey", linewidth = 0.5) +
      geom_point(aes(x = estimate, y = comb), colour = "grey", size = 2) + 
      labs(subtitle = str_to_title(biol), x = "Estimate", y = "") +
      geom_vline(xintercept = 0, linetype="dotted", 
                 color = "black", linewidth = 0.5) +
      theme_ipsum(base_size = 10,
                  axis_title_size = 10,
                  strip_text_size = 10,
                  axis_text_size = 8,
                  axis_title_just = "center",
                  base_family = "Helvetica",
                  plot_margin = margin(5, 5, 5, 5)) +
      theme(legend.position = "right",
            legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.key.size = unit(5, 'mm'),
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  } else {
    vis_tukey[[i]] <- target_res_anova %>% 
      ggplot() +
      geom_segment(aes(x = conf.low, xend = conf.high, 
                       y = comb, yend = comb, colour = signif), linewidth = 0.5) +
      geom_point(aes(x = estimate, y = comb, colour = signif), size = 2) + 
      scale_colour_manual(values = c(pal_orig[i], "grey")) +
      labs(subtitle = str_to_title(biol), x = "Estimate", y = "") +
      geom_vline(xintercept = 0, linetype="dotted", 
                 color = "black", linewidth = 0.5) +
      theme_ipsum(base_size = 10,
                  axis_title_size = 10,
                  strip_text_size = 10,
                  axis_text_size = 8,
                  axis_title_just = "center",
                  base_family = "Helvetica",
                  plot_margin = margin(5, 5, 5, 5)) +
      theme(legend.position = "right",
            legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.key.size = unit(5, 'mm'),
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  }
}

vis_tukey[[1]] + vis_tukey[[2]] + vis_tukey[[3]] +
  vis_tukey[[4]] + vis_tukey[[5]] + vis_tukey[[6]] + vis_tukey[[7]] +
  plot_layout(ncol = 2, guides = "collect", widths = c(1, 1)) & 
  theme(legend.position = 'none') # 縦横比を設定し凡例をまとめ

# Save the visualized result
ggsave("fig-supp/forestplot_turkeyHSD-significance.png",
       units = "mm", width = 170, height = 200)
ggsave("fig-supp/forestplot_turkeyHSD-significance.eps",
       units = "mm", width = 170, height = 200, device = cairo_ps)

# Topicの類似度(groupからの距離) -----------------------------------------------------

# List of biological groups
l_biol <- c("mammal", "bird", "reptile", "amphibian", 
            "fish", "invertebrate", "plant")

# Prepare empty list
vis_bray <- list()

# Loop
for (i in 1:length(l_biol)) {
  biol <- l_biol[i]
  vis_bray[[i]] <- df_bray_pairwise %>%
    filter(group_biol == biol) %>% 
    ggplot(aes(x = reorder(var1, desc(value)), y = value)) +
    geom_boxplot(aes(fill = group_biol)) +
    labs(x = "",
         y = "Bray-Curtis dissimilarity",
         subtitle = str_to_title(biol)) +
    scale_fill_manual(values = pal_orig[i]) + # cols25かalphabet2のどちらかが良さそう。
    theme_ipsum(base_size = 10,
                axis_title_size = 10,
                strip_text_size = 10,
                axis_text_size = 8,
                axis_title_just = "center",
                base_family = "Helvetica",
                plot_margin = margin(5, 5, 5, 5)) + 
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "right",
          legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          legend.key.size = unit(5, 'mm'),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
}

# Layout and save the visualized result
((vis_bray[[1]] | vis_bray[[2]]) /
    (vis_bray[[3]] | vis_bray[[4]])) & theme(legend.position = 'none')

ggsave("fig-supp/boxplot-species-dissimilarity_A.png",
       units = "mm", width = 170, height = 200)
ggsave("fig-supp/boxplot-species-dissimilarity_A.eps",
       units = "mm", width = 170, height = 200, device = cairo_ps)

((vis_bray[[5]] | vis_bray[[6]]) / 
    vis_bray[[7]]) & theme(legend.position = 'none')

ggsave("fig-supp/boxplot-species-dissimilarity_B.png",
       units = "mm", width = 170, height = 200)
ggsave("fig-supp/boxplot-species-dissimilarity_B.eps",
       units = "mm", width = 170, height = 200, device = cairo_ps)

