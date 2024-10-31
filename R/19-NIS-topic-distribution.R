#-----------------------------------------------------------------------------#
# Script Name: 19-NIS-topic-distribution.R                                  #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script analyzes topic distribution over taxonomic groups and         #
# predictability of the IAS occurrence by the topic distribution              #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Packages
pacman::p_load(tidyverse,  # for data manipulation
               hrbrthemes, # for nice visualization
               magrittr, 　# for data manipulation
               patchwork,  # to use color palette
               reshape2,
               rstatix,
               vegan,
               glue,
               ggtext
               )

# Data
# LDA output
theta <- read_csv("data/LDA-doc-topic-tweet.csv")

# Tweets counts
ias_popular <- read_csv("data/NIS-popular.csv")

# Color palette
pal_orig <- pals::cols25(7)[c(5, 4, 7, 2, 1, 6, 3)]

# Levels of taxonomic group
taxon_arrange <- c("Mammal", "Bird", "Reptile", "Amphibian", "Fish", "Invertebrate", "Plant")

# Data preparation ------------------------------------------------------------

# Count of document aligned given topics
N_doc <- theta %>% 
  group_by(name_sp) %>% 
  summarise(n_doc = n())
theta_topic <- theta %>% 
  group_by(name_sp, topic) %>% 
  summarise(n = n()) %>% 
  left_join(N_doc, by = "name_sp") %>% 
  mutate(freq = n / n_doc) %>% 
  dplyr::select(name_sp, topic, freq) %>% 
  ungroup() %>% 
  arrange(topic) %>% 
  pivot_wider(names_from = topic, values_from = freq)
theta_topic[is.na(theta_topic)] <- 0

# Merge count data to LDA output
NIS_topic <- inner_join(ias_popular, theta_topic, by = "name_sp")

# Export the association data
write_csv(NIS_topic, "data/NIS-popular-topic.csv")

# Topic distribution over biological groups -----------------------------------

NIS_topic %>% 
  pivot_longer(cols = TP01:TP25,
               names_to = "topic",
               values_to = "freq") %>% 
  group_by(topic, group_biol) %>% 
  summarise(prob = mean(freq)) %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange)) %>% 
  ggplot(aes(x = topic, y = prob)) +
  geom_bar(aes(fill = group_biol), stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = pal_orig, name = "taxonomic group") + 
  facet_wrap(. ~ group_biol, ncol = 2) +
  labs(x = "", y = "") +
  theme_ipsum(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_markdown(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.spacing=unit(0.5, "lines"))

# Save the visualized result
ggsave("fig-supp/bar-topic-distrib-taxon.png",
       units = "mm", width = 150, height = 130)
ggsave("fig-supp/bar-topic-distrib-taxon.eps",
       units = "mm", width = 150, height = 130, device = cairo_ps)

# Topic distribution over NISs ------------------------------------------------

# Bubble plot
g_bubble_topic <- NIS_topic %>% 
  arrange(desc(total)) %>% 
  group_by(group_biol) %>% 
  mutate(rank_group = row_number()) %>% 
  ungroup() %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         name_italic = str_remove_all(name_sp, c(" subspp." = "", " spp." = "")),
         name_block = if_else(str_detect(name_sp, "subspp."), "subsp.",
                              if_else(str_detect(name_sp, "spp."), "spp.", "")),
         name_block = str_replace_all(name_block, "subsp.", "subspp."),
         name_show = glue("<i>{name_italic}</i> {name_block}")) %>% 
  filter(group_biol == "Mammal" | 
           group_biol == "Bird" |
           group_biol == "Reptile"|
           group_biol == "Amphibian"|
           group_biol == "Fish") %>% 
  arrange(total) %>% 
  arrange(desc(group_biol)) %>% 
  mutate(id_reorder = row_number()) %>% 
  pivot_longer(cols = TP01:TP25, 
               names_to = "topic", 
               values_to = "value") %>% 
  ggplot(aes(x = topic, y = reorder(name_show, id_reorder), label = group_biol)) +
  geom_point(aes(size = value, color = group_biol), alpha = 0.7) +
  scale_color_manual(values = pal_orig[1:5], name = "taxonomic group") + # cols25かalphabet2のどちらかが良さそう。
  scale_size(range = c(0.05, 10)) +  # Adjust the range of points size
  scale_x_discrete(position = "top") +
  labs(x = "", y = "") +
  theme_ipsum(base_family = "Helvetica",
              base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_markdown(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

g_bar_rank <- NIS_topic %>% 
  arrange(desc(total)) %>% 
  group_by(group_biol) %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         rank_group = row_number()) %>% 
  ungroup() %>% 
  filter(group_biol == "Mammal" | 
           group_biol == "Bird" |
           group_biol == "Reptile"|
           group_biol == "Amphibian"|
           group_biol == "Fish") %>% 
  arrange(total) %>% 
  arrange(desc(group_biol)) %>% 
  mutate(id_reorder = row_number()) %>% 
  ggplot() + 
  geom_bar(aes(x = reorder(name_sp, id_reorder), 
               y = total, 
               fill = group_biol), 
           stat = "identity",
           show.legend = FALSE) +
  geom_text(aes(x = reorder(name_sp, id_reorder), 
                y = total, 
                label = total, 
                hjust = -0.2), size = 3) +
  scale_fill_manual(values = pal_orig[1:5]) + # cols25かalphabet2のどちらかが良さそう。
  labs(x = "Species", 
       y = "No. of tweets") +
  theme_ipsum(base_family = "Helvetica", 
              base_size = 8, 
              axis_text_size = 8,
              axis_title_size = 10,
              axis_title_just = "mc") + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60),
        axis.title.y = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
  coord_flip()

g_bubble_topic + 
  g_bar_rank + 
  plot_layout(guides = "collect", widths = c(7, 2)) & 
  theme(legend.position = 'bottom') # 縦横比を設定し凡例をまとめ

# Save the visualized result
ggsave("fig-supp/bubble-biol-ordered-category_A.png",
       units = "mm", width = 170, height = 230)
ggsave("fig-supp/bubble-biol-ordered-category_A.eps",
       units = "mm", width = 170, height = 230, device = cairo_ps)

g_bubble_topic <- NIS_topic %>% 
  arrange(desc(total)) %>% 
  group_by(group_biol) %>% 
  mutate(rank_group = row_number()) %>% 
  ungroup() %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         name_italic = str_remove_all(name_sp, c(" subspp." = "", " spp." = "")),
         name_block = if_else(str_detect(name_sp, "subspp."), "subsp.",
                              if_else(str_detect(name_sp, "spp."), "spp.", "")),
         name_block = str_replace_all(name_block, "subsp.", "subspp."),
         name_show = glue("<i>{name_italic}</i> {name_block}")) %>% 
  filter(group_biol == "Invertebrate" | 
           group_biol == "Plant") %>% 
  arrange(total) %>% 
  arrange(desc(group_biol)) %>% 
  mutate(id_reorder = row_number()) %>% 
  pivot_longer(cols = TP01:TP25, 
               names_to = "topic", 
               values_to = "value") %>% 
  ggplot(aes(x = topic, y = reorder(name_show, id_reorder), label = group_biol)) +
  geom_point(aes(size = value, color = group_biol), alpha = 0.7) +
  scale_color_manual(values = pal_orig[6:7], name = "taxonomic group") + # cols25かalphabet2のどちらかが良さそう。
  scale_size(range = c(0.05, 10)) +  # Adjust the range of points size
  scale_x_discrete(position = "top") +
  labs(x = "", y = "") +
  theme_ipsum(base_family = "Helvetica",
              base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_markdown(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

g_bar_rank <- NIS_topic %>% 
  arrange(desc(total)) %>% 
  group_by(group_biol) %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         rank_group = row_number()) %>% 
  ungroup() %>% 
  filter(group_biol == "Invertebrate" | 
           group_biol == "Plant") %>% 
  arrange(total) %>% 
  arrange(desc(group_biol)) %>% 
  mutate(id_reorder = row_number()) %>% 
  ggplot() + 
  geom_bar(aes(x = reorder(name_sp, id_reorder), 
               y = total, 
               fill = group_biol), 
           stat = "identity",
           show.legend = FALSE) +
  geom_text(aes(x = reorder(name_sp, id_reorder), 
                y = total, 
                label = total, 
                hjust = -0.2), size = 3) +
  scale_fill_manual(values = pal_orig[6:7]) + # cols25かalphabet2のどちらかが良さそう。
  labs(x = "Species", 
       y = "NIS name frequency") +
  theme_ipsum(base_family = "Helvetica", 
              base_size = 8, 
              axis_text_size = 8,
              axis_title_size = 10,
              axis_title_just = "mc") + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60),
        axis.title.y = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
  coord_flip()

g_bubble_topic + 
  g_bar_rank + 
  plot_layout(guides = "collect", widths = c(7, 2)) & 
  theme(legend.position = 'bottom') # 縦横比を設定し凡例をまとめ

# Save the visualized result
ggsave("fig-supp/bubble-biol-ordered-category_B.png",
       units = "mm", width = 170, height = 230)
ggsave("fig-supp/bubble-biol-ordered-category_B.eps",
       units = "mm", width = 170, height = 230, device = cairo_ps)

# ひとまとめ ------------------------------------------------------------------

NIS_topic %>% 
  arrange(desc(total)) %>% 
  group_by(group_biol) %>% 
  mutate(rank_group = row_number()) %>% 
  ungroup() %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         name_italic = str_remove_all(name_sp, c(" subspp." = "", " spp." = "")),
         name_block = if_else(str_detect(name_sp, "subspp."), "subsp.",
                              if_else(str_detect(name_sp, "spp."), "spp.", "")),
         name_block = str_replace_all(name_block, "subsp.", "subspp."),
         name_show = glue("<i>{name_italic}</i> {name_block}")) %>% 
  filter(rank_group <= 5) %>% 
  arrange(total) %>% 
  arrange(desc(group_biol)) %>% 
  mutate(id_reorder = row_number()) %>% 
  pivot_longer(cols = TP01:TP25, 
               names_to = "topic", 
               values_to = "Frequency") %>% 
  ggplot(aes(x = topic, y = reorder(name_show, id_reorder), label = group_biol)) +
  geom_point(aes(size = Frequency, color = group_biol), alpha = 0.7) +
  scale_color_manual(values = pal_orig, name = "taxonomic group") + # cols25かalphabet2のどちらかが良さそう。
  scale_size(range = c(0.05, 10)) +  # Adjust the range of points size
  scale_x_discrete(position = "top") +
  labs(x = "", y = "") +
  theme_ipsum(base_family = "Helvetica",
              base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
        axis.text.y = element_markdown(size = 8),
        axis.title = element_text(size = 8),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Save the visualized result
ggsave("fig/bubble-biol-ordered-category.png",
       units = "mm", width = 150, height = 170)
ggsave("fig/bubble-biol-ordered-category.eps",
       units = "mm", width = 150, height = 170, device = cairo_ps)
