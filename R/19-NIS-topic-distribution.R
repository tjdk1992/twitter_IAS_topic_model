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

# Summarise the distribution by taxonomic groups
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
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.spacing=unit(0.5, "lines"))

# Save the visualized result
ggsave("fig-supp/bar-topic-distrib-taxon.png",
       units = "mm", width = 150, height = 130)
ggsave("fig-supp/bar-topic-distrib-taxon.eps",
       units = "mm", width = 150, height = 130, device = cairo_ps)

# Summarise the distribution by taxonomic groups and thematic groups
NIS_theme <- NIS_topic %>% 
  pivot_longer(cols = starts_with("TP"),
               names_to = "topic",
               values_to = "freq") %>% 
  mutate(theme = case_when(str_detect(topic, "TP01|TP02|TP10|TP17|TP20|TP24") ~ "theme01",
                           str_detect(topic, "TP03|TP04|TP06|TP15") ~ "theme02",
                           str_detect(topic, "TP05|TP11|TP14|TP19|TP21|TP23|TP25") ~ "theme03",
                           str_detect(topic, "TP07|TP09|TP12|TP16") ~ "theme04",
                           str_detect(topic, "TP08|TP13|TP18|TP22") ~ "theme05",
                           TRUE ~ "undecided"),
         theme_name = str_replace_all(theme,
                                      c("theme01" = "Biol. dime. NIS", # Biological dimension of NIS
                                        "theme02" = "People recogn. NIS", # People's recognition of NIS
                                        "theme03" = "NIS manag.", # NIS management
                                        "theme04" = "Ecol. prob. NIS", # Ecological problem related to NIS
                                        "theme05" = "Inter. human. NIS" # Interaction between human and NIS
                                        )))

# Check
NIS_theme %>% 
  dplyr::select(topic, theme) %>% 
  distinct() %>% 
  table()

NIS_topic %>% 
  pivot_longer(cols = TP01:TP25,
               names_to = "topic",
               values_to = "freq") %>% 
  group_by(topic, group_biol) %>% 
  summarise(prob = mean(freq)) %>% 
  group_by(group_biol) %>% 
  summarise(prob = sum(prob))

# Visualise
n_species_group <- NIS_theme %>% 
  dplyr::select(name_sp, group_biol) %>% 
  distinct() %>% 
  group_by(group_biol) %>% 
  summarise(n_species = n())

# Stacked + percent
NIS_theme %>% 
  group_by(theme_name, group_biol) %>% 
  summarise(prob_sum = sum(freq)) %>% 
  left_join(n_species_group, by = "group_biol") %>% 
  mutate(prob_theme = prob_sum / n_species) %>% 
  # group_by(group_biol) %>% # Check
  # summarise(sum = sum(prob_theme)) # Check
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange)) %>% 
  ggplot(aes(fill = theme_name, y = prob_theme, x = group_biol)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = pal_orig) +
  theme_ipsum(base_size = 8) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
        strip.text = element_markdown(hjust = 0.5),
        strip.placement = "outside",
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.spacing=unit(0.5, "lines"))

"Biol. dime. NIS"
"People recogn. NIS"
"NIS manag."
"Ecol. prob. NIS"
"Inter. human. NIS"
# Bargraph
NIS_theme %>% 
  group_by(theme_name, group_biol) %>% 
  summarise(prob_sum = sum(freq)) %>% 
  left_join(n_species_group, by = "group_biol") %>% 
  mutate(prob_theme = prob_sum / n_species) %>% 
  # group_by(group_biol) %>% # Check
  # summarise(sum = sum(prob_theme)) # Check
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         theme_order = str_replace_all(theme_name, c("Biol. dime. NIS" = "1",
                                                     "People recogn. NIS" = "2",
                                                     "NIS manag." = "3",
                                                     "Ecol. prob. NIS" = "4",
                                                     "Inter. human. NIS"= "5")),
         theme_order = as.numeric(theme_order)) %>% 
  ggplot(aes(x = reorder(theme_name, theme_order, decreasing = TRUE), y = prob_theme)) +
  geom_bar(aes(fill = group_biol), stat = "identity", show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = pal_orig, name = "taxonomic group") + 
  #facet_grid(group_biol ~ ., switch = "both") +
  facet_wrap(group_biol ~ ., ncol = 1, strip.position = "left") +
  labs(x = "", y = "") +
  theme_ipsum(base_size = 8) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
        strip.text = element_markdown(hjust = 0.5, size = 10),
        strip.placement = "outside",
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.spacing=unit(0.5, "lines"))

# Save the visualized result
ggsave("fig/bar-theme-distrib-taxon.png",
       units = "mm", width = 74, height = 135)
ggsave("submission/manuscript-source/fig/figure4.eps",
       units = "mm", width = 74, height = 135, device = cairo_ps)

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

# Add topic labels
topic_label <- data.frame(
  topic = c("TP01", "TP02", "TP03", "TP04", "TP05", "TP06", "TP07", "TP08", 
            "TP09", "TP10", "TP11", "TP12", "TP13", "TP14", "TP15", "TP16", 
            "TP17", "TP18", "TP19", "TP20", "TP21", "TP22", "TP23", "TP24", "TP25"),
  topic_label = c("Occurrence of NIS", "Population dynamics", "Perception of NIS", 
                  "Temporal trend", "Regulation", "Detection of poisonous NIS", 
                  "Ecological damage", "Release of NIS", "Extinction of endemism", 
                  "Establishment", "Draining management", "Environmental issue", 
                  "Consumption of NIS", "Capture for conservation", 
                  "Danger information", "Hybridization and introgression", 
                  "Reproductive capacity", "Emotion", "Handling NIS", 
                  "Flower as color traits", "Control measures", "Keep as pets", 
                  "Human dimension of management", "Origin of NIS", "Destruction of NIS")) %>% 
  mutate(topic_label = str_c(topic_label, " (", topic, ")"))

# Visualization
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
  left_join(topic_label, by = "topic") %>% 
  mutate(topic_no = str_replace_all(topic, "TP", ""),
         topic_no = as.numeric(topic_no)) %>% 
  ggplot(aes(x = reorder(topic_label, topic_no), y = reorder(name_show, id_reorder), label = group_biol)) +
  geom_point(aes(size = Frequency, color = group_biol), alpha = 0.7) +
  scale_color_manual(values = pal_orig, name = "taxonomic group") + # cols25かalphabet2のどちらかが良さそう。
  scale_size(range = c(0.05, 10)) +  # Adjust the range of points size
  # scale_x_discrete(position = "top") +
  labs(x = "", y = "") +
  theme_ipsum(base_family = "Helvetica",
              base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_markdown(size = 8),
        axis.title = element_text(size = 8),
        # legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_text(size = 8),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))

# Save the visualized result
ggsave("fig/bubble-biol-ordered-category.png",
       units = "mm", width = 174, height = 200)
ggsave("fig/bubble-biol-ordered-category.eps",
       units = "mm", width = 174, height = 200, device = cairo_ps)
