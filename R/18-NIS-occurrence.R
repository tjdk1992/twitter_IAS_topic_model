#-----------------------------------------------------------------------------#
# Script Name: 18-NIS-occurrence.R                                            #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script calculates the occurrence of each IAS in the tweets.          #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Packages
pacman::p_load(tidyverse,
               readxl,
               lubridate,
               hrbrthemes,
               patchwork,
               glue,
               ggtext
)

# Data
tweet_finalized <- read_csv("data/tweet-finalized.csv")
NIS_all <- read_csv("data/NIS-compiled.csv")

# Color palette
pal_orig <- pals::cols25(7)[c(5, 4, 7, 2, 1, 6, 3)]

# Levels of taxonomic group
taxon_arrange <- c("Mammal", "Bird", "Reptile", "Amphibian", "Fish", "Invertebrate", "Plant")

# Temporal trends ---------------------------------------------------------------

# Check the annual trends of the occurrence
tweet_finalized %>% 
  group_by(name_sp, year) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  left_join(NIS_all %>% 
              dplyr::select(name_sp, group_biol) %>% 
              distinct(),
            by = "name_sp") %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange)) %>% 
  ggplot(aes(x = year, y = count)) +
  geom_line(aes(group = name_sp, colour = group_biol), 
            linewidth = 0.5,
            show.legend = FALSE) +
  facet_wrap(. ~ group_biol, ncol = 2) +
  scale_color_manual(values = pal_orig) +
  theme_ipsum(base_size = 10,
              axis_title_size = 10,
              strip_text_size = 10,
              axis_text_size = 8,
              axis_title_just = "center",
              base_family = "Helvetica")

# Extract typical value -------------------------------------------------------

# 最初に出現した年の分布
tweet_finalized %>%
  group_by(name_sp, year) %>% 
  summarise(count = n()) %>% 
  group_by(name_sp) %>% 
  summarise(year_1st_occur = min(year)) %>% 
  group_by(year_1st_occur) %>% 
  summarise(n = n())

# Calculate typical values
tweet_count <- tweet_finalized %>%
  group_by(name_sp, year) %>% 
  summarise(count = n()) %>% 
  group_by(name_sp) %>% 
  summarise(total = sum(count),
            mean = mean(count),
            sd = sd(count),
            se = sd(count) / n(),
            mean_1 = mean(count, trim = 0.1),
            mean_3 = mean(count, trim = 0.3),
            mean_5 = mean(count, trim = 0.5),
            median = median(count),
            min = min(count),
            max = max(count),
            mid_range = ((min + max)/2))

# カタカナ名による重複を除外しておく
glimpse(NIS_all)

NIS_all <- NIS_all %>% 
  dplyr::select(group_biol, name_ja, name_sp, reg1, reg2) %>% 
  distinct(.keep_all = TRUE)

# Merge
tweet_count <- left_join(tweet_count, NIS_all, by = "name_sp")

# Write data
write_csv(tweet_count, "data/NIS-count.csv")

# Write table
NIS_zero <- NIS_all %>% 
  dplyr::select(name_ja, name_sp, group_biol) %>% 
  anti_join(dplyr::select(tweet_count, name_sp), by = "name_sp")

tweet_count %>% 
  bind_rows(NIS_zero) %>% 
  mutate(total = replace_na(total, 0)) %>% 
  arrange(name_sp) %>% 
  arrange(desc(total)) %>% 
  arrange(group_biol) %>% 
  transmute(
    `No.` = row_number(),
    `Scientific name` = name_sp,
    `Japanese name`   = str_replace_all(name_ja, c("[[:punct:]]" = "",
                                                   "[[:digit:]]" = "",
                                                   "[[:lower:]]" = "",
                                                   "[[:upper:]]" = "",
                                                   "[[:space:]]" = "")),
    Taxon             = str_to_title(group_biol),
    Total             = round(total, 2)
    ) %>% 
  writexl::write_xlsx("table-supp/table-NIS-summary-unformatted.xlsx")

# Summarizing unfiltered values ------------------------------------------------------

# All
tweet_count %>% 
  summarise(mean = round(mean(total), 2),
            sd = round(sd(total), 2),
            min = min(total),
            max = max(total)) %>% 
  as.data.frame()

# Max and Minimum
tweet_count %>% 
  arrange(desc(total)) %>% 
  head(5)
tweet_count %>% 
  filter(total == 1) %>% 
  group_by(group_biol) %>% 
  summarise(n = n())

# Distribution
tweet_count %>% 
  ggplot(aes(x = total)) +
  geom_histogram(bins = 50) +
  labs(x = "The number of tweets", 
       y = "The number of species") +
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
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
  labs(x = "")

# Save the visualized result
ggsave("fig-supp/histogram-ias-count.png",
       units = "mm", width = 150, height = 100)
ggsave("fig-supp/histogram-ias-count.eps",
       units = "mm", width = 150, height = 100, device = cairo_ps)

# Popular NIS
tweet_count %>% 
  arrange(desc(total)) %>% 
  head(10) %>% 
  summarise(sum = sum(total))
tweet_count %>% 
  summarise(sum = sum(total))
round(88844/186028 * 100, 2)
tweet_count %>% 
  filter(total >= 100) %>% 
  group_by(group_biol) %>% 
  summarise(n = n())

# Top NIS
arrange(filter(tweet_count, group_biol == "mammal"), desc(total))
arrange(filter(tweet_count, group_biol == "bird"), desc(total))
arrange(filter(tweet_count, group_biol == "reptile"), desc(total))
arrange(filter(tweet_count, group_biol == "amphibian"), desc(total))
arrange(filter(tweet_count, group_biol == "fish"), desc(total))
arrange(filter(tweet_count, group_biol == "invertebrate"), desc(total))
arrange(filter(tweet_count, group_biol == "plant"), desc(total))

# Taxonomic group for popular NIS
NIS_popular <- tweet_count %>% 
  filter(total >= 100) %>% 
  dplyr::select(name_sp, name_ja, group_biol, total)

NIS_popular %>% 
  group_by(group_biol) %>% 
  summarise(mean = round(mean(total), 2),
            sd = round(sd(total), 2),
            min = min(total),
            max = max(total)) %>% 
  as.data.frame()

write_csv(NIS_popular, "data/NIS-popular.csv")

# Visualization ---------------------------------------------------------------

# Plot bargraph for all species
NIS_popular %>% 
  arrange(desc(total)) %>% 
  head(50) %>% 
  mutate(group_biol = str_to_title(group_biol),
         group_biol = factor(group_biol, levels = taxon_arrange),
         name_italic = str_remove_all(name_sp, c(" subspp." = "", " spp." = "")),
         name_block = if_else(str_detect(name_sp, "subspp."), "subsp.",
                              if_else(str_detect(name_sp, "spp."), "spp.", "")),
         name_block = str_replace_all(name_block, "subsp.", "subspp."),
         name_show = glue("<i>{name_italic}</i> {name_block}")) %>% 
  ggplot(aes(x = reorder(name_show, total), y = total)) +
  geom_bar(aes(fill = group_biol), stat = "identity") +
  scale_fill_manual(values = pal_orig, name = "taxonomic group") + # cols25かalphabet2のどちらかが良さそう。
  coord_flip() +
  labs(x = "Species", 
       y = "Frequency of NIS name") +
  theme_ipsum(axis_title_just = "center",
              base_family = "Helvetica",
              base_size = 8,
              plot_margin = margin(5, 5, 5, 5)) + 
  theme(axis.text.y = element_markdown(hjust = 1, size = 8),
        axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right",
        legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, 'mm'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
  labs(x = "")

# Save the visualized result
ggsave("fig/barplot_top-occurred-ias.png",
       units = "mm", width = 120, height = 160)
ggsave("fig/barplot_top-occurred-ias.eps",
       units = "mm", width = 120, height = 160, device = cairo_ps)

# Plot bargraph of occurrence separated by biological groups
l_biol <- c("mammal", "bird", "reptile", "amphibian", 
            "fish", "invertebrate", "plant")
vis_n_tweet <- list()
for (i in 1:length(l_biol)) {
  # Biological group
  biol <- l_biol[i]
  vis_n_tweet[[i]] <- tweet_count %>% 
    filter(total >= 100 & group_biol == biol) %>% 
    mutate(group_biol = str_to_title(group_biol),
           name_italic = str_remove_all(name_sp, c(" subspp." = "", " spp." = "")),
           name_block = if_else(str_detect(name_sp, "subspp."), "subsp.",
                                if_else(str_detect(name_sp, "spp."), "spp.", "")),
           name_block = str_replace_all(name_block, "subsp.", "subspp."),
           name_show = glue("<i>{name_italic}</i> {name_block}")) %>% 
    ggplot(aes(x = reorder(name_show, desc(total)), y = total)) +
    geom_bar(aes(fill = group_biol), stat = "identity") +
    scale_fill_manual(values = pal_orig[i]) + # cols25かalphabet2のどちらかが良さそう。
    labs(x = "", 
         y = "No. of tweets",
         subtitle = str_to_title(biol)) +
    theme_ipsum(base_size = 10,
                axis_title_size = 10,
                strip_text_size = 10,
                axis_text_size = 8,
                axis_title_just = "center",
                base_family = "Helvetica",
                plot_margin = margin(5, 5, 5, 5)) + 
    theme(axis.text.x = element_markdown(angle = 60, hjust = 1),
          legend.position = "right",
          legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          legend.key.size = unit(5, 'mm'),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
}

# Layout and save the visualized result
((vis_n_tweet[[1]] | vis_n_tweet[[2]]) /
    (vis_n_tweet[[3]] | vis_n_tweet[[4]])) & theme(legend.position = 'none')
ggsave("fig-supp/barplot_ias-tweet-count_A.png",
       units = "mm", width = 170, height = 200)
ggsave("fig-supp/barplot_ias-tweet-count_A.eps",
       units = "mm", width = 170, height = 200, device = cairo_ps)

((vis_n_tweet[[5]] | vis_n_tweet[[6]]) / 
    vis_n_tweet[[7]]) & theme(legend.position = 'none')
ggsave("fig-supp/barplot_ias-tweet-count_B.png",
       units = "mm", width = 170, height = 200)
ggsave("fig-supp/barplot_ias-tweet-count_B.eps",
       units = "mm", width = 170, height = 200, device = cairo_ps)
