#-----------------------------------------------------------------------------#
# Script Name: 17-LDA-topic-labelling.R                                       #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script labels topics by checking term frequency in each topic.       #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Packages
pacman::p_load(tidyverse,
               magrittr,
               wordcloud,
               pals,
               readxl,
               writexl)

# Color palette
pal_orig <- c(rep(pals::cols25(25), 2))

# Data
topic_term <- read_excel("data/LDA-topic-top-term.xlsx")

# Labelling topic -------------------------------------------------------------

# View the top terms
topic_term[, 1:10]
topic_term[, 11:20]
topic_term[, 21:25]

label_topic <- data.frame(
  topic = str_c("Topic ", seq(25)),
  label_topic = c("Occurrence of NIS", "Population dynamics", # 1, 2
                  "Perception of NIS", "Temporal trend", # 3, 4
                  "Regulation", "Detection of poisonous NIS", # 5, 6
                  "Ecological damage", "Release of NIS", # 7, 8
                  "Extinction of endemic things", "Establishment", # 9, 10
                  "Draining management", "Environmental issue", # 11, 12
                  "Food use", "Capture for conservation", # 13, 14
                  "Danger information", "Hybridization and introgression", # 15, 16
                  "Reproductive capacity", "Emotion", # 17, 18
                  "Handling Invasive NIS", "Flower as color traits", # 19, 20
                  "Control measures", "Keep as pets", # 21, 22
                  "Human dimension of management", "Origin of NIS", # 23, 24
                  "Destroy NIS" # 25
                  ))

# Export table ----------------------------------------------------------------

# Top 20 terms in Japanese for Supplementary materials
(topic_term_top20 <- head(topic_term, 20))

# Write as tex format
writexl::write_xlsx(topic_term_top20, "table-supp/table-top20-unformatted.xlsx")

# Write table for English translation
topic_term %>% 
  head(5) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(topic = str_c("Topic ", seq(25))) %>% 
  dplyr::select(topic, V1:V5) %>% 
  write_xlsx("table/table-top5-term-Ja.xlsx")

# Read translated table
topic_term_en <- read_excel("table/table-top5-term-En.xlsx")

# Format final table for manuscript
table_topic_en <- topic_term_en %>% 
  mutate(terms = str_c(V1, ", ", V2, ", ", V3, ", ", V4, ", ", V5)) %>%
  dplyr::select(topic, terms) %>% 
  left_join(label_topic, by = "topic") %>% 
  transmute(Topic = topic,
            `Topic label` = label_topic,
            `Top 10 terms` = terms)

# Write table as tex format
print(xtable(table_topic_en), include.rownames = FALSE,
      file = "table/table-top5-term-En.tex")
