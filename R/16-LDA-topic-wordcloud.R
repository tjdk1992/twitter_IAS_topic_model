#-----------------------------------------------------------------------------#
# Script Name: 16-LDA-topic-overview.R                                        #
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
               hrbrthemes,
               pals,
               readxl,
               writexl)

# Color palette
pal_orig <- c(rep(pals::cols25(25), 2))

# Data
beta <- read_csv("data/LDA-topic-term.csv")
theta <- read_csv("data/LDA-doc-topic-tweet.csv")

theta_raw <- read_csv("data/LDA-doc-topic.csv")

theta_raw %>% 
  mutate(id = str_c("doc", row_number())) %>% 
  pivot_longer(cols = -id,
               names_to = "topic",
               values_to = "value") %>% 
  group_by(id) %>% 
  summarise(max = max(value)) %>% 
  group_by(id) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

test <- data.frame(x = c("", "", "", ""),
           y = c(0.1, 0.1, 0.04, 0.08))

test %>% 
  summarise(max = max(y))

# Topic distribution ----------------------------------------------------------

# Check the topic distribution
theta %>% 
  group_by(topic) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = reorder(topic, as.numeric(str_remove_all(topic, "TP"))), y = n)) +
  geom_bar(stat = "identity", fill = pal_orig[1]) +
  labs(x = "Topic", y = "The number of tweets") +
  theme_ipsum(base_size = 10,
              axis_title_size = 10,
              strip_text_size = 10,
              axis_text_size = 8,
              axis_title_just = "center",
              base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90))

# Topic labelling--------------------------------------------------------------

# The number of topics
K <- nrow(beta)

# Create data frame for Wordclouds
beta_df <- as.data.frame(beta)
beta_t_matrix <- t(beta_df)
beta_t_df <- as.data.frame(beta_t_matrix)

df_topic_name <- data.frame(topic = "topic", 
                            num1 = c(rep("0", 9), 
                                     rep("", K-9)),
                            num2 = seq(1:K)) %>% 
  transmute(topic = str_c(topic, num1, num2))

colnames(beta_t_df) <- df_topic_name$topic
beta_t_df$term <- rownames(beta_t_df)

# Create wordcloud of each topic
for (i in 1:K) {
  # term-frequency
  df_wc <- beta_t_df %>% 
    transmute(term, freq = beta_t_df[, i]) %>% 
    arrange(desc(freq)) %>% 
    head(30) # top 30 words
  # term
  words <- df_wc$term
  # frequency
  freq <- as.vector(unlist(df_wc$freq))
  # family to show Japanese words
  par(family = "HiraKakuProN-W3")
  # random seed
  set.seed(135)
  # plot wordcloud
  path_wc <- str_c("fig/wordcloud/wc-topic", i, ".png")
  png(path_wc, width = 600, height = 600)
  par(family="HiraKakuPro-W3") ##日本語フォントでもok
  wordcloud(words = words, 
            freq = freq, 
            max.words = 30, 
            random.order = FALSE, 
            rot.per = 0.50, 
            colors = cols25(6))
  dev.off()
  Sys.sleep(1)
}
