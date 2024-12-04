#-----------------------------------------------------------------------------#
# Script Name: 13-LDA-tuning-automatic.R                                      #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script explores optimal number of topics by running the ldatuning(). #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse,
               magrittr,
               tidytext,
               topicmodels,
               ldatuning,
               pals,
               hrbrthemes,
               tictoc
               )

# Data
token_finalized <- read_csv("data/token-finalized.csv")

# Run ldatuning ---------------------------------------------------------------

# DTMの作成
dtm_rm_stpw_original <- token_finalized %>% 
  anti_join(token_finalized %>% 
              group_by(id_cleansed) %>% 
              summarise(n = n()) %>% 
              filter(n < 5) %>% 
              dplyr::select(id_cleansed),
            by = "id_cleansed") %>% 
  group_by(id_cleansed, term) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  tidytext::cast_dtm(document = "id_cleansed",
                     term = "term",
                     value = "count")

# Run  
res_ldatuning <- FindTopicsNumber(
  dtm_rm_stpw_original,
  topics = c(1:10 * 5, 6:11 * 10,  110, 120, 130, 140, 150, 0:3 * 50 + 150),
  metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 123),
  mc.cores = NA,
  verbose = TRUE
)

# Resultの書き出し
write_csv(as.data.frame(res_ldatuning),
          "data/LDA-tuning-automatic.csv")

# Visualize the results -------------------------------------------------------

# Read the results 
res_ldatuning <- read_csv("data/LDA-tuning-automatic.csv")

# Visualize the results
res_ldatuning %>% 
    pivot_longer(-topics, names_to = "metrics", values_to = "value") %>% 
    group_by(metrics) %>% 
    mutate(value_scaled = scale(value)[,1]) %>% 
    ungroup() %>% 
    mutate(group = if_else(metrics == "Arun2010" | 
                             metrics == "CaoJuan2009", "1", "2")) %>% 
    ggplot(aes(x = topics, y = value_scaled, group = metrics, colour = metrics)) + 
    geom_point(aes(shape = metrics)) +
    geom_line(aes(linetype = metrics)) +
    annotate("rect", 
             xmin = 10, xmax = 60, 
             ymin = -5, ymax = 5,
             alpha = 0.2, fill = "grey") +
    scale_x_continuous(breaks = seq(0, 300, 10)) +
    facet_grid(group ~ .) +
    scale_color_manual(values=as.vector(cols25(4))) +
    labs(x = "The number of topics (K)", y = "Metric score") +
  theme_ipsum(base_size = 10,
              axis_title_size = 10,
              strip_text_size = 10,
              axis_text_size = 8,
              axis_title_just = "center",
              base_family = "Helvetica") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_blank())

# Save the visualized result
ggsave("fig-supp/ldatuning-output.png",
       units = "mm", width = 180, height = 180)
ggsave("fig-supp/ldatuning-output.eps",
       units = "mm", width = 180, height = 180, device = cairo_ps)
