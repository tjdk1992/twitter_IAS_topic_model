#-----------------------------------------------------------------------------#
# Script Name: 14-LDA-tuning-automatic.R                                      #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script explores optimal number of topics by checking multiple        #
# patterns of combination of the topic number K and random seeds.             #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse, # for data manipulation
               tidytext, # to use cast_dtm function
               textmineR, # to create DTM
               ldatuning, # determine the proper number of topics inferred
               tictoc, # to calculate ran time
               hrbrthemes, # for visualization
               doParallel,
               scales,
               topicmodels
               )

# Color palette
pal_orig <- c(rep(pals::cols25(25), 2))

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

# Topic modelling with different K --------------------------------------------

# ldatuningの結果を踏まえて10〜60 BY 10でトピックを推定してみる。
for (K in seq(10, 60, by = 5)) {
  topicModel <- LDA(dtm_rm_stpw_original,
                    k = K,
                    method = "Gibbs",
                    control = list(alpha = 50/K, 
                                   iter = 1000, 
                                   verbose = 25, 
                                   seed = 123))
  name_path <- str_c("data-manual/LDA-tuning-manual-A/LDA-tuning-manual-K", 
                     K,
                     ".xlsx")
  as.data.frame(terms(topicModel, 20)) %>% 
    writexl::write_xlsx(name_path)
}

#------------------------------------------------------------------------------

# ldatuningの結果を踏まえて15〜35 BY 1でトピックを推定してみる。
for (K in seq(15, 35, 1)) {
  topicModel <- LDA(dtm_rm_stpw_original,
                    k = K,
                    method = "Gibbs",
                    control = list(alpha = 50/K, 
                                   iter = 1000, 
                                   verbose = 25, 
                                   seed = 123))
  name_path <- str_c("data-manual/LDA-tuning-manual-B/LDA-tuning-manual-K", 
                     K,
                     ".xlsx")
  as.data.frame(terms(topicModel, 20)) %>% 
    writexl::write_xlsx(name_path)
} # K = 25で決定

#------------------------------------------------------------------------------

K = 25

for (rn in c(123, 135, 159, 246, 369)) {
    topicModel <- LDA(dtm_rm_stpw_original,
                      k = K,
                      method = "Gibbs",
                      control = list(alpha = 50/K, 
                                     iter = 1000, 
                                     verbose = 25, 
                                     seed = rn))
    name_path <- str_c("data-manual/LDA-tuning-manual-C/LDA-tuning-manual-K25_seed", 
                       rn,
                       ".xlsx")
    as.data.frame(terms(topicModel, 20)) %>% 
      writexl::write_xlsx(name_path)
  }
