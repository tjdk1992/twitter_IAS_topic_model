#-----------------------------------------------------------------------------#
# Script Name: 12-LDA-stopword-examination.R                                  #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script is for test run checking validity of the original stop words. #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse,
               magrittr,
               readxl,
               tidytext,
               topicmodels)

# Data
token_rm_stpw_basic <- read_csv("data-proc/token-04-basic-stopword-removed.csv")
token_rm_stpw_general <- read_csv("data-proc/token-05-general-stopword-removed.csv")
token_rm_stpw_original <- read_csv("data-proc/token-06-original-stopword-removed.csv")

# Final decision --------------------------------------------------------------

# Re-write the selected corpus for LDA inference
write_csv(token_rm_stpw_original, "data/token-finalized.csv")

# Helper function -------------------------------------------------------------

# Function for creating doc-term matrix
makeDTMatrix <- function(df_token) {
  df_token %>% 
    anti_join(df_token %>% 
                group_by(id_cleansed) %>% 
                summarise(n = n()) %>% 
                # 単語が1個/1Tweetの場合は除外
                filter(n < 5) %>% 
                dplyr::select(id_cleansed),
              by = "id_cleansed") %>% 
    group_by(id_cleansed, term) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% 
    tidytext::cast_dtm(document = "id_cleansed",
                       term = "term",
                       value = "count")
}

# Function for test running of LDA
runLDAtest <- function(dtm, 
                       K_from = 15, K_to = 60, K_by = 15, 
                       seed = 123, 
                       path_base = "data-manual/LDA-stopword-examination/", 
                       file_name = "topic-list-stopword-rm") {
  for (K in seq(K_from, K_to, K_by)) {
    topicModel <- topicmodels::LDA(dtm,
                                   k = K,
                                   method = "Gibbs",
                                   control = list(alpha = 50/K, 
                                                  iter = 1000, 
                                                  verbose = 25, 
                                                  seed = seed))
    name_path <- str_c(path_base,
                       Sys.time() %>% 
                         str_remove_all("-") %>% 
                         str_replace_all(c(" " = "_",
                                           ":" = "")) %>% 
                         str_remove(".{2}$"),
                       "_",
                       file_name,
                       "_seed",
                       seed,
                       "_K",
                       K,
                       ".xlsx"
    )
    as.data.frame(terms(topicModel, 20)) %>% 
      writexl::write_xlsx(name_path)
  }
}

#------------------------------------------------------------------------------

# DTMの作成
dtm_rm_stpw_basic <- makeDTMatrix(token_rm_stpw_basic)
dtm_rm_stpw_general <- makeDTMatrix(token_rm_stpw_general)
dtm_rm_stpw_original <- makeDTMatrix(token_rm_stpw_original)

#------------------------------------------------------------------------------

# 2023-05-29
runLDAtest(dtm_rm_stpw_basic, file_name = "stopword-basic")
runLDAtest(dtm_rm_stpw_general, file_name = "stopword-general")
runLDAtest(dtm_rm_stpw_original, file_name = "stopword-original")

# Trial and error of original stop words---------------------------------------
# K = 30でチェックしてみる。

token_rm_stpw_original <- read_csv("data/token-05-general-stopword-removed.csv")
stopword_noun_original <- 
  read_excel("data-manual/token-stopword-selected.xlsx", sheet = 1)
stopword_noun_original %<>% filter(!(is.na(judge)))
stopword_noun_original %>% 
  group_by(judge, sub) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

## Pattern I:
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-I")

## Pattern II: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-II")

## Pattern III: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-III")

## Pattern IV: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       sub == "top") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-IV")

## Pattern V: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-V")

## Pattern VI: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       sub == "top" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VI")

## Pattern VII: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VII")

## Pattern VIII: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VIII")

## Pattern VIX: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VIX")

## Pattern VX: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       sub == "top" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VX")

## Pattern VXI: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VXI")

## Pattern VXII: 
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VXII")

## Pattern VXIII: 決定！
dtm_rm_stpw_original <- token_rm_stpw_original %>% 
  anti_join(stopword_noun_original %>% 
              filter(judge == "general" |
                       judge == "unclear" |
                       sub == "searched" |
                       sub == "lower" |
                       sub == "middle" |
                       sub == "upper" |
                       sub == "top" |
                       judge == "place") %>% 
              dplyr::select(term),
            by = "term") %>%
  makeDTMatrix()

runLDAtest(dtm_rm_stpw_original, K_from = 30, K_to = 30, K_by = 30,
           file_name = "stopword-original-VXIII")
