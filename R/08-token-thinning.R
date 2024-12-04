#-----------------------------------------------------------------------------#
# Script Name: 08-token-thinning.R                                            #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script cuts major and minor terms in frequency for LDA inference.    #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse,
               magrittr
               )

# Data
token_thinning <- read_csv("data-proc/tokens-02-cleaned.csv")

# 出現頻度による単語の除外 ----------------------------------------------------

# 出現頻度のチェック
N_doc <- length(unique(token_thinning$id_cleansed))
tokens_summary_check <- token_thinning %>% 
  group_by(term) %>% 
  summarise(n = n(), freq = n / N_doc * 100) %>% 
  as.data.frame()

# 出現頻度の高い単語を確認
arrange(filter(tokens_summary_check, freq > 50), desc(freq))
arrange(filter(tokens_summary_check, freq > 25), desc(freq))
arrange(filter(tokens_summary_check, freq > 20), desc(freq))
arrange(filter(tokens_summary_check, freq > 15), desc(freq)) # この辺かな…
arrange(filter(tokens_summary_check, freq > 10), desc(freq))
arrange(filter(tokens_summary_check, freq > 5), desc(freq))

# 出現頻度の低い単語を確認
arrange(filter(tokens_summary_check, freq < 0.01), desc(freq))
arrange(filter(tokens_summary_check, freq < 0.05), desc(freq))
arrange(filter(tokens_summary_check, freq < 0.001), desc(freq)) # この辺かな
arrange(filter(tokens_summary_check, freq < 0.0005), desc(freq))
arrange(filter(tokens_summary_check, freq < 0.0006), desc(freq))
arrange(filter(tokens_summary_check, freq < 0.0007), desc(freq))
arrange(filter(tokens_summary_check, freq < 0.0008), desc(freq))
arrange(filter(tokens_summary_check, freq < 0.0009), desc(freq))

# 出現頻度の低い単語を除外する（n < 5）。
token_thinning %<>% 
  anti_join(tokens_summary_check %>% 
              # filter(freq > 15 | freq < 0.001) %>% 
              filter(freq > 10 | freq < 0.01) %>% 
              dplyr::select(term), 
            by = "term")

# 前処理が完了したデータの書き出し
write_csv(token_thinning, "data-proc/token-03-thinned.csv")
