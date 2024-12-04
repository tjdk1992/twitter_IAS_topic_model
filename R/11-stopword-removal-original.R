#-----------------------------------------------------------------------------#
# Script Name: 11-stopword-removal-original.R                                 #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script removes stop words that selected manually from terms          #
# which are frequently occur in the retrieved tweets.                         #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Packages
pacman::p_load(tidyverse,
               magrittr,
               writexl,
               readxl)

# Data
token_rm_stpw_original <- read_csv("data-proc/token-05-general-stopword-removed.csv")
token_created <- read_csv("data/token-01-created.csv")
#------------------------------------------------------------------------------

# Original stopwordsの除外（名詞）

## ストップワード抽出用の単語リストを作成
## どの程度がいいんかな…
## 実際に分析に含まれたツイート数の取得
len_tweet <- nrow(token_rm_stpw_original %>% 
                    distinct(.keep_all = TRUE) %>% 
                    dplyr::select(id_cleansed) %>% 
                    distinct()) # *0.001

# 分布はどんな感じ？
token_rm_stpw_original %>% 
  distinct(.keep_all = TRUE) %>% 
  group_by(term) %>% 
  summarise(freq = n()) %>% 
  arrange(desc(freq)) %>% 
  mutate(freq_doc = freq / len_tweet * 100) %>% 
  filter(freq_doc > 0.1)

# 5%だったら上位22個
# 1%だったら上位224個
# 0.5%だったら上位561個
# 0.1%だったら上位2386個

# 分布を可視化してみる
token_rm_stpw_original %>% 
  distinct(.keep_all = TRUE) %>% 
  group_by(term) %>% 
  summarise(freq = n()) %>% 
  arrange(desc(freq)) %>% 
  mutate(freq_doc = freq / len_tweet * 100) %>% 
  filter(freq_doc > 0.5) %>% # 5%以上
  ggplot(aes(x = reorder(term, freq), y = freq)) +
  geom_bar(stat = "identity") +
  theme(axis.text = element_text(family = "HiraKakuPro-W3",
                                 size = 8,
                                 angle = 90))

# トップ2854位までの単語に含まれる名詞
token_rm_stpw_original %>% 
  distinct(.keep_all = TRUE) %>% 
  group_by(term) %>% 
  summarise(freq = n()) %>% 
  arrange(desc(freq)) %>% 
  mutate(freq_doc = freq / len_tweet * 100) %>% 
  filter(freq_doc > 0.1) %>%
  left_join(token_created %>% 
              dplyr::select(-id_cleansed) %>% 
              distinct(.keep_all = TRUE), by = "term") %>% 
  filter(hinshi == "名詞") %>% # 動詞と形容詞は自動でジャッジする
  mutate(judge = "") %>% 
  write_xlsx("data-manual/token-stopword-unselected.xlsx")

# Check the selected stop words
stopword_noun_original <- 
  read_excel("data-manual/token-stopword-selected.xlsx", sheet = 1)

# Summarise
stopword_noun_original %>% 
  group_by(judge, sub) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

# Pattern VXIII: lda-test-runの結果、決定
token_rm_stpw_original %<>% 
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
            by = "term")

# プラスアルファの試行錯誤
token_rm_stpw_original %<>% 
  filter(term != "今日") %>% 
  filter(term != "昨日") %>% 
  filter(term != "明日") %>% 
  filter(term != "本日") %>% 
  filter(term != "笑") %>% 
  filter(term != "もも") # 〜も〜も

# データの書き出し
write_csv(token_rm_stpw_original, "data-proc/token-06-original-stopword-removed.csv")
