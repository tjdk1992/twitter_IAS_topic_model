#-----------------------------------------------------------------------------#
# Script Name: 07-token-cleaning.R                                            #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script cleanse tokens by removing unrelated parts of documents.      #
#-----------------------------------------------------------------------------#


# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse,
               magrittr,
               stringi,
               zipangu
               )

# Data
token_cleaning <- read_csv("data-proc/tokens-01-created.csv")

# Helper function -------------------------------------------------------------

# View top-occurred terms
viewTopTerm <- function(df_tokens, head = 5000) {
  df_tokens %>% 
    group_by(term) %>% 
    summarise(n = n()) %>% 
    arrange(term) %>% 
    head(head) %>% 
    View()
}

# Automatic cleaning ----------------------------------------------------------

# 品詞の選択
token_cleaning %<>% 
  filter(hinshi == "名詞" |
           hinshi == "動詞" |
           hinshi == "形容詞")

# 日本語を正規化する。
token_cleaning %<>% 
  mutate(term = stri_trans_general(term, "Halfwidth-Fullwidth"),
         term = str_jnormalize(term))

# 細々した不用品のクリーニング
token_cleaning %<>% 
  mutate(
    # 半角にしておく
    term = stri_trans_general(term, "Fullwidth-Halfwidth"),
    # 記号の除外
    term = str_replace_all(term, "[[:punct:]]", ""),
    # 数字の除外
    term = str_replace_all(term, "[[:digit:]]", ""),
    # 英語の除外
    term = str_replace_all(term, "[[:lower:]]", ""),
    term = str_replace_all(term, "[[:upper:]]", ""),
    # 空白、タブ、改行の除外
    term = str_replace_all(term, "[[:space:]]", ""),
    # 絵文字の除外
    term = stri_replace_all(term, regex = "\\p{Emoji}", "")) %>%
  filter(!(term == "")) %>% 
  filter(!(term == " ")) %>% 
  filter(!(term == "　"))

# 日本語を扱う上で半角カタカナが邪魔なので全角に変換しておく
token_cleaning %<>% 
  mutate(term = stri_trans_general(term, "Halfwidth-Fullwidth"))

# 先頭が"ー"のものを削除
token_cleaning %<>% mutate(term = str_replace_all(term, "^ー*", ""))

# Manual cleaning (transformative) --------------------------------------------

# 出現頻度が1のものは必ずノイズになるので除外する。
token_cleaning %<>% 
  anti_join(token_cleaning %>% 
              group_by(term) %>% 
              summarise(n = n()) %>% 
              filter(n == 1) %>% 
              dplyr::select(term),
            by = "term")

## 記号は意味の解釈が困難なので除外
token_cleaning %<>% 
  anti_join(token_cleaning %>% 
              group_by(term) %>% 
              summarise(n = n()) %>% 
              arrange(term) %>% 
              head(407) %>% 
              dplyr::select(term), 
            by = "term")

## 漢字以外の1文字のtermは除外する。
# all type of kana
kana_all <- unique(c(kana(type = "hira", dakuon = TRUE),
                     kana(type = "hira", handakuon = TRUE),
                     kana(type = "hira", kogaki = TRUE),
                     kana(type = "kata", dakuon = TRUE),
                     kana(type = "kata", handakuon = TRUE),
                     kana(type = "kata", kogaki = TRUE)))
# 全角の濁音がついてるやつ
kana_all <- c(kana_all, str_c(kana_all, "゛"), str_c(kana_all, "゙"))

# 漢字以外の1文字のトークンは除外する。
for (i in 1:length(kana_all)) {
  token_cleaning %<>% 
    filter(term != kana_all[i])
}

# 繰り返しは鳴き声だったりして意味がないので除外
# チェック
filter(token_cleaning, str_detect(term, "(.)\\1{2,}")) %>% pull(term)
filter(token_cleaning, str_detect(term, "(..)\\1{2,}")) %>% pull(term)
filter(token_cleaning, str_detect(term, "(...)\\1{2,}")) %>% pull(term)
filter(token_cleaning, str_detect(term, "(....)\\1{2,}")) %>% pull(term)

# 2回の繰り返しは意味が取れる可能性があるのでここでは残す（例：オオコウモリ）。
token_cleaning %<>% 
  filter(!(str_detect(term, "(.)\\1{2,}"))) %>% 
  filter(!(str_detect(term, "(..)\\1{2,}"))) %>%
  filter(!(str_detect(term, "(...)\\1{2,}"))) %>%
  filter(!(str_detect(term, "(....)\\1{2,}")))

# クリーニングが完了したデータを書き出し
write_csv(token_cleaning, "data-proc/tokens-02-cleaned.csv")
