#-----------------------------------------------------------------------------#
# Script Name: 10-stopword-removal-general.R                                 #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script removes stop words which can be automatically selected by     #
# using general tweets data.                                                  #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
library(tidyverse,       # for data manipulation
        academictwitteR, # for handling twitter
        RMeCab           # for tokenization of Japanese sentence
        )

# Color palette
pal_orig <- c(rep(pals::cols25(25), 2))

# Data
token_rm_stpw_basic <- read_csv("data-proc/token-04-basic-stopword-removed.csv")

# Twitter general stopword selection ---------------------------------------

# Retrieve random tweets
period_get_tweet <- tibble(y_start = c(2008:2022),
                             y_end = c(2009:2023),
                             time = "-01-01T00:00:00Z") %>% 
  transmute(start = str_c(y_start, time),
            end = str_c(y_end, time))

# "、"を含んでいれば文節が複数あるので動詞を含む可能性が高い。
for (i in 1:length(period_get_tweet$start)) {
  start_tweets <- unlist(period_get_tweet[i, 1], use.names = FALSE)
  end_tweets <- unlist(period_get_tweet[i, 2], use.names = FALSE)
  data_path <- paste0("data-raw/tweet-stopword/doc_", i)
  t_tweets_verb_est <-
    academictwitteR::get_all_tweets(query = "、",
                                    start_tweets = start_tweets,
                                    end_tweets = end_tweets,
                                    bearer_token = bearer_token,
                                    data_path = data_path,
                                    bind_tweets = FALSE,
                                    n = 10000)
}

# bind retrieved tweets for verb test
tweet_stopword <- tibble()
for (i in 1:length(period_get_tweet$start)) {
  data_path <- paste0("data-raw/tweet-stopword/doc_", i)
  t_tweet_stopword_binded <- bind_tweets(data_path = data_path, 
                                     output_format = "tidy")
  tweet_stopword <- rbind(tweet_stopword, t_tweet_stopword_binded)
}

# 形態素解析（DFじゃないと形態素解析してくれない）
df_tweet_stopword <- tweet_stopword %>% 
  transmute(id = row_number(), text) %>% 
  as.data.frame()

# MeCabによる形態素解析で単語のリストを取得
tic() # コードの実行時間を計測するコマンド
tweet_stopword_list <- RMeCabDF(df_tweet_stopword, "text", 1)
toc() # 100.964 sec elapsed

# リストをトークンのDFに変換
tic()
tweet_stopword_token <- purrr::pmap_df(list(nv = tweet_stopword_list,
                                                  title = df_tweet_stopword$id),
                                             function(nv, title){
                                               tibble(title = title,
                                                      term = nv,
                                                      hinshi = names(nv))
                                             })
toc() # 

# データの書き出し
write_csv(tweet_stopword_token, "data/token-general-tweet.csv")

#------------------------------------------------------------------------------

# My general stopwords
stopword_test_tokens <- read_csv("data/token-general-tweet.csv")

# 動詞・形容詞のストップワードの作成

# 実際に解析に含まれたツイート数
len_tweet_stopword <- nrow(stopword_test_tokens %>% 
                                  distinct(.keep_all = TRUE) %>% 
                                  dplyr::select(title) %>% 
                                  distinct())

# 動詞、名詞、形容詞のみ抽出したものを書き出しておく。
# 他の研究の解析にも使える。
stopword_test_tokens %>% 
  filter(hinshi == "動詞" | hinshi == "形容詞" | hinshi == "名詞") %>% 
  write_csv("data/df-id-tokens-general-tweet.csv")

# 動詞・形容詞の抽出
df_id_stop_verb_test_tokens <- stopword_test_tokens %>% 
  distinct(.keep_all = TRUE) %>% # 同一ツイート中の複数の同じ単語の除外
  filter(hinshi == "動詞")
df_id_stop_abj_test_tokens <- stopword_test_tokens %>% 
  distinct(.keep_all = TRUE) %>% 
  filter(hinshi == "形容詞")

# 動詞
stopword_verb_original <- df_id_stop_verb_test_tokens %>% 
  group_by(term) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  mutate(id = row_number(),
         freq = n / len_tweet_stopword * 100) %>% # n = 2599
  filter(freq > 0.05) %>% 
  dplyr::select(term)
View(stopword_verb_original)

# 形容詞
stopword_adj_original <- df_id_stop_abj_test_tokens %>% 
  group_by(term) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  mutate(id = row_number(),
         freq = n / len_tweet_stopword * 100) %>% # n = 2599
  filter(freq > 0.05) %>% 
  dplyr::select(term)
View(stopword_adj_original)

df_id_stop_verb_test_tokens %>% 
  filter(term == "です")

# 出現頻度 > 0.05の動詞・形容詞と
stopword_additional_original <- data.frame(
  term = c(
    # 敬語表現、
    "らっしゃる", "られる", "ます", "なさる", "さる", 
    # 敬称・接尾語、
    "くん", "君", "ちゃん", "氏", "さん",
    "さま", "様", "殿", "さん", "っち" ,
    # 代名詞
    "わたし", "私", "ぼく", "僕", "あなた", 
    "きみ", "かれ", "彼", "かのじょ", "彼女",
    "だれ", "誰", "どなた", "これ", "ここ",
    "こちら", "それ", "そこ", "そちら", "あれ",
    "あそこ", "あちら", "どれ", "どこ", "どちら"
  )
)

# ストップワードを除外したデータの書き出し
token_rm_stpw_basic %<>% 
  anti_join(stopword_verb_original, by = "term") %>% 
  anti_join(stopword_adj_original, by = "term") %>% 
  anti_join(stopword_additional_original, by = "term")

# Write data
write_csv(token_rm_stpw_basic, "data-proc/token-05-general-stopword-removed.csv")
