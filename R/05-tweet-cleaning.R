#-----------------------------------------------------------------------------#
# Script Name: 05-tweet-cleaning.R                                            #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script cleanses tweets by removing some specific proper terms.       #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Packages
pacman::p_load(tidyverse, # for data manipulation
               magrittr,  # to overwrite data
               lubridate, # to handle date class
               stringi,   # to handle strings
               writexl,   # to write excel sheets
               readxl     # to read excel sheets
               )

# Data
tweet_cleaning <- read_csv("data-proc/tweet-03-screened.csv")
NIS_all <- read_csv("data/NIS-compiled.csv")

# Prepare basic data-----------------------------------------------------------

# Select necessary columns
tweet_cleaning %<>% 
  dplyr::select(id_screened,
                name_ja, 
                name_sp,
                date, 
                year, 
                text)

# 扱いやすくするために可能なものは半角にしてスペースは除外する。
tweet_cleaning %<>% 
  mutate(text = stringi::stri_trans_general(text, "Fullwidth-Halfwidth"))

# Cleansing tweet--------------------------------------------------------------

# tweetの確認
tweet_cleaning %>% 
  sample_n(1000) %>% 
  pull(text)

# -----------------------------------------------------------------------------
# 除外すべきは以下のものたち
# 日付、メタ文字、URL、メディア名、検索ワード、ハッシュタグ、メンションタグ
# -----------------------------------------------------------------------------

# Remove hash tags, mention tags, meta-chr and spaces
tweet_cleaning %<>% 
  mutate(text = str_replace_all(text, c("#[^\\s#]*" = "", # Hash tag
                                        "@[^\\s@]*" = "", # Mention
                                        "\n" = "",
                                        "\r" = "",
                                        "\t" = "",
                                        "[:blank:]" = "")))
tweet_cleaning
# Remove date and other general expression
tweet_cleaning %<>% 
  mutate(text = str_replace_all(
    text, c("\\([月火水木金土日]\\)" = "", # 曜日
            "[月火水木金土日]曜日" = "",
            "[0-9]{4}年([1-9]|0[1-9]|1[0-2])月([1-9]|0[1-9]|[12][0-9]|3[01])日" = "",
            "([1-9]|0[1-9]|1[0-2])月([1-9]|0[1-9]|[12][0-9]|3[01])日" = ""
    )))  # 改行文字

# Remove searched query
term_query <- c(NIS_all %>% 
                  mutate(n_chr = nchar(KATAKANA)) %>% 
                  arrange(desc(n_chr)) %>% 
                  mutate(KATAKANA = 
                           stri_trans_general(KATAKANA, 
                                              "Fullwidth-Halfwidth")) %>% 
                  pull(KATAKANA),
                "特定外来生物", "要注意外来生物", "侵略的外来種",
                "外来種", "移入種", "帰化種",
                "外来", "侵入", "移入", "帰化")

# 削除
for (i in 1:length(term_query)) {
  term_rm <- term_query[i]
  tweet_cleaning %<>% mutate(text = str_replace_all(text, term_rm, ""))
}

# メディアの引用が多いので、この時点で除外しておく。
## 記事等のタイトルの除外について
## Quotedと同じニュース記事の引用をみつける
## 引用記事タイトルはその人の発言としてみなしづらいので除外。
## URLを除外してsummarizeした時に複数出てくる文字列を抽出して、
## 複数あるもの（個別のコメントが含まれる場合にはgroupingされないとみなす）
## をテキスト中から除外する。

## URLを貼らずに記事のタイトルだけを貼ったようなtweetが散見されるため。
## 抽出したstringsの除外はtweet_normal_wo_URLにも適用する。
str_rm <- tweet_cleaning %>% 
  filter(str_detect(text, "https?://.*[a-zA-Z0-9]")) %>% 
  mutate(text = str_replace_all(text, 
                                c("https?://.*[a-zA-Z0-9]" = "",
                                  "^[:“”-]+" = "",
                                  "[…:“”\\|\\.]+$" = "",
                                  "から$" = "",
                                  "より$" = ""))) %>% 
  filter(text != "") %>% 
  group_by(text) %>% 
  summarise(n = n())

# View the results
str_rm %>% 
  mutate(n_chr = nchar(text)) %>% 
  arrange(n_chr) %>% 
  dplyr::select(n, n_chr, text) # %>% View()

# メディア名と本文を分解してよりessentialな状態にする。

## |で分割
str_rm %<>% 
  dplyr::select(-n) %>% 
  # まずは|で分割
  separate(text, sep = "\\|", into = str_c("SEP", seq(50)), fill = "right") %>%
  pivot_longer(cols = SEP1:SEP50,
               names_to = "SEP",
               values_to = "text") %>%
  mutate(text = if_else(SEP != "SEP1", str_c("|", text), text),
         text = str_replace_all(text, "\\|$", "")) %>% 
  dplyr::select(-SEP) %>%
  filter(!(is.na(text))) %>%
  # 続いて-で分割
  separate(text, sep = "-", into = str_c("SEP", seq(50)), fill = "right") %>% 
  pivot_longer(cols = SEP1:SEP50,
               names_to = "SEP",
               values_to = "text") %>% 
  mutate(text = if_else(SEP != "SEP1", str_c("-", text), text),
         text = str_replace_all(text, "-$", "")) %>% 
  dplyr::select(-SEP) %>% 
  filter(!(is.na(text))) %>% 
  # :で分割
  separate(text, sep = ":", into = str_c("SEP", seq(50)), fill = "right") %>% 
  pivot_longer(cols = SEP1:SEP50,
               names_to = "SEP",
               values_to = "text") %>% 
  mutate(text = if_else(SEP != "SEP1", str_c(":", text), text),
         text = str_replace_all(text, ":$", "")) %>% 
  dplyr::select(-SEP) %>% 
  filter(!(is.na(text))) %>% 
  # 続いて後方から最短一致で()内を分割
  mutate(text = str_replace_all(text, "\\((.*?)\\)$", "＿\\(\\1\\)")) %>% 
  mutate(text = str_replace_all(text, "\\[(.*?)\\]$", "＿\\[\\1\\]")) %>% 
  mutate(text = str_replace_all(text, "【(.*?)】$", "＿【\\1】")) %>% 
  mutate(text = str_replace_all(text, "^\\((.*?)\\)", "\\(\\1\\)＿")) %>% 
  mutate(text = str_replace_all(text, "^\\[(.*?)\\]", "\\[\\1\\]＿")) %>% 
  mutate(text = str_replace_all(text, "^【(.*?)】", "【\\1】＿")) %>% 
  separate(text, sep = "＿", into = str_c("SEP", seq(50)), fill = "right") %>% 
  pivot_longer(cols = SEP1:SEP50,
               names_to = "SEP",
               values_to = "text") %>% 
  dplyr::select(-SEP) %>% 
  filter(!(is.na(text)))

# Manualで選ぶ
str_rm %>% 
  mutate(text = stri_trans_general(text, "Halfwidth-Fullwidth"),
         n_chr = nchar(text)) %>% 
  filter(n_chr > 2, n_chr < 20) %>% 
  filter(!(str_detect(text, "。"))) %>%
  filter(!(str_detect(text, "、"))) %>%
  filter(!(str_detect(text, "「"))) %>%
  filter(!(str_detect(text, "」"))) %>%
  filter(!(str_detect(text, "った"))) %>%
  filter(!(str_detect(text, "って"))) %>%
  filter(!(str_detect(text, "です"))) %>%
  filter(!(str_detect(text, "ます"))) %>%
  filter(!(str_detect(text, "ました"))) %>%
  filter(!(str_detect(text, "ません"))) %>%
  filter(!(str_detect(text, "ください"))) %>%
  filter(!(str_detect(text, "である"))) %>%
  filter(!(str_detect(text, "^が"))) %>%
  filter(!(str_detect(text, "^の"))) %>%
  filter(!(str_detect(text, "^を"))) %>%
  filter(!(str_detect(text, "[都道府県市区町村]$"))) %>%
  mutate(text = stri_trans_general(text, "Fullwidth-Halfwidth")) %>% 
  group_by(text) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  arrange(desc(n)) %>%
  mutate(text2 = text,
         text = str_replace_all(text, c("^-" = "",
                                        "^:" = "",
                                        "^\\|" = "",
                                        "^[\\(\\)【】\\[\\]]" = "",
                                        "[\\(\\)【】\\[\\]]$" = "")),
         text_check = str_replace_all(text, c("[:punct:]" = "",
                                              "[:blank:]" = "")),
         judge = "NA") %>% 
  filter(!(str_detect(text, "^[0-9a-zA-Z]*$"))) %>% 
  dplyr::select(-text_check) %>% 
  mutate(n_chr = nchar(text)) %>% 
  filter(n != 1) %>%
  write_xlsx("data-manual/tweet-cleaning-unselected.xlsx")
  # 記事関連を全て除外

str_rm <- read_excel("data-manual/tweet-cleaning-selected.xlsx")

list_str_rm <- rbind(
  str_rm %>% 
    filter(judge == 1) %>% 
    transmute(text), 
  str_rm %>% 
    filter(judge == 2) %>% 
    transmute(text = text2)) %>% 
  mutate(n_chr = nchar(text)) %>% 
  arrange(desc(n_chr)) %>% 
  pull(text)

for (i in 1:length(list_str_rm)) {
  pattern = list_str_rm[i]
  tweet_cleaning %<>% 
    mutate(text = str_replace_all(text, pattern = stringr::fixed(pattern), ""))
}

# データの書き出し
tweet_cleaning %>% 
  mutate(
    # URLの削除
    text = str_replace_all(text, "https?://.*[a-zA-Z0-9]", ""),
    # 全角に戻しておく
    text = stringi::stri_trans_general(text, "Halfwidth-Fullwidth"),
    # idの付与
    id_cleansed = row_number()
  ) %>% 
  # 書き出し
  write_csv("data-proc/tweet-04-cleaned.csv")
