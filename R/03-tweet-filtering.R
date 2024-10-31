#-----------------------------------------------------------------------------#
# Script Name: 03-tweet-filtering.R                                           #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script screens tweets for LDA by their attributes and cleanse tweet. #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse, # for data manipulation
               magrittr,  # to overwrite data
               lubridate, # to handle date class
               stringi,   # to handle strings
               zipangu,
               RMeCab
               )

# Data
tweet_filtering <- read_csv("data-proc/tweet-01-binded.csv")
NIS_all <- read_csv("data/NIS-compiled.csv")

# Prepare raw data ------------------------------------------------------------

# Data structure
str(tweet_filtering)

# Merge IAS species information
tweet_filtering %<>% 
  mutate(code_ias = str_replace_all(MyFile, 
                                    "data-raw/doc-tweet-ias/doc-tweet-", 
                                    "")) %>% 
  left_join(NIS_all %>% 
              dplyr::select(code_ias,
                            name_ja,
                            name_sp), 
            by = "code_ias")
# Check the result of merging
filter(tweet_filtering, is.na(name_sp)) # OK

# Select necessary columns
tweet_filtering %<>% 
  transmute(tweet_id = id, 
            author_id, 
            type,
            lang,
            code_ias,
            name_ja,
            name_sp,
            date = date(created_at), 
            year = year(date), 
            text)

# Screen the tweet ------------------------------------------------------------

# Tweetの確認
tweet_filtering 

# Remove tweets including ""invad*
tweet_filtering %<>% 
  filter(str_detect(text, "外来|移入|帰化"))

# Remove duplicate within individual species
tweet_filtering %<>% 
  mutate(rm_dup_sp = str_c(tweet_id, "_", name_sp)) %>% 
  distinct(rm_dup_sp, .keep_all = TRUE) %>% 
  dplyr::select(-rm_dup_sp) 

# Screen by language
tweet_filtering %<>% dplyr::select(-lang) # limit to Japanese

# Here is a raw number of tweets for this study
(n_tweet_original <- nrow(tweet_filtering))

# Screen by the Type
# Check the type
table(tweet_filtering$type, useNA = "always")

# Replace NA
tweet_filtering %<>% mutate(type = replace_na(type, "normal"))

# Remove retweeted tweet
tweet_filtering %<>% filter(type != "retweeted") 

# The number of removed retweets at this step
(n_tweet_second <- n_tweet_original - nrow(tweet_filtering))

#-----------------------------------------------------------------------------#
# If there are many identical tweets from a single user account, it is likely 
# that it is a promotion or a bot. In this case, only the very first tweet (in years) 
# is included in the analysis to avoid overestimation due to multiple tweets.
#-----------------------------------------------------------------------------#

# Summarize by user and text
tweet_filtering %>% 
  mutate(tweet_rm = str_c(author_id, "_", name_sp, "_", text)) %>% 
  group_by(tweet_rm) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  head(20)

# Remove duplicate
tweet_filtering %<>% 
  mutate(tweet_rm = str_c(author_id, "_", name_sp, "_", text)) %>% 
  group_by(tweet_rm, year) %>%
  mutate(id_rm = row_number()) %>% 
  ungroup() %>% 
  filter(id_rm == 1) %>% 
  dplyr::select(-id_rm) 

# The number of removed retweets at this step
(n_tweet_third <- n_tweet_second - nrow(tweet_filtering))

# Export data
tweet_filtering %>% 
  mutate(id_filtered = row_number()) %>% 
  write_csv("data-proc/tweet-02-filtered.csv")
