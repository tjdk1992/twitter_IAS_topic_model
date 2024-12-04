#-----------------------------------------------------------------------------#
# Script Name: 09-stopword-removal-basic.R                                    #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script removes basic stop words prepared on the Internet in advance. #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse,
               magrittr)

# Data
token_rm_stpw_basic <- read_csv("data-proc/tokens-03-thinned.csv")

#------------------------------------------------------------------------------

# General stopwordsの除外
## General stopwordのダウンロード
## (svn.sourceforge.jp/svnroot/slothlib/CSharp/Version1/
## SlothLib/NLP/Filter/StopWord/word/Japanese.txt)
stopword_general <- read.csv("resource/jpn_stopword.csv", 
                             header = T, fileEncoding = "UTF-8-BOM")

## General stopwordを除外
token_rm_stpw_basic %<>% 
  anti_join(stopword_general, by = "term")

## データの書き出し
write_csv(token_rm_stpw_basic, "data-proc/token-04-basic-stopword-removed.csv")
