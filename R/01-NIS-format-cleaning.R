#-----------------------------------------------------------------------------#
# Script Name: 01-nis-list-compilation.R                                      #
#                                                                             #
# Author: Daiki Tomojiri                                                      #
# Email: tomojiri.daiki@gmail.com                                             #
#                                                                             #
# This R script organize the format of NIS dataset.                           #
#-----------------------------------------------------------------------------#

# Setup -----------------------------------------------------------------------

# Initialization
rm(list = ls())
gc(); gc();

# Package
pacman::p_load(tidyverse,
               readxl,
               stringi,
               writexl
)

# Data
NIS_raw <- read_excel("data-raw/NIS-DB-finalized.xlsx", sheet = "target")

# Data organization -----------------------------------------------------------

# Rename columns and create katakana
NIS_all <- NIS_raw %>% 
  mutate(num_code1 = c(rep("0", 9), rep("", count(NIS_raw)-9)),
         num_code2 = c(rep("0", 99) , rep("", count(NIS_raw)-99))) %>% 
  transmute(code_ias = str_c("IAS", num_code1, num_code2, row_number()),
            group_biol = str_replace_all(Taxon,
                                         c("\\(insect\\)" = "",
                                           "\\(other\\)" = "")),
            taxon = 高次分類群,
            name_ja = 和名,
            KATAKANA = カタカナ名,
            name_ja_katakana = stri_trans_general(name_ja, "Fullwidth-Halfwidth"),
            katakana = stri_trans_general(KATAKANA, "Fullwidth-Halfwidth"),
            name_sp = 学名,
            reg1 = 法的扱い,
            reg2 = 外来生物法)

# Export cleaned data ---------------------------------------------------------

# Write in csv format
write_csv(NIS_all, "data/NIS-compiled.csv")
