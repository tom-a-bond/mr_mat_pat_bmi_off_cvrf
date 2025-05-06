
# Tom Bond; tom.bond@bristol.ac.uk
# 6th May 2025
# Code to accompany: Bond TA, Bhatta L, Yang Q, Moen G-H, Wang G, Beaumont R, et al. Parental body mass index and offspring cardiovascular risk factors in adulthood: an intergenerational Mendelian randomization study. medrXiv. 2025.

# 1: setup ####

library(TwoSampleMR)
library(readxl)
library(dplyr)
rm(list = ls())
path <- '/path/to/input/data/'

# function to make clearer outcomes names etc. for MR results
format_res <- function(res){
  res$mat_pat_off = unlist(lapply(str_split(res$outcome, fixed('_')), '[', 2))
  res$marg_cond = unlist(lapply(str_split(res$outcome, fixed('_')), '[', 3))
  res$outcome_name = unlist(lapply(str_split(res$outcome, fixed('_')), '[', 1))
  res$adj_unadj <- fct_recode(res$marg_cond, adjusted = 'cond', unadjusted = 'marg')
  res$marg_cond <- NULL
  res$mat_pat_off <- fct_recode(res$mat_pat_off, maternal = 'mat', paternal = 'pat', offspring = 'off')
  res
}

# 2: load harmonised GWAS data ####
# data are as per the primary analyses in the manuscript: parental genetic effects
# on offspring outcomes were adjusted for offspring genotype via the duos WLM,
# apart from for paternal BMI --> offspring birth weight, for which we used the
# trios WLM to simultaneously adjust for maternal genotype and offspring
# genotype; this was necessary to avoid collider bias

d <- read_xlsx(paste0(path, 'parental_bmi_mr_SUPPLEMENTARY_TABLES.xlsx'),
               sheet = 'Supplementary_table_S3', col_names = TRUE) %>%
  as.data.frame()
nrow(d) # 45116

# 3: run main MR analysis: maternal BMI -> offspring outcomes, adjusting for offspring genotype ####

mat_adj <- d[which(d$mat_pat_off == 'maternal' & d$adj_unadj == 'adjusted'), ]
res_mat_adj <- mr(mat_adj, method_list = 'mr_ivw')
res_mat_adj <- format_res(res_mat_adj)

# 4: run main MR analysis: paternal BMI -> offspring outcomes, adjusting for offspring genotype ####

pat_adj <- d[which(d$mat_pat_off == 'paternal' & d$adj_unadj == 'adjusted'), ]
res_pat_adj <- mr(pat_adj, method_list = 'mr_ivw')
res_pat_adj <- format_res(res_pat_adj)

# 5: rerun (3), without adjusting for offspring genotype ####

mat_unadj <- d[which(d$mat_pat_off == 'maternal' & d$adj_unadj == 'unadjusted'), ]
res_mat_unadj <- mr(mat_unadj, method_list = 'mr_ivw')
res_mat_unadj <- format_res(res_mat_unadj)








