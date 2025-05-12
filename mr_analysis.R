
# Tom Bond; tom.bond@bristol.ac.uk
# 6th May 2025
# Code to repeat main MR analysis in Bond TA, Bhatta L, Yang Q, Moen G-H, Wang G, Beamont RN, et al. Parental body mass index and offspring cardiovascular risk factors in adulthood: an intergenerational Mendelian randomization study. medRxiv. 2025:2025.05.06.25327111.

# 1: setup ####

library(TwoSampleMR)
library(dplyr)
library(forcats)
library(openxlsx)
library(stringr)
rm(list = ls())
path <- '/path/to/input/data/'

# function to make clearer outcomes names etc. for MR results
format_res <- function(res){
  res$mat_pat_off2 = str_split_i(res$outcome, '_', 2)
  res$marg_cond = str_split_i(res$outcome, '_', 3)
  res$outcome_name = str_split_i(res$outcome, '_', 1)
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

d <- read.xlsx(paste0(path, 'parental_bmi_mr_SUPPLEMENTARY_TABLES.xlsx'),
               sheet = 'Supplementary_table_4', colNames = TRUE) %>%
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

# 6: convert results from (3) to SD scale ####

sds <- read.xlsx('https://github.com/tom-a-bond/mr_mat_pat_bmi_off_cvrf/raw/refs/heads/main/weighted_mean_pheno_sd_CORRECTED_GITHUB.xlsx')
# function to standardise MR results
standardise_res <- function(res){
  for(i in 1:nrow(res)){
    out <- res[i, 'outcome_name']
    if(!grepl('bw', out)){ # birth weight was already standardised prior to GWAS so should not be standardised again
      sd_use <- sds[which(sds$phen == out), 'weighted_mean_sd']
      res[i, 'b'] <- res[i, 'b'] / sd_use
      res[i, 'se'] <- res[i, 'se'] / sd_use
    }
  }
  res
}
res_mat_adj_sd <- standardise_res(res_mat_adj)

