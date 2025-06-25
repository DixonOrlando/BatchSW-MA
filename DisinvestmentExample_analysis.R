###############Disinvestment example###############

###Loading the packages###
library("lme4")
library("glmm")
library("swdpwr")
library("MASS")
library("glmmTMB")
library("meta")
library("tidyverse")
library("readxl")
library("lmerTest")
library(tidyverse)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(Matrix)
library(patchwork)
library(insight)

###Reading the data###
dis_data = read_excel("pmed.1002412.s002.xlsx")


###Preparation before analysis###

dis_data_final = dis_data %>%
  dplyr::select(no_we_exposure,
         sw_step,
         study1,
         dhvicavglos,
         site,
         index_ward)

dis_data_final = dis_data_final %>%
  filter(sw_step < 8,
         study1 == 1,
         index_ward %in% c(1,2,3,4,5,6,11,12,13,14,15,16))

#Adjusting the sw_step variable so it matches the calendar time based on the study.
dis_data_final = dis_data_final %>%
  mutate(sw_step_ori = sw_step,
         sw_step = ifelse(site == 2, sw_step +2, sw_step))

#Creating cluster by time variable so it can be used later for modelling purposes.
dis_data_final = dis_data_final %>%
  unite("clusbytime", index_ward, sw_step, sep = "c_t", remove = FALSE)

#Convert the variables into factor

dis_data_final = dis_data_final %>%
  mutate(sw_step = as.factor(sw_step),
         index_ward = as.factor(index_ward),
         site = as.factor(site),
         clusbytime = as.factor(clusbytime))

###Code for analysis###

for (k in 1:2) {
  if (k == 1) {
    dis_data_final = dis_data_final %>%
      mutate(outcome = dhvicavglos)
    
    row = 1
    
    analysis_results = data.frame(model = c("Model A", "Model B", "Model C", "Model D", "me.fixed", "me.random.reml.classic", "me.random.reml.HK", "me.random.reml.KR", "me.random.SJ.classic", "me.random.SJ.HK", "me.random.SJ.KR"), 
                                  transform_outcome = "No", 
                                  corr_struct = c(rep("Exchangeable", 11), rep("Block-exchangeble", 11), rep("DTD", 11)), 
                                  TE_estimate = 1, 
                                  lower_ci = 1, 
                                  upper_ci = 1, 
                                  p_val = 1)
    
    analysis_results = analysis_results %>%
      rbind(analysis_results %>%
              mutate(transform_outcome = "Yes"))
  }
  
  if (k != 1) {
    dis_data_final = dis_data_final %>%
      mutate(outcome = log(dhvicavglos+1))
  }
  
  #Type of linear mixed models:
  #A: shared period effects across batches & no treatment effect heterogeneity across batches
  #B: separate period effects across batches & no treatment effect heterogeneity across batches
  #C: shared period effects across batches & treatment effect heterogeneity across batches
  #D: separate period effects across batches & treatment effect heterogeneity across batches
  
  #Parameter information:
  #dhvicavglos: outcome
  #sw_step: time effect
  #index_ward: clusters
  #site: batch
  #no_we_exposure: exposure/treatment
  #clusbytime: cluster by time variable
  
  ###Model with exchangeable within-cluster correlation structure.
  
  #Model A:
  mod_A = lmer(outcome~no_we_exposure +sw_step + (1|index_ward), dis_data_final)
  analysis_results[row,4] = fixef(mod_A)[2]
  analysis_results[row, 5] = fixef(mod_A)[2] - coef(summary(mod_A))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_A)[2] + coef(summary(mod_A))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_A))["no_we_exposure",][5]
  row = row + 1
  
  #Model B:
  mod_B = lmer(outcome~no_we_exposure +  sw_step*site + (1|index_ward), dis_data_final)
  analysis_results[row,4] = fixef(mod_B)[2]
  analysis_results[row, 5] = fixef(mod_B)[2] - coef(summary(mod_B))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_B)[2] + coef(summary(mod_B))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_B))["no_we_exposure",][5]
  row = row + 1
  
  #Model C:
  mod_C = lmer(outcome~no_we_exposure +sw_step + (1|index_ward) + (-1+no_we_exposure| site), dis_data_final)
  analysis_results[row,4] = fixef(mod_C)[2]
  analysis_results[row, 5] = fixef(mod_C)[2] - coef(summary(mod_C))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_C)[2] + coef(summary(mod_C))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_C))["no_we_exposure",][5]
  row = row + 1
  
  #Model D:
  mod_D = lmer(outcome~no_we_exposure + sw_step*site + (1|index_ward) + (-1+no_we_exposure| site), dis_data_final)
  analysis_results[row,4] = fixef(mod_D)[2]
  analysis_results[row, 5] = fixef(mod_D)[2] - coef(summary(mod_D))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_D)[2] + coef(summary(mod_D))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_D))["no_we_exposure",][5]
  row = row + 1
  
  batchresults <- matrix(data=NA, nrow = 2, ncol = 2)
  
  if (k == 1) {
    var_results_exch = c(B = get_variance(mod_B), D = get_variance(mod_D))
  }
  
  if (k != 1) {
    pre = c(B = get_variance(mod_B), D = get_variance(mod_D))
    
    var_results_exch = rbind(var_results_exch, pre)
  }
  
  
  for (i in 1:2) {
    fit_batch <- lmer(outcome~no_we_exposure +sw_step + (1|index_ward), 
                      dis_data_final[dis_data_final$site == i,])
    
    batchresults[i,1] <- lme4::fixef(fit_batch)[2]
    batchresults[i,2] <- vcov(fit_batch)[2,2]
  }
  
  batchresults_exch = batchresults
  
  m.gen.fixed = metagen(TE = batchresults[,1],
                        seTE = sqrt(batchresults[,2]),
                        sm = "MD",
                        random = FALSE)
  analysis_results[row, 4] = m.gen.fixed$TE.common
  analysis_results[row, 5] = m.gen.fixed$lower.common
  analysis_results[row, 6] = m.gen.fixed$upper.common
  analysis_results[row, 7] = m.gen.fixed$pval.common
  row = row + 1
  
  
  m.gen.random1 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "classic" )
  analysis_results[row, 4] = m.gen.random1$TE.random
  analysis_results[row, 5] = m.gen.random1$lower.random
  analysis_results[row, 6] = m.gen.random1$upper.random
  analysis_results[row, 7] = m.gen.random1$pval.random
  row = row + 1
  
  m.gen.random2 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "HK" )
  analysis_results[row, 4] = m.gen.random2$TE.random
  analysis_results[row, 5] = m.gen.random2$lower.random
  analysis_results[row, 6] = m.gen.random2$upper.random
  analysis_results[row, 7] = m.gen.random2$pval.random
  row = row + 1
  
  m.gen.random3 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "KR" )
  analysis_results[row, 4] = m.gen.random3$TE.random
  analysis_results[row, 5] = m.gen.random3$lower.random
  analysis_results[row, 6] = m.gen.random3$upper.random
  analysis_results[row, 7] = m.gen.random3$pval.random
  row = row + 1
  
  m.gen.random4 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "classic" )
  analysis_results[row, 4] = m.gen.random4$TE.random
  analysis_results[row, 5] = m.gen.random4$lower.random
  analysis_results[row, 6] = m.gen.random4$upper.random
  analysis_results[row, 7] = m.gen.random4$pval.random
  row = row + 1
  
  m.gen.random5 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "HK" )
  analysis_results[row, 4] = m.gen.random5$TE.random
  analysis_results[row, 5] = m.gen.random5$lower.random
  analysis_results[row, 6] = m.gen.random5$upper.random
  analysis_results[row, 7] = m.gen.random5$pval.random
  row = row + 1
  
  m.gen.random6 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "KR" )
  analysis_results[row, 4] = m.gen.random6$TE.random
  analysis_results[row, 5] = m.gen.random6$lower.random
  analysis_results[row, 6] = m.gen.random6$upper.random
  analysis_results[row, 7] = m.gen.random6$pval.random
  row = row + 1
  
  #Parameter information:
  #dhvicavglos: outcome
  #sw_step: time effect
  #index_ward: clusters
  #site: batch
  #no_we_exposure: exposure/treatment
  #clusbytime: cluster by time variable
  
  ###Model with block-exchangeable within-cluster correlation structure.
  
  #Model A:
  mod_A = lmer(outcome~ no_we_exposure +sw_step + (1|index_ward) +(1|clusbytime), dis_data_final)
  analysis_results[row,4] = fixef(mod_A)[2]
  analysis_results[row, 5] = fixef(mod_A)[2] - coef(summary(mod_A))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_A)[2] + coef(summary(mod_A))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_A))["no_we_exposure",][5]
  row = row + 1
  
  #Model B:
  mod_B = lmer(outcome~ no_we_exposure + sw_step*site + (1|index_ward) +(1|clusbytime), dis_data_final)
  analysis_results[row,4] = fixef(mod_B)[2]
  analysis_results[row, 5] = fixef(mod_B)[2] - coef(summary(mod_B))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_B)[2] + coef(summary(mod_B))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_B))["no_we_exposure",][5]
  row = row + 1
  
  #Model C:
  mod_C = lmer(outcome~no_we_exposure +sw_step + (1|index_ward) +(1|clusbytime) + (-1+no_we_exposure| site), dis_data_final)
  analysis_results[row,4] = fixef(mod_C)[2]
  analysis_results[row, 5] = fixef(mod_C)[2] - coef(summary(mod_C))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_C)[2] + coef(summary(mod_C))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_C))["no_we_exposure",][5]
  row = row + 1
  
  #Model D:
  mod_D = lmer(outcome~no_we_exposure + sw_step*site + (1|index_ward) +(1|clusbytime) + (-1+no_we_exposure| site), dis_data_final)
  analysis_results[row,4] = fixef(mod_D)[2]
  analysis_results[row, 5] = fixef(mod_D)[2] - coef(summary(mod_D))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_D)[2] + coef(summary(mod_D))["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = coef(summary(mod_D))["no_we_exposure",][5]
  row = row + 1
  
  if (k == 1) {
    var_results_BE = c(B = get_variance(mod_B), D = get_variance(mod_D))
  }
  
  if (k != 1) {
    pre = c(B = get_variance(mod_B), D = get_variance(mod_D))
    
    var_results_BE = rbind(var_results_BE, pre)
  }
  
  batchresults <- matrix(data=NA, nrow = 2, ncol = 2)
  
  for (i in 1:2) {
    fit_batch <- lmer(outcome~no_we_exposure +sw_step + (1|index_ward) +(1|clusbytime), 
                      dis_data_final[dis_data_final$site == i,])
    
    batchresults[i,1] <- lme4::fixef(fit_batch)[2]
    batchresults[i,2] <- vcov(fit_batch)[2,2]
  }
  
  batchresults_BE = batchresults
  
  m.gen.fixed = metagen(TE = batchresults[,1],
                        seTE = sqrt(batchresults[,2]),
                        sm = "MD",
                        random = FALSE)
  analysis_results[row, 4] = m.gen.fixed$TE.common
  analysis_results[row, 5] = m.gen.fixed$lower.common
  analysis_results[row, 6] = m.gen.fixed$upper.common
  analysis_results[row, 7] = m.gen.fixed$pval.common
  row = row + 1
  
  
  m.gen.random1 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "classic" )
  analysis_results[row, 4] = m.gen.random1$TE.random
  analysis_results[row, 5] = m.gen.random1$lower.random
  analysis_results[row, 6] = m.gen.random1$upper.random
  analysis_results[row, 7] = m.gen.random1$pval.random
  row = row + 1
  
  m.gen.random2 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "HK" )
  analysis_results[row, 4] = m.gen.random2$TE.random
  analysis_results[row, 5] = m.gen.random2$lower.random
  analysis_results[row, 6] = m.gen.random2$upper.random
  analysis_results[row, 7] = m.gen.random2$pval.random
  row = row + 1
  
  m.gen.random3 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "KR" )
  analysis_results[row, 4] = m.gen.random3$TE.random
  analysis_results[row, 5] = m.gen.random3$lower.random
  analysis_results[row, 6] = m.gen.random3$upper.random
  analysis_results[row, 7] = m.gen.random3$pval.random
  row = row + 1
  
  m.gen.random4 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "classic" )
  analysis_results[row, 4] = m.gen.random4$TE.random
  analysis_results[row, 5] = m.gen.random4$lower.random
  analysis_results[row, 6] = m.gen.random4$upper.random
  analysis_results[row, 7] = m.gen.random4$pval.random
  row = row + 1
  
  m.gen.random5 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "HK" )
  analysis_results[row, 4] = m.gen.random5$TE.random
  analysis_results[row, 5] = m.gen.random5$lower.random
  analysis_results[row, 6] = m.gen.random5$upper.random
  analysis_results[row, 7] = m.gen.random5$pval.random
  row = row + 1
  
  m.gen.random6 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "KR" )
  analysis_results[row, 4] = m.gen.random6$TE.random
  analysis_results[row, 5] = m.gen.random6$lower.random
  analysis_results[row, 6] = m.gen.random6$upper.random
  analysis_results[row, 7] = m.gen.random6$pval.random
  row = row + 1
  
  #Parameter information:
  #dhvicavglos: outcome
  #sw_step: time effect
  #index_ward: clusters
  #site: batch
  #no_we_exposure: exposure/treatment
  #clusbytime: cluster by time variable
  
  ###Model with discrete time decay within cluster correlation structure. 
  
  #Model A:
  mod_A = glmmTMB(outcome ~ no_we_exposure + sw_step + ar1(sw_step + 0 | index_ward), REML=T, 
                  data = dis_data_final, family = gaussian)
  analysis_results[row,4] = fixef(mod_A)$cond[2]
  analysis_results[row, 5] = fixef(mod_A)$cond[2] - (coef(summary(mod_A))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_A)$cond[2] + (coef(summary(mod_A))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = (coef(summary(mod_A))$cond)["no_we_exposure",][4]
  row = row + 1
  
  #Model B:
  mod_B = glmmTMB(outcome ~ no_we_exposure + sw_step*site + ar1(sw_step + 0 | index_ward), REML=T, 
                  data = dis_data_final, family = gaussian)
  analysis_results[row,4] = fixef(mod_B)$cond[2]
  analysis_results[row, 5] = fixef(mod_B)$cond[2] - (coef(summary(mod_B))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_B)$cond[2] + (coef(summary(mod_B))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = (coef(summary(mod_B))$cond)["no_we_exposure",][4]
  row = row + 1
  
  #Model C:
  mod_C = glmmTMB(outcome ~ no_we_exposure + sw_step + ar1(sw_step + 0 | index_ward)
                  + (no_we_exposure + 0| site), REML=T, 
                  data = dis_data_final, family = gaussian)
  analysis_results[row,4] = fixef(mod_C)$cond[2]
  analysis_results[row, 5] = fixef(mod_C)$cond[2] - (coef(summary(mod_C))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_C)$cond[2] + (coef(summary(mod_C))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = (coef(summary(mod_C))$cond)["no_we_exposure",][4]
  row = row + 1
  
  #Model D:
  mod_D = glmmTMB(outcome ~ no_we_exposure + sw_step*site + ar1(sw_step + 0 | index_ward)
                  + (no_we_exposure + 0| site), REML=T, 
                  data = dis_data_final, family = gaussian)
  
  analysis_results[row,4] = fixef(mod_D)$cond[2]
  analysis_results[row, 5] = fixef(mod_D)$cond[2] - (coef(summary(mod_D))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row,6] = fixef(mod_D)$cond[2] + (coef(summary(mod_D))$cond)["no_we_exposure",][2] * (-qnorm(0.025))
  analysis_results[row, 7] = (coef(summary(mod_D))$cond)["no_we_exposure",][4]
  row = row + 1
  
  if (k == 1) {
    var_results_DTD = c(B = get_variance(mod_B), D = get_variance(mod_D))
  }
  
  if (k != 1) {
    pre = c(B = get_variance(mod_B), D = get_variance(mod_D))
    
    var_results_DTD = rbind(var_results_DTD, pre)
  }
  
  batchresults <- matrix(data=NA, nrow = 2, ncol = 2)
  
  for (i in 1:2) {
    fit_batch <- glmmTMB(outcome ~ no_we_exposure + sw_step + ar1(sw_step + 0 | index_ward), REML=T, 
                         data = dis_data_final[dis_data_final$site == i,], family = gaussian)
    
    batchresults[i,1] <- fixef(fit_batch)$cond[2]
    batchresults[i,2] <- vcov(fit_batch)$cond[2,2]
  }
  
  batchresults_DTD = batchresults
  
  m.gen.fixed = metagen(TE = batchresults[,1],
                        seTE = sqrt(batchresults[,2]),
                        sm = "MD",
                        random = FALSE)
  analysis_results[row, 4] = m.gen.fixed$TE.common
  analysis_results[row, 5] = m.gen.fixed$lower.common
  analysis_results[row, 6] = m.gen.fixed$upper.common
  analysis_results[row, 7] = m.gen.fixed$pval.common
  row = row + 1
  
  
  m.gen.random1 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "classic" )
  analysis_results[row, 4] = m.gen.random1$TE.random
  analysis_results[row, 5] = m.gen.random1$lower.random
  analysis_results[row, 6] = m.gen.random1$upper.random
  analysis_results[row, 7] = m.gen.random1$pval.random
  row = row + 1
  
  m.gen.random2 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "HK" )
  analysis_results[row, 4] = m.gen.random2$TE.random
  analysis_results[row, 5] = m.gen.random2$lower.random
  analysis_results[row, 6] = m.gen.random2$upper.random
  analysis_results[row, 7] = m.gen.random2$pval.random
  row = row + 1
  
  m.gen.random3 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "REML",
                           method.random.ci = "KR" )
  analysis_results[row, 4] = m.gen.random3$TE.random
  analysis_results[row, 5] = m.gen.random3$lower.random
  analysis_results[row, 6] = m.gen.random3$upper.random
  analysis_results[row, 7] = m.gen.random3$pval.random
  row = row + 1
  
  m.gen.random4 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "classic" )
  analysis_results[row, 4] = m.gen.random4$TE.random
  analysis_results[row, 5] = m.gen.random4$lower.random
  analysis_results[row, 6] = m.gen.random4$upper.random
  analysis_results[row, 7] = m.gen.random4$pval.random
  row = row + 1
  
  m.gen.random5 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "HK" )
  analysis_results[row, 4] = m.gen.random5$TE.random
  analysis_results[row, 5] = m.gen.random5$lower.random
  analysis_results[row, 6] = m.gen.random5$upper.random
  analysis_results[row, 7] = m.gen.random5$pval.random
  row = row + 1
  
  m.gen.random6 <- metagen(TE = batchresults[,1],
                           seTE = sqrt(batchresults[,2]),
                           sm = "MD",
                           common = FALSE,
                           method.tau = "SJ",
                           method.random.ci = "KR" )
  analysis_results[row, 4] = m.gen.random6$TE.random
  analysis_results[row, 5] = m.gen.random6$lower.random
  analysis_results[row, 6] = m.gen.random6$upper.random
  analysis_results[row, 7] = m.gen.random6$pval.random
  row = row + 1
  
}

write_csv(analysis_results, "disinvestment_analysis_results_19-06-2025.csv")


#Plotting the results

analysis_results = read_csv("disinvestment_analysis_results_19-06-2025.csv")

analysis_results = analysis_results %>%
  mutate(model = case_when(
    model == "me.fixed" ~ "ma.fixed",
    model == "me.random.reml.classic" ~ "ma.random.reml.classic",
    model == "me.random.reml.HK" ~ "ma.random.reml.HK",
    model == "me.random.reml.KR" ~ "ma.random.reml.KR",
    model == "me.random.SJ.classic" ~ "ma.random.SJ.classic",
    model == "me.random.SJ.HK" ~ "ma.random.SJ.HK",
    model == "me.random.SJ.KR" ~ "ma.random.SJ.KR",
    TRUE ~ model
  ))

#analysis_results$model <- factor(analysis_results$model, levels = c("Model A", "Model B", "Model C", "Model D", "ma.fixed", "ma.random.reml.classic", "ma.random.reml.HK", "ma.random.reml.KR", "ma.random.SJ.classic", "ma.random.SJ.HK", "ma.random.SJ.KR"))

#Sorting the levels so that when plotted it can be ordered nicely.
analysis_results$model = factor(analysis_results$model, levels = sort(unique(analysis_results$model), decreasing = T))

analysis_results = analysis_results %>%
  dplyr::mutate(corr_struct = if_else(corr_struct == "Block-exchangeble", "Block-exchangeable", corr_struct))

analysis_results$corr_struct = factor(analysis_results$corr_struct, levels = c("Exchangeable", "Block-exchangeable", "DTD"))


###Plot to show###

#Without log-transformation
analysis_results %>%
  filter(transform_outcome == "No") %>%
  ggplot(aes(y = model, x = TE_estimate)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci)) +
  facet_grid(corr_struct~.) +
  labs(x = "Treatment Effect Estimate", y = "Model Type") #+
  #ggtitle("Disinvesment example: analysis using non-log-transformed outcome")

#With log-transformation
analysis_results %>%
  filter(transform_outcome == "Yes") %>%
  ggplot(aes(y = model, x = TE_estimate)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci)) +
  facet_grid(corr_struct~.) +
  labs(x = "Treatment Effect Estimate", y = "Model Type") #+
  #ggtitle("Disinvesment example: analysis using log-transformed outcome")




############Creating schema for disinvestment example to show participants per cluster period############

SWdesmat <- function(Ts) {
  Xsw <- matrix(data=0, ncol = (Ts), nrow = (Ts-1))
  for(i in 1:(Ts-1)) {
    Xsw[i,(i+1):Ts] <- 1
  }
  return(Xsw)
}
#Batched SW schematic
batchSWscheme <- function(olap, S, Ts,  Kseq){
  #Creates a batched SW design schematic with
  #S batches of Ts-period stepped wedge designs, 
  #with overlap between batches of olap periods
  #and K clusters in each sequence
  
  g<- Ts- olap #gap between each batch and the next
  
  SWbatch <- matrix(data=NA, nrow=S*(Ts-1), ncol=Ts + (S-1)*(Ts-olap))  
  
  for(i in 1:S){
    SWbatch[((i-1)*(Ts-1) + 1):(i*(Ts-1)) , 
            ((i-1)*g +1):((i-1)*g + Ts)] <- SWdesmat(Ts)
    
  }
  
  return( SWbatch[sort(rep(1:nrow(SWbatch), Kseq)), ])
  
  
}



Xdes00 <- batchSWscheme(olap=5, S=2, Ts=7,  Kseq=1)
melted_myx00 <- melt(Xdes00)
names(melted_myx00)[names(melted_myx00)=="Var1"] <- "Sequence"
names(melted_myx00)[names(melted_myx00)=="Var2"] <- "Period"

melted_myx00 = melted_myx00 %>%
  mutate(site = if_else(Sequence >=7, 2, 1),
         site = factor(as.character(site)))


data_join = dis_data_final %>% ###CORRECT ONE!
  mutate(sw_step = as.integer(as.character(sw_step)),
         index_ward = as.integer(as.character(index_ward)),
         site = as.integer(as.character(site))) %>% 
  mutate(index_ward_new = case_when(
    index_ward == 5 ~ 1,
    index_ward == 1 ~ 2,
    index_ward == 6 ~ 3,
    index_ward == 2 ~ 4,
    index_ward == 4 ~ 5,
    index_ward == 3 ~ 6,
    index_ward == 11 ~ 7,
    index_ward == 12 ~ 8,
    index_ward == 13 ~ 9,
    index_ward == 14 ~ 10,
    index_ward == 15 ~ 11,
    index_ward == 16 ~ 12,
  )) %>%
  dplyr::count(sw_step, index_ward_new, site, no_we_exposure) 
  
  ##%>% 
  #mutate(index_ward = if_else(site == 2, index_ward - 4, index_ward))

final_join = melted_myx00 %>% ####CORRECT ONE!
  dplyr::select(-site) %>%
  left_join(data_join %>% dplyr::select(-site), by = c("Period" = "sw_step", "Sequence" = "index_ward_new", "value" = "no_we_exposure")) 


color_palette1 <-colorRampPalette(c( "white", "grey"))(2)
ggplot(data =final_join, aes(x=Period, y=Sequence, fill = factor(value))) + ###FINAL PLOT!
  geom_tile( colour = "grey50") +
  scale_y_reverse(breaks=c(1:12)) +
  scale_x_continuous(breaks=c(1:12)) +
  theme(panel.grid.minor = element_blank()) +
  geom_text(aes(Period, Sequence, label = round(n,4)), color = "black", size = 5) +
  scale_fill_manual(values = color_palette1, na.value = theme_grey()$panel.background$fill) + 
  theme(legend.position="none") +
  labs(y = "Ward")







