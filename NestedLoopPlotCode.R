###Analysing the results of the simulation

library(tidyverse)
library(looplot)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

#Reading the results from the simulation study and preparing it for analysis.
all_contsresults = read_csv("all_contsresults_BE&E_final.csv")

all_contsresults = all_contsresults %>%
  dplyr::select(-`...1`)

all_contsresults_bias = all_contsresults %>%
  dplyr::select(model, Bias, Bias_minusSE, Bias_plusSE, T, nbatch, nolap, K, M, sharedtime, effsize, ICC, CAC, corrstruct, treateffvar, comb_ID) %>%
  mutate(with_corr = case_when(
    CAC == 1 ~ "Exch",
    CAC != 1 & corrstruct == 0 ~ "BE",
    CAC != 1 & corrstruct == 1 ~ "DTD"
  )) %>%
  pivot_wider(names_from = model,
              values_from = c(Bias, Bias_minusSE, Bias_plusSE)) 

all_contsresults_EmpSE = all_contsresults %>%
  dplyr::select(model, Emp_SE, EmpSE_minusSE, EmpSE_plusSE, T, nbatch, nolap, K, M, sharedtime, effsize, ICC, CAC, corrstruct, treateffvar, comb_ID) %>%
  mutate(with_corr = case_when(
    CAC == 1 ~ "Exch",
    CAC != 1 & corrstruct == 0 ~ "BE",
    CAC != 1 & corrstruct == 1 ~ "DTD"
  )) %>%
  pivot_wider(names_from = model,
              values_from = c(Emp_SE, EmpSE_minusSE, EmpSE_plusSE))

all_contsresults_coverage = all_contsresults %>%
  dplyr::select(model, Coverage, Coverage_minusSE, Coverage_plusSE, T, nbatch, nolap, K, M, sharedtime, effsize, ICC, CAC, corrstruct, treateffvar, comb_ID) %>%
  mutate(with_corr = case_when(
    CAC == 1 ~ "Exch",
    CAC != 1 & corrstruct == 0 ~ "BE",
    CAC != 1 & corrstruct == 1 ~ "DTD"
  )) %>%
  pivot_wider(names_from = model,
              values_from = c(Coverage, Coverage_minusSE, Coverage_plusSE))


########################All the meta-analysis plots first########################

#Prespecifying the steps for the nested loop plots.
steps = c("treateffvar", "M", "K",  "T", "CAC")
steps_names = c(expression(sigma[eta]^2), "m", "Clusters", "Period", "CAC")

######Bias plot######

bias_MA_separate_SJ = nested_loop_plot(resdf = all_contsresults_bias %>%
                                         mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                nolap = as.factor(as.character(nolap))) %>%
                                         filter(sharedtime == 1 & effsize == 0) %>%
                                         mutate(ICC = factor(ICC)),
                                       x = "ICC", 
                                       methods = c("Bias_MA.random.SJ.classic", "Bias_MA.random.SJ.HK", "Bias_MA.random.SJ.HK.SEadhoc", "Bias_MA.random.SJ.HK.IQadhoc", "Bias_MA.random.SJ.KR"),
                                       steps = steps,
                                       colors = brewer.pal(5, name = "Dark2"),
                                       steps_y_base = -0.017,
                                       steps_y_height = 0.005,
                                       steps_y_shift = 0.005, #Separate the distance between step
                                       steps_annotation_nudge = 0.1,
                                       steps_annotation_size  = 3,
                                       steps_values_annotate = T,
                                       steps_names = steps_names,
                                       grid_rows = "nbatch",
                                       grid_cols = "nolap",
                                       x_name = "ICC",
                                       y_name = "Bias",
                                       hline_intercept = 0,
                                       spu_x_shift = 2,
                                       legend_name = "Model",
                                       legend_labels = c("Random MA; Tau:SJ; CI:classic", "Random MA; Tau:SJ; CI:HK", "Random MA; Tau:SJ; CI:HK.SEadhoc", "Random MA; Tau:SJ; CI:HK.IQadhoc", "Random MA; Tau:SJ; CI:KR"),
                                       post_processing = list(
                                         add_custom_theme = list(
                                           legend.position="bottom",
                                           axis.text.x = element_text(angle = -90, 
                                                                      vjust = 0.5, 
                                                                      size = 4) 
                                         )))

bias_MA_shared_SJ = nested_loop_plot(resdf = all_contsresults_bias %>%
                                         mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                nolap = as.factor(as.character(nolap))) %>%
                                         filter(sharedtime != 1 & effsize == 0) %>%
                                         mutate(ICC = factor(ICC)),
                                       x = "ICC", 
                                       methods = c("Bias_MA.random.SJ.classic", "Bias_MA.random.SJ.HK", "Bias_MA.random.SJ.HK.SEadhoc", "Bias_MA.random.SJ.HK.IQadhoc", "Bias_MA.random.SJ.KR"),
                                       steps = steps,
                                       colors = brewer.pal(5, name = "Dark2"),
                                       steps_y_base = -0.03,
                                       steps_y_height = 0.005,
                                       steps_y_shift = 0.005, #Separate the distance between step
                                       steps_annotation_nudge = 0.1,
                                       steps_annotation_size  = 3,
                                       steps_values_annotate = T,
                                       steps_names = steps_names,
                                       grid_rows = "nbatch",
                                       grid_cols = "nolap",
                                       x_name = "ICC",
                                       y_name = "Bias",
                                      hline_intercept = 0,
                                     spu_x_shift = 2,
                                       legend_name = "Model",
                                       legend_labels = c("Random MA; Tau:SJ; CI:classic", "Random MA; Tau:SJ; CI:HK", "Random MA; Tau:SJ; CI:HK.SEadhoc", "Random MA; Tau:SJ; CI:HK.IQadhoc", "Random MA; Tau:SJ; CI:KR"),
                                       post_processing = list(
                                         add_custom_theme = list(
                                           legend.position="bottom",
                                           axis.text.x = element_text(angle = -90, 
                                                                      vjust = 0.5, 
                                                                      size = 4) 
                                         )))
#Saving the results.
bias_MA_separate_SJ
ggsave("bias_MA_separate_SJ.png", height = 13, width = 15)
bias_MA_shared_SJ
ggsave("bias_MA_shared_SJ.png", height = 13, width = 15)

bias_MA_separate_reml = nested_loop_plot(resdf = all_contsresults_bias %>%
                                           mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                  nolap = as.factor(as.character(nolap))) %>%
                                           filter(sharedtime == 1 & effsize == 0) %>%
                                           mutate(ICC = factor(ICC)),
                                         x = "ICC", 
                                         methods = c("Bias_MA.random.reml.classic", "Bias_MA.random.reml.HK", "Bias_MA.random.SJ.HK.SEadhoc", "Bias_MA.random.reml.HK.IQadhoc", "Bias_MA.random.reml.KR"),
                                         steps = steps,
                                         colors = brewer.pal(5, name = "Dark2"),
                                         steps_y_base = -0.017,
                                         steps_y_height = 0.005,
                                         steps_y_shift = 0.005, #Separate the distance between step
                                         steps_annotation_nudge = 0.1,
                                         steps_annotation_size  = 3,
                                         steps_values_annotate = T,
                                         steps_names = steps_names,
                                         grid_rows = "nbatch",
                                         grid_cols = "nolap",
                                         x_name = "ICC",
                                         y_name = "Bias",
                                         hline_intercept = 0,
                                         legend_name = "Model",
                                         spu_x_shift = 2,
                                         legend_labels = c("Random MA; Tau:REML; CI:classic", "Random MA; Tau:REML; CI:HK", "Random MA; Tau:REML; CI:HK.SEadhoc", "Random MA; Tau:REML; CI:HK.IQadhoc", "Random MA; Tau:REML; CI:KR"),
                                         post_processing = list(
                                           add_custom_theme = list(
                                             legend.position="bottom",
                                             axis.text.x = element_text(angle = -90, 
                                                                        vjust = 0.5, 
                                                                        size = 4) 
                                           )))

bias_MA_shared_reml = nested_loop_plot(resdf = all_contsresults_bias %>%
                                           mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                  nolap = as.factor(as.character(nolap))) %>%
                                           filter(sharedtime != 1 & effsize == 0) %>%
                                           mutate(ICC = factor(ICC)),
                                         x = "ICC", 
                                         methods = c("Bias_MA.random.reml.classic", "Bias_MA.random.reml.HK", "Bias_MA.random.SJ.HK.SEadhoc", "Bias_MA.random.reml.HK.IQadhoc", "Bias_MA.random.reml.KR"),
                                         steps = steps,
                                         colors = brewer.pal(5, name = "Dark2"),
                                         steps_y_base = -0.017,
                                         steps_y_height = 0.005,
                                         steps_y_shift = 0.005, #Separate the distance between step
                                         steps_annotation_nudge = 0.1,
                                         steps_annotation_size  = 3,
                                         steps_values_annotate = T,
                                       spu_x_shift = 2,
                                         steps_names = steps_names,
                                         grid_rows = "nbatch",
                                         grid_cols = "nolap",
                                         x_name = "ICC",
                                         y_name = "Bias",
                                         hline_intercept = 0,
                                         legend_name = "Model",
                                         legend_labels = c("Random MA; Tau:REML; CI:classic", "Random MA; Tau:REML; CI:HK", "Random MA; Tau:REML; CI:HK.SEadhoc", "Random MA; Tau:REML; CI:HK.IQadhoc", "Random MA; Tau:REML; CI:KR"),
                                         post_processing = list(
                                           add_custom_theme = list(
                                             legend.position="bottom",
                                             axis.text.x = element_text(angle = -90, 
                                                                        vjust = 0.5, 
                                                                        size = 4) 
                                           )))
#Saving the results
bias_MA_separate_reml
ggsave("bias_MA_separate_REML.png", height = 13, width = 15)
bias_MA_shared_reml
ggsave("bias_MA_shared_REML.png", height = 13, width = 15)

######Coverage######

coverage_MA_separate_SJ = nested_loop_plot(resdf = all_contsresults_coverage %>%
                                         mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                nolap = as.factor(as.character(nolap))) %>%
                                         filter(sharedtime == 1 & effsize == 0) %>%
                                         mutate(ICC = factor(ICC)),
                                       x = "ICC", 
                                       methods = c("Coverage_MA.random.SJ.classic", "Coverage_MA.random.SJ.HK", "Coverage_MA.random.SJ.HK.SEadhoc", "Coverage_MA.random.SJ.HK.IQadhoc", "Coverage_MA.random.SJ.KR"),
                                       steps = steps,
                                       colors = brewer.pal(5, name = "Dark2"),
                                      steps_y_base = 0.65,
                                       steps_y_height = 0.04,
                                       steps_y_shift = 0.04,
                                       spu_x_shift = 2,
                                       steps_annotation_nudge = 0.1,
                                       steps_annotation_size  = 3,
                                       steps_values_annotate = T,
                                       steps_names = steps_names,
                                       grid_rows = "nbatch",
                                       grid_cols = "nolap",
                                       x_name = "ICC",
                                       y_name = "Coverage",
                                       hline_intercept = 0.95,
                                       legend_name = "Model",
                                       legend_labels = c("Random MA; Tau:SJ; CI:classic", "Random MA; Tau:SJ; CI:HK", "Random MA; Tau:SJ; CI:HK.SEadhoc", "Random MA; Tau:SJ; CI:HK.IQadhoc", "Random MA; Tau:SJ; CI:KR"),
                                       point_alpha = 0.3,
                                       post_processing = list(
                                         add_custom_theme = list(
                                           legend.position="bottom",
                                           axis.text.x = element_text(angle = -90, 
                                                                      vjust = 0.5, 
                                                                      size = 4) 
                                         )))

coverage_MA_shared_SJ = nested_loop_plot(resdf = all_contsresults_coverage %>%
                                       mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                              nolap = as.factor(as.character(nolap))) %>%
                                       filter(sharedtime != 1 & effsize == 0) %>%
                                       mutate(ICC = factor(ICC)),
                                     x = "ICC", 
                                     methods = c("Coverage_MA.random.SJ.classic", "Coverage_MA.random.SJ.HK", "Coverage_MA.random.SJ.HK.SEadhoc", "Coverage_MA.random.SJ.HK.IQadhoc", "Coverage_MA.random.SJ.KR"),
                                     steps = steps,
                                     colors = brewer.pal(5, name = "Dark2"),
                                     steps_y_base = 0.65,
                                     steps_y_height = 0.04,
                                     steps_y_shift = 0.04,
                                     spu_x_shift = 2,
                                     steps_annotation_nudge = 0.1,
                                     steps_annotation_size  = 3,
                                     steps_values_annotate = T,
                                     steps_names = steps_names,
                                     grid_rows = "nbatch",
                                     grid_cols = "nolap",
                                     x_name = "ICC",
                                     y_name = "Coverage",
                                     hline_intercept = 0.95,
                                     legend_name = "Model",
                                     legend_labels = c("Random MA; Tau:SJ; CI:classic", "Random MA; Tau:SJ; CI:HK", "Random MA; Tau:SJ; CI:HK.SEadhoc", "Random MA; Tau:SJ; CI:HK.IQadhoc", "Random MA; Tau:SJ; CI:KR"),
                                     point_alpha = 0.3,
                                     post_processing = list(
                                       add_custom_theme = list(
                                         legend.position="bottom",
                                         axis.text.x = element_text(angle = -90, 
                                                                    vjust = 0.5, 
                                                                    size = 4) 
                                       )))

#Saving the results
coverage_MA_separate_SJ
ggsave("coverage_MA_separate_SJ.png", height = 13, width = 15)
coverage_MA_shared_SJ
ggsave("coverage_MA_shared_SJ.png", height = 13, width = 15)

coverage_MA_separate_reml = nested_loop_plot(resdf = all_contsresults_coverage %>%
                                           mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                  nolap = as.factor(as.character(nolap))) %>%
                                           filter(sharedtime == 1 & effsize == 0) %>%
                                           mutate(ICC = factor(ICC)),
                                         x = "ICC", 
                                         methods = c("Coverage_MA.random.reml.classic", "Coverage_MA.random.reml.HK", "Coverage_MA.random.reml.HK.SEadhoc", "Coverage_MA.random.reml.HK.IQadhoc", "Coverage_MA.random.reml.KR"),
                                         steps = steps,
                                         colors = brewer.pal(5, name = "Dark2"),
                                         steps_y_base = 0.65,
                                         steps_y_height = 0.04,
                                         steps_y_shift = 0.04,
                                         spu_x_shift = 2,
                                         steps_annotation_nudge = 0.1,
                                         steps_annotation_size  = 3,
                                         steps_values_annotate = T,
                                         steps_names = steps_names,
                                         grid_rows = "nbatch",
                                         grid_cols = "nolap",
                                         x_name = "ICC",
                                         y_name = "Coverage",
                                         hline_intercept = 0.95,
                                         legend_name = "Model",
                                         legend_labels = c("Random MA; Tau:REML; CI:classic", "Random MA; Tau:REML; CI:HK", "Random MA; Tau:REML; CI:HK.SEadhoc", "Random MA; Tau:REML; CI:HK.IQadhoc", "Random MA; Tau:REML; CI:KR"),
                                         point_alpha = 0.3,
                                         post_processing = list(
                                           add_custom_theme = list(
                                             legend.position="bottom",
                                             axis.text.x = element_text(angle = -90, 
                                                                        vjust = 0.5, 
                                                                        size = 4) 
                                           )))

coverage_MA_shared_reml = nested_loop_plot(resdf = all_contsresults_coverage %>%
                                         mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                nolap = as.factor(as.character(nolap))) %>%
                                         filter(sharedtime != 1 & effsize == 0) %>%
                                         mutate(ICC = factor(ICC)),
                                       x = "ICC", 
                                       methods = c("Coverage_MA.random.reml.classic", "Coverage_MA.random.reml.HK", "Coverage_MA.random.reml.HK.SEadhoc", "Coverage_MA.random.reml.HK.IQadhoc", "Coverage_MA.random.reml.KR"),
                                       steps = steps,
                                       colors = brewer.pal(5, name = "Dark2"),
                                       steps_y_base = 0.65,
                                       steps_y_height = 0.04,
                                       steps_y_shift = 0.04,
                                       steps_annotation_nudge = 0.1,
                                       steps_annotation_size  = 3,
                                       steps_values_annotate = T,
                                       spu_x_shift = 2,
                                       steps_names = steps_names,
                                       grid_rows = "nbatch",
                                       grid_cols = "nolap",
                                       x_name = "ICC",
                                       y_name = "Coverage",
                                       hline_intercept = 0.95,
                                       legend_name = "Model",
                                       legend_labels = c("Random MA; Tau:REML; CI:classic", "Random MA; Tau:REML; CI:HK", "Random MA; Tau:REML; CI:HK.SEadhoc", "Random MA; Tau:REML; CI:HK.IQadhoc", "Random MA; Tau:REML; CI:KR"),
                                       point_alpha = 0.3,
                                       post_processing = list(
                                         add_custom_theme = list(
                                           legend.position="bottom",
                                           axis.text.x = element_text(angle = -90, 
                                                                      vjust = 0.5, 
                                                                      size = 4) 
                                         )))

#Saving the results
coverage_MA_separate_reml
ggsave("coverage_MA_separate_REML.png", height = 13, width = 15)
coverage_MA_shared_reml
ggsave("coverage_MA_shared_REML.png", height = 13, width = 15)

######Empirical standard error######

EmpSE_MA_separate_SJ = nested_loop_plot(resdf = all_contsresults_EmpSE %>%
                                             mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                    nolap = as.factor(as.character(nolap))) %>%
                                             filter(sharedtime == 1 & effsize == 0) %>%
                                             mutate(ICC = factor(ICC)),
                                           x = "ICC", 
                                           methods = c("Emp_SE_MA.random.SJ.classic", "Emp_SE_MA.random.SJ.HK", "Emp_SE_MA.random.SJ.HK.SEadhoc", "Emp_SE_MA.random.SJ.HK.IQadhoc", "Emp_SE_MA.random.SJ.KR"),
                                           steps = steps,
                                           colors = brewer.pal(5, name = "Dark2"),
                                           steps_y_base = -0.014,
                                        steps_y_height = 0.01,
                                        steps_y_shift = 0.03,
                                        spu_x_shift = 2,
                                           steps_annotation_nudge = 0.1,
                                           steps_annotation_size  = 3,
                                           steps_values_annotate = T,
                                           steps_names = steps_names,
                                           grid_rows = "nbatch",
                                           grid_cols = "nolap",
                                           x_name = "ICC",
                                           y_name = "Empirical standard error",
                                           hline_intercept = 0,
                                           legend_name = "Model",
                                           legend_labels = c("Random MA; Tau:SJ; CI:classic", "Random MA; Tau:SJ; CI:HK", "Random MA; Tau:SJ; CI:HK.SEadhoc", "Random MA; Tau:SJ; CI:HK.IQadhoc", "Random MA; Tau:SJ; CI:KR"),
                                           post_processing = list(
                                             add_custom_theme = list(
                                               legend.position="bottom",
                                               axis.text.x = element_text(angle = -90, 
                                                                          vjust = 0.5, 
                                                                          size = 4) 
                                             )))

EmpSE_MA_shared_SJ = nested_loop_plot(resdf = all_contsresults_EmpSE %>%
                                           mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                  nolap = as.factor(as.character(nolap))) %>%
                                           filter(sharedtime != 1 & effsize == 0) %>%
                                           mutate(ICC = factor(ICC)),
                                         x = "ICC", 
                                      methods = c("Emp_SE_MA.random.SJ.classic", "Emp_SE_MA.random.SJ.HK", "Emp_SE_MA.random.SJ.HK.SEadhoc", "Emp_SE_MA.random.SJ.HK.IQadhoc", "Emp_SE_MA.random.SJ.KR"),
                                         steps = steps,
                                         colors = brewer.pal(5, name = "Dark2"),
                                          steps_y_base = -0.014,
                                      steps_y_height = 0.01,
                                      steps_y_shift = 0.03,
                                      spu_x_shift = 2,
                                         steps_annotation_nudge = 0.1,
                                         steps_annotation_size  = 3,
                                         steps_values_annotate = T,
                                         steps_names = steps_names,
                                         grid_rows = "nbatch",
                                         grid_cols = "nolap",
                                         x_name = "ICC",
                                         y_name = "Empirical standard error",
                                         hline_intercept = 0,
                                         legend_name = "Model",
                                         legend_labels = c("Random MA; Tau:SJ; CI:classic", "Random MA; Tau:SJ; CI:HK", "Random MA; Tau:SJ; CI:HK.SEadhoc", "Random MA; Tau:SJ; CI:HK.IQadhoc", "Random MA; Tau:SJ; CI:KR"),
                                         post_processing = list(
                                           add_custom_theme = list(
                                             legend.position="bottom",
                                             axis.text.x = element_text(angle = -90, 
                                                                        vjust = 0.5, 
                                                                        size = 4) 
                                           )))

#Saving the results
EmpSE_MA_separate_SJ
ggsave("EmpSE_MA_separate_SJ.png", height = 13, width = 15)
EmpSE_MA_shared_SJ
ggsave("EmpSE_MA_shared_SJ.png", height = 13, width = 15)

EmpSE_MA_separate_reml = nested_loop_plot(resdf = all_contsresults_EmpSE %>%
                                               mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                      nolap = as.factor(as.character(nolap))) %>%
                                               filter(sharedtime == 1 & effsize == 0) %>%
                                               mutate(ICC = factor(ICC)),
                                             x = "ICC", 
                                             methods = c("Emp_SE_MA.random.reml.classic", "Emp_SE_MA.random.reml.HK", "Emp_SE_MA.random.reml.HK.SEadhoc", "Emp_SE_MA.random.reml.HK.IQadhoc", "Emp_SE_MA.random.reml.KR"),
                                             steps = steps,
                                             colors = brewer.pal(5, name = "Dark2"),
                                             steps_y_base = -0.014,
                                          steps_y_height = 0.01,
                                          steps_y_shift = 0.03,
                                          spu_x_shift = 2,
                                             steps_annotation_nudge = 0.1,
                                             steps_annotation_size  = 3,
                                             steps_values_annotate = T,
                                             steps_names = steps_names,
                                             grid_rows = "nbatch",
                                             grid_cols = "nolap",
                                             x_name = "ICC",
                                             y_name = "Empirical standard error",
                                             hline_intercept = 0,
                                             legend_name = "Model",
                                             legend_labels = c("Random MA; Tau:REML; CI:classic", "Random MA; Tau:REML; CI:HK", "Random MA; Tau:REML; CI:HK.SEadhoc", "Random MA; Tau:REML; CI:HK.IQadhoc", "Random MA; Tau:REML; CI:KR"),
                                             post_processing = list(
                                               add_custom_theme = list(
                                                 legend.position="bottom",
                                                 axis.text.x = element_text(angle = -90, 
                                                                            vjust = 0.5, 
                                                                            size = 4) 
                                               )))

EmpSE_MA_shared_reml = nested_loop_plot(resdf = all_contsresults_EmpSE %>%
                                             mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                    nolap = as.factor(as.character(nolap))) %>%
                                             filter(sharedtime != 1 & effsize == 0) %>%
                                             mutate(ICC = factor(ICC)),
                                           x = "ICC", 
                                        methods = c("Emp_SE_MA.random.reml.classic", "Emp_SE_MA.random.reml.HK", "Emp_SE_MA.random.reml.HK.SEadhoc", "Emp_SE_MA.random.reml.HK.IQadhoc", "Emp_SE_MA.random.reml.KR"),
                                           steps = steps,
                                           colors = brewer.pal(5, name = "Dark2"),
                                           steps_y_base = -0.014,
                                        steps_y_height = 0.01,
                                        steps_y_shift = 0.03,
                                        spu_x_shift = 2,
                                           steps_annotation_nudge = 0.1,
                                           steps_annotation_size  = 3,
                                           steps_values_annotate = T,
                                           steps_names = steps_names,
                                           grid_rows = "nbatch",
                                           grid_cols = "nolap",
                                           x_name = "ICC",
                                           y_name = "Empirical standard error",
                                           hline_intercept = 0,
                                           legend_name = "Model",
                                           legend_labels = c("Random MA; Tau:REML; CI:classic", "Random MA; Tau:REML; CI:HK", "Random MA; Tau:REML; CI:HK.SEadhoc", "Random MA; Tau:REML; CI:HK.IQadhoc", "Random MA; Tau:REML; CI:KR"),
                                           post_processing = list(
                                             add_custom_theme = list(
                                               legend.position="bottom",
                                               axis.text.x = element_text(angle = -90, 
                                                                          vjust = 0.5, 
                                                                          size = 4) 
                                             )))
#Saving the results
EmpSE_MA_separate_reml
ggsave("EmpSE_MA_separate_REML.png", height = 13, width = 15)
EmpSE_MA_shared_reml
ggsave("EmpSE_MA_shared_REML.png", height = 13, width = 15)

########################LINEAR-MIXED MODELS vs META-ANALYSIS PLOTS########################

#Prespecifying steps and legend labels for the nested loop plots.
steps = c("treateffvar", "M", "K",  "T", "CAC")
steps_names = c(expression(sigma[eta]^2), "m", "Clusters", "Period", "CAC")
legend_labels = c("Model A", "Model B", "Model C", "Model D", "Fixed-effects meta-analysis", "Random MA; Tau:REML; CI:HK")

######Bias######

#Bias for separate period effects
bias_final_separate = nested_loop_plot(resdf = all_contsresults_bias %>%
                                         mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                nolap = as.factor(as.character(nolap))) %>%
                                         filter(sharedtime == 1 & effsize == 0) %>%
                                         mutate(ICC = factor(ICC)),
                                       x = "ICC", 
                                       methods = c("Bias_mod_A", "Bias_mod_B", "Bias_mod_C", "Bias_mod_D", "Bias_MA.fixed", "Bias_MA.random.reml.HK"),
                                       steps = steps,
                                       colors = brewer.pal(6, name = "Dark2"),
                                       steps_y_base = -1,
                                       spu_x_shift = 2,
                                       steps_y_height = 0.2,
                                       steps_y_shift = 0.2, #Separate the distance between step
                                       steps_annotation_nudge = 0.1,
                                       steps_annotation_size  = 3,
                                       steps_values_annotate = T,
                                       steps_names = steps_names,
                                       grid_rows = "nbatch",
                                       grid_cols = "nolap",
                                       x_name = "ICC",
                                       y_name = "Bias",
                                       hline_intercept = 0,
                                       legend_name = "Model",
                                       legend_labels = legend_labels,
                                       post_processing = list(
                                         add_custom_theme = list(
                                           legend.position="bottom",
                                           axis.text.x = element_text(angle = -90, 
                                                                      vjust = 0.5, 
                                                                      size = 4) 
                                         )))

#Bias for shared period effects
bias_final_shared = nested_loop_plot(resdf = all_contsresults_bias %>%
                                       mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                              nolap = as.factor(as.character(nolap))) %>%
                                       filter(sharedtime != 1 & effsize == 0) %>%
                                       mutate(ICC = factor(ICC)),
                                     x = "ICC", 
                                     methods = c("Bias_mod_A", "Bias_mod_B", "Bias_mod_C", "Bias_mod_D", "Bias_MA.fixed", "Bias_MA.random.reml.HK"),
                                     steps = steps,
                                     colors = brewer.pal(6, name = "Dark2"),
                                     line_linetypes = c(2, 2, 2, 2, 1, 1),
                                     point_shapes = c(15,15,15,15,17,17,17, 17, 17, 17),
                                     spu_x_shift = 2,
                                     steps_y_base = -0.027,
                                     steps_y_height = 0.005,
                                     steps_y_shift = 0.005, #Separate the distance between step
                                     steps_annotation_nudge = 0.1,
                                     steps_annotation_size  = 3,
                                     steps_values_annotate = T,
                                     steps_names = steps_names,
                                     grid_rows = "nbatch",
                                     grid_cols = "nolap",
                                     x_name = "ICC",
                                     y_name = "Bias",
                                     hline_intercept = 0,
                                     legend_name = "Model",
                                     legend_labels = legend_labels,
                                     post_processing = list(
                                       add_custom_theme = list(
                                         legend.position="bottom",
                                         axis.text.x = element_text(angle = -90, 
                                                                    vjust = 0.5, 
                                                                    size = 4) 
                                       )))

#Bias for separate period effects without models A and C
bias_final_separate_without = nested_loop_plot(resdf = all_contsresults_bias %>%
                                                 mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                        nolap = as.factor(as.character(nolap))) %>%
                                                 filter(sharedtime == 1 & effsize == 0) %>%
                                                 mutate(ICC = factor(ICC)),
                                               x = "ICC", 
                                               methods = c("Bias_mod_B", "Bias_mod_D", "Bias_MA.fixed", "Bias_MA.random.reml.HK"),
                                               steps = steps,
                                               colors = brewer.pal(4, name = "Dark2"),
                                               line_linetypes = c(1,1,2, 2),
                                               steps_y_base = -0.019,
                                               steps_y_height = 0.005,
                                               steps_y_shift = 0.005, #Separate the distance between step
                                               steps_annotation_nudge = 0.08,
                                               steps_annotation_size  = 2.5,
                                               steps_values_annotate = T,
                                               steps_names = steps_names,
                                               spu_x_shift = 2,
                                               grid_rows = "nbatch",
                                               grid_cols = "nolap",
                                               x_name = "ICC",
                                               y_name = "Bias",
                                               hline_intercept = 0,
                                               legend_name = "Model",
                                               legend_labels = c("Model B", "Model D", "Fixed-effects meta-analysis", "Random MA; Tau:REML; CI:HK"),
                                               post_processing = list(
                                                 add_custom_theme = list(
                                                   legend.position="bottom",
                                                   axis.text.x = element_text(angle = -90, 
                                                                              vjust = 0.5, 
                                                                              size = 4) 
                                                 )))

#Saving the results
bias_final_separate
ggsave("bias_final_separate.png", height = 13, width = 15)

bias_final_shared
ggsave("bias_final_shared.png", height = 13, width = 15)

bias_final_separate_without
ggsave("bias_final_separate_without_modelA&C.png", height = 13, width = 15)


######Coverage######

#Coverage for separate period effects
coverage_final_separate = nested_loop_plot(resdf = all_contsresults_coverage %>%
                                             mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                    nolap = as.factor(as.character(nolap))) %>%
                                             filter(sharedtime == 1 & effsize == 0) %>%
                                             mutate(ICC = factor(ICC)),
                                           x = "ICC", 
                                           methods = c("Coverage_mod_B", "Coverage_mod_D", "Coverage_MA.fixed", "Coverage_MA.random.reml.HK"),
                                           steps = steps,
                                           colors = brewer.pal(6, name = "Dark2"),
                                           line_linetypes = c(1, 1, 2, 2),
                                           point_shapes = c(15, 15, 17, 17),
                                           steps_y_base = 0.15,
                                           steps_y_height = 0.08,
                                           steps_y_shift = 0.08,
                                           steps_annotation_nudge = 0.1,
                                           steps_annotation_size  = 3,
                                           steps_values_annotate = T,
                                           spu_x_shift = 2,
                                           steps_names = steps_names,
                                           grid_rows = "nbatch",
                                           grid_cols = "nolap",
                                           x_name = "ICC",
                                           y_name = "Coverage",
                                           hline_intercept = 0.95,
                                           legend_name = "Model",
                                           legend_labels = c("Model B", "Model D", "Fixed-effects meta-analysis", "Random MA; Tau:REML; CI:HK"),
                                           post_processing = list(
                                             add_custom_theme = list(
                                               legend.position="bottom",
                                               axis.text.x = element_text(angle = -90, 
                                                                          vjust = 0.5, 
                                                                          size = 4) 
                                             )))

#Coverage for shared period effects
coverage_final_shared = nested_loop_plot(resdf = all_contsresults_coverage %>%
                                           mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                  nolap = as.factor(as.character(nolap))) %>%
                                           filter(sharedtime != 1 & effsize == 0) %>%
                                           mutate(ICC = factor(ICC)),
                                         x = "ICC", 
                                         methods = c("Coverage_mod_A", "Coverage_mod_B", "Coverage_mod_C", "Coverage_mod_D", "Coverage_MA.fixed", "Coverage_MA.random.reml.HK"),
                                         steps = steps,
                                         colors = brewer.pal(6, name = "Dark2"),
                                         line_linetypes = c(2, 2, 2, 2, 1, 1),
                                         point_shapes = c(15, 15, 15, 15, 17, 17),
                                         steps_y_base = 0.15,
                                         steps_y_height = 0.08,
                                         steps_y_shift = 0.08,
                                         steps_annotation_nudge = 0.1,
                                         spu_x_shift = 2,
                                         steps_annotation_size  = 3,
                                         steps_values_annotate = T,
                                         steps_names = steps_names,
                                         grid_rows = "nbatch",
                                         grid_cols = "nolap",
                                         x_name = "ICC",
                                         y_name = "Coverage",
                                         hline_intercept = 0.95,
                                         legend_name = "Model",
                                         legend_labels = legend_labels,
                                         post_processing = list(
                                           add_custom_theme = list(
                                             legend.position="bottom",
                                             axis.text.x = element_text(angle = -90, 
                                                                        vjust = 0.5, 
                                                                        size = 4) 
                                           )))

#Saving the results
coverage_final_separate
ggsave("coverage_final_separate.png", height = 13, width = 15)

coverage_final_shared
ggsave("coverage_final_shared.png", height = 13, width = 15)


######Empirical standard error######

#Empirical standard error for separate period effects
EmpSE_final_separate = nested_loop_plot(resdf = all_contsresults_EmpSE %>%
                                          mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                 nolap = as.factor(as.character(nolap))) %>%
                                          filter(sharedtime == 1 & effsize == 0) %>%
                                          mutate(ICC = factor(ICC)),
                                        x = "ICC", 
                                        methods = c("Emp_SE_mod_B",  "Emp_SE_mod_D", "Emp_SE_MA.fixed", "Emp_SE_MA.random.reml.HK"),
                                        steps = steps,
                                        colors = brewer.pal(6, name = "Dark2"),
                                        line_linetypes = c(2,2,1, 1),
                                        point_shapes = c(15, 15, 17, 17),
                                        steps_y_base = -0.03,
                                        steps_y_height = 0.03,
                                        steps_y_shift = 0.03,
                                        steps_annotation_nudge = 0.1,
                                        steps_annotation_size  = 3,
                                        steps_values_annotate = T,
                                        spu_x_shift = 2,
                                        steps_names = steps_names,
                                        grid_rows = "nbatch",
                                        grid_cols = "nolap",
                                        x_name = "ICC",
                                        y_name = "Empirical standard error",
                                        hline_intercept = 0,
                                        legend_name = "Model",
                                        legend_labels = c("Model B", "Model D", "Fixed-effects meta-analysis", "Random MA; Tau:REML; CI:HK"),
                                        post_processing = list(
                                          add_custom_theme = list(
                                            legend.position="bottom",
                                            axis.text.x = element_text(angle = -90, 
                                                                       vjust = 0.5, 
                                                                       size = 4) 
                                          )))

#Empirical standard error for shared period effects
EmpSE_final_shared = nested_loop_plot(resdf = all_contsresults_EmpSE %>%
                                        mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                               nolap = as.factor(as.character(nolap))) %>%
                                        filter(sharedtime != 1 & effsize == 0) %>%
                                        mutate(ICC = factor(ICC)),
                                      x = "ICC", 
                                      methods = c("Emp_SE_mod_A", "Emp_SE_mod_B", "Emp_SE_mod_C", "Emp_SE_mod_D", "Emp_SE_MA.fixed", "Emp_SE_MA.random.reml.HK"),
                                      steps = steps,
                                      colors = brewer.pal(6, name = "Dark2"),
                                      line_linetypes = c(2, 2, 2, 2, 1, 1),
                                      point_shapes = c(15, 15, 15, 15, 17, 17),
                                      steps_y_base = -0.03,
                                      steps_y_height = 0.03,
                                      steps_y_shift = 0.03,
                                      steps_annotation_nudge = 0.1,
                                      steps_annotation_size  = 3,
                                      steps_values_annotate = T,
                                      spu_x_shift = 2,
                                      steps_names = steps_names,
                                      grid_rows = "nbatch",
                                      grid_cols = "nolap",
                                      x_name = "ICC",
                                      y_name = "Empirical standard error",
                                      hline_intercept = 0,
                                      legend_name = "Model",
                                      legend_labels = legend_labels,
                                      post_processing = list(
                                        add_custom_theme = list(
                                          legend.position="bottom",
                                          axis.text.x = element_text(angle = -90, 
                                                                     vjust = 0.5, 
                                                                     size = 4) 
                                        )))

EmpSE_final_separate
ggsave("EmpSE_final_separate.png", height = 13, width = 15)

EmpSE_final_shared
ggsave("EmpSE_final_shared.png", height = 13, width = 15)

######Rate of non-convergence######

#Rate of non-convergence for separate period effects
nonconv_separate = nested_loop_plot(resdf = all_contsresults_nonconv %>%
                                          mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                                 nolap = as.factor(as.character(nolap))) %>%
                                          filter(sharedtime == 1 & effsize == 0) %>%
                                          mutate(ICC = factor(ICC)),
                                        x = "ICC", 
                                        methods = c("mod_A", "mod_B", "mod_C", "mod_D", "MA.random.reml.HK", "MA.fixed"),
                                        steps = steps,
                                        colors = brewer.pal(6, name = "Dark2"),
                                    line_linetypes = c(2, 2, 2, 2, 1, 1),
                                    point_shapes = c(15, 15, 15, 15, 17, 17),
                                        steps_y_base = -0.05,
                                        steps_y_height = 0.05,
                                        steps_y_shift = 0.05,
                                        steps_annotation_nudge = 0.1,
                                        steps_annotation_size  = 3,
                                        steps_values_annotate = T,
                                        spu_x_shift = 2,
                                        steps_names = steps_names,
                                        grid_rows = "nbatch",
                                        grid_cols = "nolap",
                                        x_name = "ICC",
                                        y_name = "Rate of non-convergence",
                                        hline_intercept = 0,
                                        legend_name = "Model",
                                        legend_labels = c("Model A", "Model B", "Model C", "Model D", "Random MA; Tau:REML; CI:HK", "Fixed-effects meta-analysis"),
                                        post_processing = list(
                                          add_custom_theme = list(
                                            legend.position="bottom",
                                            axis.text.x = element_text(angle = -90, 
                                                                       vjust = 0.5, 
                                                                       size = 4) 
                                          )))

#Rate of non-convergence for shared period effects
nonconv_shared = nested_loop_plot(resdf = all_contsresults_nonconv %>%
                                      mutate(nolap = ifelse(nolap == 2 | nolap == 4, "T-2", nolap),
                                             nolap = as.factor(as.character(nolap))) %>%
                                      filter(sharedtime == 0 & effsize == 0) %>%
                                      mutate(ICC = factor(ICC)),
                                    x = "ICC", 
                                    methods = c("mod_A", "mod_B", "mod_C", "mod_D", "MA.random.reml.HK", "MA.fixed"),
                                    steps = steps,
                                    colors = brewer.pal(6, name = "Dark2"),
                                    line_linetypes = c(2, 2, 2, 2, 1, 1),
                                    point_shapes = c(15, 15, 15, 15, 17, 17),
                                    steps_y_base = -0.05,
                                    steps_y_height = 0.05,
                                    steps_y_shift = 0.05,
                                    steps_annotation_nudge = 0.1,
                                    steps_annotation_size  = 3,
                                    steps_values_annotate = T,
                                    spu_x_shift = 2,
                                    steps_names = steps_names,
                                    grid_rows = "nbatch",
                                    grid_cols = "nolap",
                                    x_name = "ICC",
                                    y_name = "Rate of non-convergence",
                                    hline_intercept = 0,
                                    legend_name = "Model",
                                    legend_labels = c("Model A", "Model B", "Model C", "Model D", "Random MA; Tau:REML; CI:HK", "Fixed-effects meta-analysis"),
                                    post_processing = list(
                                      add_custom_theme = list(
                                        legend.position="bottom",
                                        axis.text.x = element_text(angle = -90, 
                                                                   vjust = 0.5, 
                                                                   size = 4) 
                                      )))

#Saving the results
nonconv_separate
ggsave("nonconv_overall_separate.png", height = 13, width = 15)

nonconv_shared
ggsave("nonconv_overall_shared.png", height = 13, width = 15)


