rm(list=ls())

#### Libraries ####
library(dplyr) # easier data wrangling 
library(randomForest) # randomForest()
library(caret) # train()
library(gbm) # gbm()
library(kernlab) # kernelMatrix()
library(personalized) # weighted.ksvm()
library(ggplot2) # ggplot

#### Set Directories ####
jobname = "TRule_simdata_ITRpred"
rseed = 72; rsd = paste0("seed",rseed)
fitseed = 2020 # random seed used to fit machine learning models
cs = paste0("case",0)
grp = 3; group = paste0("group",grp) # subgroup index in which learn ITR
casename = paste(jobname,cs,rsd,group,sep="_")

maindir = "/pine/scr/j/i/jitong/thesis/TRule_model" # project directory
# maindir = "D:/bios994/TRule_model" # project directory
codesdir = file.path(maindir,"codes") # code directory
wsdir = file.path(maindir,"rdata") # R workspace directory

#### Load data ####
source(file=file.path(codesdir, "functions_ITRlearn.R"))
load(file=file.path(wsdir,paste("TRule_simdata_ITRlearn",cs,rsd,group,"workspace.RData",sep="_")))
load(file=file.path(wsdir,paste("TRule_simdata_ITRlearn",cs,rsd,group,"est.RData",sep="_")))

#### Simulate another dataset ####
N_pred = 500
set.seed(fitseed)

covariate_pred = mapply(
    FUN = rnorm, n = N_pred, 
    mean = apply(dat_incluster[,c(demo_name,bm_name)], 2, mean),
    sd = apply(dat_incluster[,c(demo_name,bm_name)], 2, sd)
  ) %>%
  as.data.frame
colnames(covariate_pred) = c(demo_name,bm_name)

prop_pred = custom$modelInfo$prob(custom, newdata=covariate_pred)
prop_pred_adj = prop_pred %>%
  apply(1, function(row){
    (row - apply(pat_prop_trunc_norm, 2, mean))/apply(pat_prop_trunc_norm, 2, sd) 
  }) %>%
  t() %>%
  as.data.frame()
prop_pred_adj = prop_pred_adj %>%
  dplyr::select(-starts_with(trt_ref)) # remove the propensity score for the reference
colnames(prop_pred_adj) = paste0(colnames(prop_pred_adj), "_prop")

prog_pred = predict(prognostic.gbm, newdata=covariate_pred, n.trees=5000, interaction.depth=4)
prog_pred_adj = (prog_pred - mean(prognostic.gbm$fit))/sd(prognostic.gbm$fit) # do standardization later
  
dat_pred = data.frame(
    covariate_pred,
    prop_pred_adj,
    prog = prog_pred_adj
  ) %>%
  as.matrix # columns are all covariates, propensity scores, prognostic scores

ITR_models = ITR_OVO[["fits"]][[1]][["best_fit"]][["model"]] # a list of K(K-1)/2 binary classifiers

vote = predict.ITR(object = ITR_models, newdata = dat_pred, treatments = trts)

pred_info_ITR = data.frame(
    ID = paste0("test_",1:N_pred),
    group = paste("Group", grp),
    treatment = sample(trts, N_pred, replace = TRUE),
    covariate_pred,
    eps = rnorm(N_pred, 0, 1/3),
    vote = vote
  ) %>%
  mutate(lc = 1 + age - sex + 2*V1 - 2*V2) %>%
  mutate(reward = lc + eps)

df_evf_pred = map_propsc(dat_pred=pred_info_ITR, propsc=prop_pred) %>%
  dplyr::select(ID, group, reward, treatment, vote, pi)

df_evf_summary_TRT = df_evf_pred %>% 
  group_by(treatment) %>% 
  summarise(evf_ipw=sum(reward/pi)/sum(1/pi)) %>%
  as.data.frame()
df_evf_summary_TRT

attach(df_evf_pred)
ITR_evf_ipw = sum((reward/pi)[treatment==vote])/sum((1/pi)[treatment==vote])
detach(df_evf_pred)
ITR_evf_ipw

#### Barplot of proportion of treatments in ITR and TRT #### 
df_evf_dist_TRT = df_evf_pred %>%
  count(group, treatment) %>%
  group_by(group) %>%
  summarise(prop=n/sum(n), drug=treatment, type="Observed Treatments")

df_evf_dist_ITR = df_evf_pred %>%
  count(group, vote) %>%
  group_by(group) %>%
  summarise(prop=n/sum(n), drug=vote, type="Optimal ITRs")

df_evf_dist = rbind(df_evf_dist_TRT, df_evf_dist_ITR)

## draw figures
windowsFonts(Times=windowsFont("Times New Roman"))

plot_TRT_ITR_dist = ggplot(df_evf_dist, aes(x=drug, y=prop, fill=type)) +
  geom_bar(stat="identity", width=0.5, position=position_dodge()) +
  theme_bw() +
  facet_wrap(~group, scales="fixed",nrow=1) +
  scale_y_continuous(labels = scales::percent, 
                     expand = c(0,0), limits = c(0,max(df_evf_dist$prop)+0.1)) +
  scale_fill_manual(values = c("red","blue")) +
  labs(y = "Proportion", x="Treatment class", fill="",
       caption = paste0("EVF of ITRs: ", round(ITR_evf_ipw,3), ". ",
                        "EVF of universal rules: ", paste(trts, round(df_evf_summary_TRT$evf_ipw,3), sep=": ", collapse=", "), ".")
  ) +
  theme(legend.position = "top",
        legend.title = element_text(size = 12, family = "Times"),
        legend.text = element_text(size = 12, family = "Times"),
        plot.caption = element_text(size = 8, family = "Times")) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   size = 12, hjust = 1, family = "Times"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12, family = "Times"),
        axis.title.y = element_text(size = 12 , family = "Times"), # formats of y axis title
        # axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12, family = "Times"), # formats of y axis label
        strip.text.x = element_text(size = 12, colour = "black", 
                                    angle = 0, family = "Times") # formats of facet grid title
  )

ggsave(filename="TRule_simdata_ITRpredict_distribution.pdf", path=maindir, 
       plot=plot_TRT_ITR_dist,
       width = 8.5, height = 5.3, units = "in", dpi=72)

  


