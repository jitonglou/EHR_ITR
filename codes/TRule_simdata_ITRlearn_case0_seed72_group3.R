rm(list=ls())

#### Libraries ####
library(dplyr) # easier data wrangling 
library(randomForest) # randomForest()
library(caret) # train()
library(gbm) # gbm()
library(kernlab) # kernelMatrix()
library(personalized) # weighted.ksvm()

#### Set Directories ####
jobname = "TRule_simdata_ITRlearn"
rseed = 72; rsd = paste0("seed",rseed)
fitseed = 2020 # random seed used to fit machine learning models
cs = paste0("case",0)
grp = 3; group = paste0("group",grp) # subgroup index in which learn ITR
casename = paste(jobname,cs,rsd,group,sep="_")

maindir = "/pine/scr/j/i/jitong/thesis/TRule_model" # project directory
codesdir = file.path(maindir,"codes") # code directory
wsdir = file.path(maindir,"rdata") # R workspace directory

#### Load data ####
source(file=file.path(codesdir, "functions_ITRlearn.R"))
load(file=file.path(wsdir,paste0("TRule_simdata_ITRpreprocess_",cs,"_",rsd,"_df.RData")))

#### Functions and settings ####
trts = c("A","B","C") # treatment types
trt_ref = "A" # reference treatment in the estimation of propensity scores
g_func = function(x){abs(x)} # function used in equation (6) in our manuscript
SNN = 1 # size of the matched set (Size of Near Neighbor)
p = 2 # number of health markers

#### Calculate propensity score (within the current subgroup) ####
colnames(pat_info_ITR)
demo_name = c("age","sex") # names of demographic variables
bm_name = paste0("V", 1:p) # names of variables reflecting the recent pattern of health markers

## Extract related variables and recreate a data frame
dat_feature = pat_info_ITR %>%
  filter(cluster == grp & treatment %in% trts) %>%
  dplyr::select(ID, treatment, reward, all_of(demo_name), all_of(bm_name))

dat_feature$treatment = relevel(as.factor(dat_feature$treatment), ref=trt_ref)
dat_prop = dat_feature %>% dplyr::select(-ID, -reward)

## Train random forest models to estimate propensity scores
## Tuning parameters: 1. number of features selected at each node (.mtry). 2. number of trees (.ntree) 
## Use 3 repeats of 5-fold cross-validations to select best tuning parameters 
metric = "Accuracy"
tunegrid = expand.grid(.mtry=seq(2,ncol(dat_prop)-1, 1), .ntree=seq(500,2500,500))
control = trainControl(method="repeatedcv", number=5, repeats=3)
set.seed(fitseed)
custom = train(treatment~., data=dat_prop, 
               method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
print(custom) # cross-validation errors
custom$bestTune # best tuning parameters

table(pred=custom$modelInfo$predict(custom, newdata=dat_prop), 
      true=dat_prop$treatment) # classification results using the whole dataset

## propensity scores (standardized)
pat_prop = data.frame(ID=dat_feature$ID, custom$modelInfo$prob(custom, newdata=dat_prop))
prop_qt = quantile(c(as.matrix(pat_prop[,-1])), probs=c(0.01,0.025,0.05,seq(0.1,0.9,0.1),0.95,0.975,0.99))
prop_qt # display quantiles of propensity scores to see whether truncation is needed

pat_prop_trunc = pat_prop # no extreme propensity scores, do not truncate
pat_prop_trunc_norm = pat_prop_trunc[,-1] # remove column ID

pat_prop_adj = apply(pat_prop_trunc_norm, 2, function(x){
  return(scale(x))
}) %>% as.data.frame # scale propensity scores to estimate reward hat and calculate d(Hi,Hj)
colnames(pat_prop_adj) = paste0(colnames(pat_prop[,-1]), "_prop")
pat_prop_adj = pat_prop_adj %>%
  dplyr::select(-starts_with(paste0(trt_ref,"_prop"))) # remove the propensity score for the reference


#### Calculate prognostic score (within the current subgroup) ####
## Train gradient boosting machines to estimate prognostic scores
## Number of trees = 5000 (n.trees), maximum depth of each tree = 4 (interaction.depth)
set.seed(fitseed)
prognostic.gbm = gbm(reward~., data=dat_feature %>% dplyr::select(-ID, -treatment), 
                     distribution="gaussian", n.trees=5000, interaction.depth=4)
(prognostic.gbm$fit - dat_feature$reward)^2 %>% mean # mean squared error
## prognostic scores (standardized)
pat_prog_adj = scale(prognostic.gbm$fit)

#### Create dataset recording rewards R_i, treatments A_i, and features H_i (within the current subgroup) ####
pat_H_incluster = data.frame((pat_info_ITR %>%
                                filter(cluster == grp & treatment %in% trts) %>%
                                dplyr::select(ID, cluster, treatment, reward)),
                             dat_feature[,c(demo_name, bm_name)],
                             pat_prop_adj,
                             prog = pat_prog_adj
)
summary(pat_H_incluster)

dat_incluster = pat_H_incluster %>%
  mutate(reward_res = reward)

## Calculate the Euclidean distance between each subjects
dist_incluster = dat_incluster %>%
  dplyr::select(-ID, -cluster, -treatment, -reward, -reward_res) %>%
  dist() %>%
  as.matrix()
dim(dist_incluster)

idx_incluster = data.frame(ID=dat_incluster$ID, 
                           index=1:nrow(dat_incluster)) # index is used to select the row of the distance matrix

if(fitseed==2020){
  save(list=setdiff(ls(), c("rseed","rsd","casename","cs","jobname",
                            "maindir","codesdir","wsdir","fitseed")),
       file=file.path(wsdir,paste0(casename,"_workspace.RData")))
}


#### Learn ITR ####
ksvm.grid = expand.grid(C = 2^(-15:15)) # Potential values of the cost parameter in SVM
nfolds_obj = 1 # Number of folds in the cross-validations, 1 means using the whole dataset to learn ITR
nfolds = 2 # Number of folds in the cross-validations used to select the best tuning parameters in SVM
eps = 1e-16

set.seed(fitseed)
ITR_OVO = ksvm_cv_wvf_obj_knownfold_v2(
  dat=dat_incluster, idx=idx_incluster, trts=trts, SNN=SNN, nfolds_obj=nfolds_obj,
  kernel="rbfdot", kpar="automatic", eps=eps, nfolds=nfolds, tuneGrid=ksvm.grid,
  delta=max(dist_incluster), propensity=pat_prop_trunc_norm
)

## A data frame recording rewards, propensity scores, treatments, and ITRs
df_evf = data.frame(map_propsc(ITR_OVO$pred, pat_prop_trunc_norm), iter=fitseed) # column "vote" is the optimal individualized treatments

#### Save results ####
save(ITR_OVO, df_evf, file=file.path(wsdir,paste0(casename,"_est.RData")))

#### print results ####
## contingency table of treatment vs ITR
table(True=df_evf$treatment, ITR=df_evf$vote)

## empirical value function for one-size-fits-all rules
df_evf %>% 
  group_by(treatment) %>% 
  summarise(evf_ipw=sum(reward/pi)/sum(1/pi)) %>%
  as.data.frame() # column "evf_ipw"

## empirical value function for our ITR, if this value is greater than the above outputs, then our ITR should be better than above rules
df_evf_summary = df_evf %>%
  group_by(fold) %>%
  summarise(evf_ipw=sum((reward/pi)[treatment==vote])/sum((1/pi)[treatment==vote]),
            n=n()) %>%
  as.data.frame()
mean(df_evf_summary$evf_ipw)

