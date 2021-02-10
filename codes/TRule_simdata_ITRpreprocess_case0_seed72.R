rm(list=ls())
#### Libraries ####
library(dplyr) # easier data wrangling 
library(purrr) # reduce()

#### Set directories ####
jobname =  "TRule_simdata_ITRpreprocess"
rseed = 72; rsd = paste0("seed",rseed) # random seed used to generate one data set
cs = paste0("case", 0)
casename = paste(jobname,cs,rsd,sep="_")

maindir = "/pine/scr/j/i/jitong/thesis/TRule_model" # project directory
codesdir = file.path(maindir,"codes") # code directory
wsdir = file.path(maindir,"rdata") # R workspace directory

#### Load files ####
load(file.path(wsdir,paste0("TRule_simdata_GenerateData_",rsd,"_workspace.RData")))
load(file.path(wsdir,paste0("TRule_simdata_DistMat_",cs,"_",rsd,"_dendrogram.RData")))

#### Model settings ####
Obs_allbm = time_Xikj_all_df_list # measurements and time points
m = ncol(X)-1 # number of baseline variables
p = length(Obs_allbm) # number of health markers

#### Identify latent subject subgroups ####
num = 4 # number of group
clus = cutree(hc, num) # group index of each subject
table(clus)

#### Create a data frame recording group indices to subjects ####
pat_cluster = time_Xikj_all_df_list %>% 
  lapply(function(df){
    df$i %>% unique
  }) %>% 
  do.call(union, .) %>%
  sort %>%
  data.frame(ID = .,
             cluster = clus) %>%
  filter(cluster != 4)
pat_cluster$ID %>% is.unsorted # check if subject IDs are ordered, should be FALSE

#### Add demographic variables ####
pat_demo = data.frame(ID = 1:nrow(X), X) %>%
  inner_join(pat_cluster, by="ID")

#### Calculate recent pattern of health markers ####
Obs_allbm_recent = Obs_allbm %>%
  lapply(function(df){
    df %>%
      semi_join(pat_demo, by=c("i"="ID")) %>%
      filter(time > 3 & time < 11) %>%
      dplyr::select(i, time, value) # raw_data
  })
names(Obs_allbm_recent) = paste0("V", 1:p)

pat_bmpatterns = Obs_allbm_recent %>%
  lapply(function(df){
    df %>%
      group_by(i) %>%
      summarise(mean=mean(value)) # calculate the average measurements, between t=3 and t=11, for each subject and all health markers
  }) %>%
  reduce(inner_join, by="i")
colnames(pat_bmpatterns) = c("ID",names(Obs_allbm_recent)) # set the variable names to V1 and V2
summary(pat_bmpatterns)

#### Assign treatments and rewards to subjects ####
pat_trt = pat_demo %>%
  inner_join(pat_bmpatterns, by="ID") %>%
  mutate(lc = 1 + age - sex + 2*V1 - 2*V2)

## In all groups, randomly assign one of the three treatments(A/B/C) to each subject
## For each subject, reward = 1 + age - sex + 2*V1 -2*V2 + random error, and random error follows a normal distribution with mean 0 and standard deviation 1/3.
set.seed(2020)
pat_info_ITR = pat_trt %>%
  mutate(eps = rnorm(nrow(pat_trt), 0, 1/3),
         treatment = sample(c("A","B","C"), nrow(pat_trt), replace = TRUE)) %>%
  mutate(reward = lc + eps)

table(pat_info_ITR$cluster, pat_info_ITR$treatment) # number of subjects by treatments and groups

pat_info_ITR %>%
  group_by(cluster, treatment) %>%
  summarise(reward_avg=mean(reward)) # mean reward by treatments and groups

save(pat_info_ITR, file=file.path(wsdir,paste0(casename,"_df.RData")))

# save.image(file=file.path(wsdir,paste0(casename,"_workspace.RData")))