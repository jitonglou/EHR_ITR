#### Customized Function for Random Forest ####
customRF = list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters = data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid = function(x, y, len = NULL, search = "grid"){}
customRF$fit = function(x, y, wts, param, lev, last, weights, classProbs, ...){
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict = function(modelFit, newdata, preProc = NULL, submodels = NULL){
  predict(modelFit, newdata)
}
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
  predict(modelFit, newdata, type = "prob")
}
customRF$sort <- function(x){x[order(x[,1]),]}
customRF$levels <- function(x){x$classes}


#### Matched set ####
create_matchedset = function(dat_MS, idx_MS, SNN, delta){
  pat_matchedset = matrix(NA, nrow=nrow(dat_MS), ncol=SNN)
  ## for each patient, find his/her matched set in patients having the same cluster/sex as him/her
  for (i in 1:nrow(dat_MS)){
    pat_info = data.frame(dat_MS[i,], distance = 0)
    ms_temp = data.frame(dat_MS, distance = dist_incluster[idx_MS$index,i]) %>%
      filter(sex == pat_info$sex & treatment != pat_info$treatment & distance < delta)
    if(nrow(ms_temp)>0){
      pat_matchedset[i,] = (ms_temp %>% arrange(distance))[1:SNN,] %>%
        mutate(reward_res_diff = pat_info$reward_res - reward_res) %>%
        select(reward_res_diff) %>% t()
    }
  }
  return(c(t(pat_matchedset))) # return a vector of length SNN*nrow(dat_MS)
} # easier case, fix-size-matched set

#### Map prediction and propensity score to the original data ####
map_pred = function(dat, pred){
  # dat, pred are data frames
  out = dat %>% select(ID, reward, treatment)
  out$vote = apply(pred, 1, function(z){
    names(table(z))[which.max(table(z))]
  })
  return(out)
} # used in the SVM

map_propsc = function(dat_pred, propsc){
  # dat_pred, propsc are data frames
  dat_pred$pi = 1
  for(i in 1:nrow(dat_pred)){
    dat_pred$pi[i] = propsc[i,dat_pred$treatment[i]]
  }
  return(dat_pred)  
} # used to create df_evf


#### Empirical value function ####
ValueFunc_propsc = function(dat, pred, propsc){
  fit_vote = apply(pred, 1, function(z){
    names(table(z))[which.max(table(z))]
  })
  row_idx = which(dat$treatment == fit_vote) # TRUE/FALSE
  if (length(row_idx)>0){
    denom = rep(1, length(row_idx)) # denom has the same length as row_idx
    for (i in 1:length(row_idx)){
      denom[i] = propsc[i,fit_vote[i]]
    }
    return(sum(dat$reward[row_idx]/denom)/sum(1/denom)) # DO NOT USE MEAN()
  } else {return(NA)}
}


#### Binary classifier (matched learning, given hyperparameters)####
my_ksvm = function(dat, idx, which.test=NULL, trts, SNN,
                   kernel="rbfdot", kpar="automatic", C, eps, delta){
  K = length(trts)
  ITR_OVO = vector("list", K*(K-1)/2) # a list of K(K-1)/2 classifiers
  
  if (!is.null(which.test)){ # do cross-validation
    pred_test = matrix(nrow = length(which.test), ncol = K*(K-1)/2)
  } else { # use whole data
    pred_test = matrix(nrow = nrow(dat), ncol = K*(K-1)/2)
  }

  for (l in 1:(K-1)){
    for (j in ((l+1):K)){
      idxc = l*(K-1)-l*(l-1)/2-(K-j)
      names(ITR_OVO)[idxc] = paste(trts[l], "vs", trts[j], sep="_")
      
      if (!is.null(which.test)){ # do cross-validation
        dat_forsvm = dat[-which.test,] %>%
          filter(treatment == trts[l] | treatment == trts[j]) %>%
          mutate(trt_bi = case_when(
            treatment == trts[l] ~ 1,
            TRUE ~ -1))
      } else { # use whole dataset
        dat_forsvm = dat %>%
          filter(treatment == trts[l] | treatment == trts[j]) %>%
          mutate(trt_bi = case_when(
            treatment == trts[l] ~ 1,
            TRUE ~ -1))
      }

      idx_forsvm = idx %>%
        semi_join(dat_forsvm, by="ID")
        
      ## difference of reward residuals for patients and their matched sets
      reward_res_diff = create_matchedset(dat_MS = dat_forsvm, 
                                          idx_MS = idx_forsvm,
                                          SNN = SNN,
                                          delta = delta) # length of 2*nrow(dat_forsvm)
      
      idx_reward_res_diff = rep(1:nrow(dat_forsvm),each=SNN)
      
      ## exclude patients who couldn't find a matchup or the difference between reward residuals is 0
      flag = (is.na(reward_res_diff)) | (reward_res_diff==0)
      
      dat_forsvm = dat_forsvm[idx_reward_res_diff[!flag],]
      idx_forsvm = idx_forsvm[idx_reward_res_diff[!flag],]
      reward_res_diff = reward_res_diff[!flag]

      ## Set sigma parameter(inverse bandwidth)
      if (kpar=="practical"){
        prac_sigma = dat_forsvm %>% 
          select(-ID, -cluster, -treatment, -reward, -reward_res, -trt_bi) %>%
          as.matrix() %>%
          dist() %>%
          median()
        inv_width=list(sigma=1/(2*prac_sigma^2)*0.5)
      } else {inv_width=kpar}
      
      ## save the binary classifier (Matched Learning)
      ITR_OVO[[idxc]] = weighted.ksvm(
        y=(dat_forsvm$trt_bi*sign(reward_res_diff)),
        x=(dat_forsvm %>% 
             select(-ID, -cluster, -treatment, -reward, -reward_res, -trt_bi) %>%
             as.matrix()),
        weights=g_func(reward_res_diff)/SNN,
        kernel=kernel, kpar=inv_width,
        C=C, eps=eps,
        nfolds=1, foldid = NULL
      )
      
      ## save the classifications for samples in the validate fold
      if (!is.null(which.test)){ # do cross-validation
        pred_test[,idxc] = ifelse(predict(ITR_OVO[[idxc]],
                                          newx=(dat[which.test,] %>%
                                                  select(-ID, -cluster, -treatment, -reward, -reward_res) %>%
                                                  as.matrix()
                                          ))>=0, trts[l], trts[j])
      } else { # use whole dataset
        pred_test[,idxc] = ifelse(predict(ITR_OVO[[idxc]],
                                          newx=(dat %>%
                                                  select(-ID, -cluster, -treatment, -reward, -reward_res) %>%
                                                  as.matrix()
                                          ))>=0, trts[l], trts[j])
      }
      
    } # loop j ends
  } # loop l ends
    
  return(list(model=ITR_OVO, predicted=pred_test))
}


#### CV for hyperparameters in the binary classifier ####
my_ksvm_cv_wvf = function(dat, idx, trts, SNN,
                          kernel="rbfdot", kpar="automatic", eps,
                          nfolds=3, tuneGrid, delta, propensity){
  best_idx = 1L
  best_param = tuneGrid[1,]
  cv_mat = cv_est <- NULL
  
  if (nfolds > 1 & nrow(tuneGrid) > 1){
    cv_mat = matrix(0, nrow = nrow(tuneGrid), ncol = nfolds)
    foldid = sample(rep(seq(nfolds), length = nrow(dat)))
    
    for (k in 1:nfolds){
      which.test = which(foldid == k) # a numeric vector recording index
      print(paste("Fold",k,"has",length(which.test),"test samples."))
      
      for (i in 1:nrow(tuneGrid)){
        fit = try(
          my_ksvm(dat=dat, idx=idx, which.test=which.test, trts, SNN,
                  kernel, kpar, 
                  C=tuneGrid[i,1], eps, delta),
          silent = TRUE
        ) # a list of length 2, $model and $predicted; or "typr-error" class; not every C can yield converged results
        
        if (class(fit) == "list"){
          cv_mat[i,k] = ValueFunc_propsc(dat[which.test,], fit$predicted, 
                                         propsc = propensity[which.test,])
        } else {
          cv_mat[i,k] = NA
        }
      }
      print(paste("Fold",k,"is done."))
    }
    cv_est <- rowMeans(cv_mat, na.rm = TRUE)
    best_idx <- which.max(cv_est)
    best_param <- tuneGrid[best_idx,]
  }
  
  print("C is selected.")
  print(paste0("Best C is 2^", log2(best_param)))
  
  best_fit = try(
    my_ksvm(dat=dat, idx=idx, which.test=NULL, trts, SNN,
            kernel, kpar, 
            C=best_param, eps, delta),
    silent = TRUE
  ) # even each fold converged, whole dataset might not converge
  
  if (class(best_fit) == "try-error") {
    best_idy = which.max(cv_mat[best_idx, ])
    
    best_fit = my_ksvm(dat=dat, idx=idx, which.test=which(foldid == best_idy), trts, SNN,
                       kernel, kpar, 
                       C=best_param, eps, delta)
  }
  
  print("best_fit is done.")
  print(summary(best_fit))
  
  return(list(best_fit = best_fit, params = tuneGrid, best_param = best_param, best_idx = best_idx, 
              cv_mat = cv_mat, cv_est = cv_est)
  )
}


#### CV for ITR ####
ksvm_cv_wvf_obj_knownfold_v2 = function(dat, idx, trts, SNN, nfolds_obj=3,
                                     kernel="rbfdot", kpar="automatic", eps,
                                     nfolds=3, tuneGrid, delta, propensity, 
                                     foldid_outer=NULL, foldid_inner=NULL){
  K = length(trts)
  fit_list = NULL
  dat_pred = NULL

  if (nfolds_obj == 1) {
    fit = my_ksvm_cv_wvf(dat, idx, trts, SNN, kernel, kpar, eps, 
                         nfolds, tuneGrid, delta, propensity)
    
    dat_pred = data.frame(map_pred(dat, fit$best_fit$predicted), fold=1)
    fit_list = list(fit)
    foldid_obj = rep(1, nrow(dat))
  } else if (nfolds_obj > 1){
    fit_list = vector("list", nfolds_obj)

    foldid_obj <- sample(rep(seq(nfolds_obj), length = nrow(dat)))
    if (!is.null(foldid_outer$fold)){
      if (all(sort(unique(foldid_outer$fold)) == seq(nfolds_obj))){
        foldid_obj = foldid_outer$fold
      }
    }
    
    for (k in 1:nfolds_obj){
      obj.test <- which(foldid_obj == k) # a numeric vector recording index
      print(paste("Obj fold",k,"has",length(obj.test),"test samples."))
      
      foldid_inner_k = do.call("c",lapply(foldid_inner, function(z){z[[k]]}))
      print(table(foldid_inner_k))
      fit = my_ksvm_cv_wvf_knownfold(dat=dat[-obj.test,], idx=idx[-obj.test,], 
                           trts, SNN, kernel, kpar, eps, 
                           nfolds, tuneGrid, delta, propensity=propensity[-obj.test,],
                           foldid_inner_k=foldid_inner_k)
      
      print("fit is done.")
      print(summary(fit))
      
      pred.test = sapply(fit$best_fit$model, function(z){
        predict(z, newx = (dat[obj.test,] %>% 
                             select(-ID, -cluster, -treatment, -reward, -reward_res) %>% as.matrix))
      })
      
      print("pred is done.")
      
      for (l in 1:(K-1)){
        for (j in ((l+1):K)){
          idxc = l*(K-1)-l*(l-1)/2-(K-j)
          pred.test[,idxc] = ifelse(pred.test[,idxc]>=0, trts[l], trts[j])
        }
      }
      
      dat_pred = dat_pred %>%
        rbind(data.frame(map_pred(dat[obj.test,], pred.test), fold=k)) 
      
      fit_list[[k]] = fit

      print(paste("Obj fold",k,"is done."))
    }
  }
  
  return(list(fits = fit_list, foldid_obj = foldid_obj, pred = dat_pred %>% arrange(ID))
  )
}


#### CV for Qlearning ####
qlearn_cv_wvf_obj = function(dat, idx, trts, nfolds_obj=3, nfolds=3,
                                mMain=".", mCont=NULL,
                                method="rf", scaled=FALSE,
                                tuneLength=NULL, tuneGrid=NULL, tuneMethod="cv",
                                foldid_outer=NULL, foldid_inner=NULL){
  # K = length(trts)
  fit_list <- NULL
  dat_pred = NULL

  if (nfolds_obj == 1) {
    # fit = qlearn_rf_cv_wvf(dat, idx, trts, nfolds, tuneGrid)
    # 
    # dat_pred = data.frame(map_pred(dat, fit$best_fit$predicted), fold=1)
    # fit_list = list(fit)
  } else if (nfolds_obj > 1){
    fit_list = vector("list", nfolds_obj)

    foldid_obj <- sample(rep(seq(nfolds_obj), length = nrow(dat)))
    if (!is.null(foldid_outer$fold)){
      if (all(sort(unique(foldid_outer$fold)) == seq(nfolds_obj))){
        foldid_obj = foldid_outer$fold
      }
    }

    for (k in 1:nfolds_obj){
      obj.test <- which(foldid_obj == k) # a numeric vector recording index
      print(paste("Obj fold",k,"has",length(obj.test),"test samples."))

      # if (!is.null(mMain)){
      #   if ((mMain %in% c("0","1",".")) | all(mMain %in% colnames(dat))){
      #     # inputs are vector of variable names...
      #   }
      # }
      
      form = paste0("reward~",mMain,"+treatment*(",mCont,")") %>% as.formula()
      
      fitControl <- trainControl(
        method = tuneMethod , # k-fold cross validation
        number = nfolds, # number of folds
        savePredictions = 'final'
      )
      
      fit = train(form = form,
                        data = dat[-obj.test,c(3:17)],
                        method = method, scaled=scaled,
                        tuneLength = tuneLength, tuneGrid = tuneGrid,
                        trControl = fitControl)
      
      print("fit is done.")
      print(fit$coefnames)
      # print(str(fit)) # long
      
      newdata = dat[rep(obj.test, each=length(trts)), c(3:17)]
      newdata$treatment = rep(trts, length(obj.test))
      print(dim(newdata))
      
      pred.test = predict(fit, newdata=newdata) %>%
        matrix(ncol = length(trts), byrow = TRUE)
      colnames(pred.test) = trts
      
      # pred.test = rep(0, length(obj.test))
      # for (i in seq(length(obj.test))){
      #   pred.test[i] = trts[which.max(result[(4*i-3):(4*i)])]
      # }
      
      print(dim(pred.test))
      print("pred is done.")

      # dat_pred = dat_pred %>%
      #   rbind(data.frame(map_pred(dat[obj.test,], pred.test), fold=k))
      
      dat_pred = dat_pred %>%
        rbind(data.frame(dat[obj.test,] %>%
                           select(ID, reward, treatment) %>%
                           mutate(vote=trts[apply(pred.test,1,which.max)]), 
                         fold=k))
      

      fit_list[[k]] = fit

      print(paste("Obj fold",k,"is done."))
    }
  }

  return(list(fits = fit_list, foldid_obj = foldid_obj, pred = dat_pred %>% arrange(ID))
  )
}

#### Applying learned ITRs ####
predict.ITR = function(object, newdata, treatments){
  ## create two treatment vectors trt1, trt2, which record the first and second treatment of each binary classifier, respectively 
  K = length(trts)
  trt1 = trt2 = rep(NA, K)
  for (l in 1:(K-1)){
    for (j in ((l+1):K)){
      idxc = l*(K-1)-l*(l-1)/2-(K-j)
      trt1[idxc] = trts[l]
      trt2[idxc] = trts[j]
    }
  } # the order of elements in trt1 and trt2 should match the order of names(object)
  
  ## provide predicted optimal treatments for newdata
  object %>% sapply(function(classifier){
      return(predict(classifier, newx=newdata))
      # return a N_pred*(K(K-1)/2) matrix, where N_pred = nrow(newdata)
      # each entity is either 1 or -1, indicating which treatment is recommended in the binary comparison
    }) %>% 
    apply(1, function(row){
      trt_vec = ifelse(row==1, trt1, trt2) # a vector of length K(K-1)/2, converting 1/-1 to treatment names
      return(names(table(trt_vec))[which.max(table(trt_vec))])
      # count which treatment wins the most binary comparisons
      # return a vector of length N_pred, of which each element is the optimal treatment for each subject
    }) %>%
    return()
}
