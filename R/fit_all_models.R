
# Load Packages -----------------------------------------------------------

#Set required packages
.cran_packages <- c("tidyverse",
                    "tidymodels",
                    "patchwork", 
                    "vegan", 
                    "seqinr",
                    "ape", 
                    "RColorBrewer",
                    "devtools",
                    "data.table",
                    "future",
                    "furrr")
.bioc_packages <- c("dada2",
                    "phyloseq", 
                    "ALDEx2",
                    "Biostrings")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

sapply(c(.cran_packages,.bioc_packages), require, character.only = TRUE)

library(speedyseq)
library(metacal)

source("helper_functions.R")

# Define functions --------------------------------------------------------


# Metacal spec ------------------------------------------------------------

metacal_spec <- function(x, denom="all", in_scale = NULL, out_scale="linear"){
  
  #input checks
  if(!is(x, "data.frame")){return(NULL)}    
  if(!denom %in% c("all", unique(x$taxon))){
    stop("denom must either be all, or a taxon name")
  }
  
  #Transform x to matrix
  mat <- x %>%
    dplyr::select(sample_id, taxon, error)%>%
    pivot_wider(id_cols= sample_id,
                names_from = "taxon",
                values_from = "error",
                values_fill = list(error=NaN))%>%    
    column_to_rownames("sample_id") %>%
    as.matrix() 
  
  # Check if matrix is log or linear
  if(is.null(in_scale)){
    curr_nas <- sum(is.na(mat))
    new_nas <- sum(is.na(log(mat)))
    if(curr_nas == new_nas){
      in_scale <- "linear"
      warning("No input scale provided, guessed ", in_scale)
    } else if(!curr_nas == new_nas){
      in_scale <- "log"
      warning("No input scale provided, guessed ", in_scale)
    }
  }
  
  # Drop missing data
  if (in_scale == "linear"){
    rows_to_keep <- rowSums(mat, na.rm=TRUE) > 0
    cols_to_keep <- colSums(mat, na.rm=TRUE) > 0
    
  } else  if(in_scale=="log"){
    rows_to_keep <- rowSums(is.na(mat)) < ncol(mat)
    cols_to_keep <- colSums(is.na(mat)) < nrow(mat)
  }
  mat <- mat[rows_to_keep, cols_to_keep]
  
  # Transform to log scale
  if (in_scale == "linear"){
    mat <- log(mat)
  }
  
  # log center proj
  K <- ncol(mat)
  mat0 <- mat
  mat0[is.nan(mat0)] <- 0
  
  P_sum <- diag(0, nrow = K)
  v_sum <- rep(0, K)
  weights = rep(1, nrow(mat))
  
  proj_mat <- function(K, M = c()) {
    mat <- diag(nrow = K) - 1/(K - length(M))
    mat[M,] <- 0
    mat[,M] <- 0
    mat
  }
  
  for (i in seq(nrow(mat))) {
    M <- which(is.nan(mat[i,]))
    v <- mat0[i,]
    P <- proj_mat(K, M) * weights[i]
    P_sum <- P_sum + P
    v_sum <- v_sum + P %*% v
  }
  
  #CLR transform
  clr_center <- (MASS::ginv(P_sum) %*% v_sum) %>% c
  names(clr_center) <- colnames(mat)
  
  out <- clr_center
  
  if (!denom=="all"){
    out <- out - mean(out[denom])
  }
  if (out_scale == "linear"){
    out <- exp(out)
  }
  out <- tibble::enframe(out, "taxon", ".pred")    
}


# Modelling func ----------------------------------------------------------
modelling_func <- function(train, test, cv, error_type, model_type, val_err = TRUE){
  if (error_type  %in% c("error_alr", "error_clr")){
    type <- "ratio"
  } else if(!error_type  %in% c("error_alr", "error_clr")){
    type <- "prop"
  }
  print(error_type)
  
  #Create model specification for each different type
  if(model_type == "lm"){
    model_spec <- 
      linear_reg() %>% 
      set_mode("regression") %>% 
      set_engine("lm")
  } else if (model_type == "lasso"){
    model_spec <- 
      linear_reg(penalty = tune(), mixture = tune()) %>% 
      set_mode("regression") %>% 
      set_engine("glmnet") %>%
      set_args(family="gaussian")
  } else if (model_type == "rf"){
    model_spec <- 
      rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
      set_mode("regression") %>% 
      set_engine("ranger", importance = "impurity")
  } else if (model_type == "xgboost"){
    model_spec <- 
      boost_tree(
        trees = 1000, 
        tree_depth = tune(), min_n = tune(), 
        loss_reduction = tune(),                     ## first three: model complexity
        sample_size = tune(), mtry = tune(),         ## randomness
        learn_rate = tune(),                         ## step size
      ) %>% 
      set_mode("regression") %>% 
      set_engine("xgboost") 
  } else if (model_type == "svm"){
    model_spec <- 
      svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>% 
      set_mode("regression") %>% 
      set_engine("kernlab") 
  } else if (model_type == "metacal"){
    if(error_type == "error_alr"){
      denom <- "Carpophilus_hemipterus"
    } else {
      denom <- "all"
    }
    model_obj <- metacal_spec(train, denom=denom, out_scale = "log")
    
    pred_train <- train %>%
      left_join(model_obj, by="taxon") %>%
      dplyr::select(.pred)
    
    pred_test <- test %>%
      left_join(model_obj, by="taxon") %>%
      dplyr::select(.pred)
    
    if(val_err){
      cv_error <- cv %>%
        mutate(.metrics = purrr::map(splits, function(x){
          cv_train <- x %>%
            analysis() %>%
            metacal_spec(denom=denom, out_scale="linear")
          x %>%
            assessment()%>%
            left_join(cv_train, by="taxon") %>%
            mutate(error_type = error_type,
                   estimated = .pred,
                   observed.prop = case_when(
                     error_type == "error_abs" ~ exp(observed.prop),
                     error_type  == "error_prop" ~ expit(observed.prop),
                     TRUE ~ observed.prop),
                   expected.prop = case_when(
                     error_type == "error_abs" ~ exp(expected.prop),
                     error_type  == "error_prop" ~ expit(expected.prop),
                     TRUE ~ expected.prop),
                   predicted = case_when(
                     error_type %in% c("error_prop", "error_abs") ~ estimated,
                     TRUE ~ expected.prop * estimated)
            ) %>%
            group_by(sample_id) %>%
            mutate_at(vars(predicted, observed.prop, expected.prop), ~ . / sum(.) ) %>%
            ungroup()  %>% 
            rmse(truth = expected.prop, estimate = predicted)
        }))%>% 
        unnest(.metrics) %>% 
        dplyr::select(fold= id, RMSE = .estimate)
    } else (cv_error <- NULL)
    
    out <- tibble(model_obj = list(model_obj),
                  pred_train = list(pred_train),
                  pred_test = list(pred_test),
                  cv_error = list(cv_error),
                  tune_res = NULL)
    return(out)
  } else if (model_type == "uncorrected"){
    model_obj <- train %>%
      dplyr::select(taxon) %>%
      distinct() %>%
      mutate(.pred = 0)
    pred_train <- train %>%
      left_join(model_obj, by="taxon") %>%
      dplyr::select(.pred)
    
    pred_test <- test %>%
      left_join(model_obj, by="taxon") %>%
      dplyr::select(.pred)
    
    if(val_err){
      cv_error <- cv %>%
        mutate(.metrics = purrr::map(splits, function(x){
          x %>%
            assessment()%>%
            rmse(truth = expected.prop, estimate = observed.prop)
        }))%>% 
        unnest(.metrics) %>% 
        dplyr::select(fold= id, RMSE = .estimate)
    } else (cv_error <- NULL)
    
    out <- tibble(model_obj = list(model_obj),
                  pred_train = list(pred_train),
                  pred_test = list(pred_test),
                  cv_error = list(cv_error),
                  tune_res = NULL
    )
    return(out)
  }
  
  # Create workflow
  if(type == "ratio"){
    # Tune ratio workflow
    model_workflow <- 
      workflow() %>% 
      add_recipe(bias_recipe) %>% 
      add_model(model_spec) 
  } else if(type == "prop"){
    # Tune proporiton workflow
    model_workflow <- 
      workflow() %>% 
      add_recipe(prop_recipe) %>%  #Proportions recipe
      add_model(model_spec) 
  } else (stop("type must be 'ratio' or 'prop' "))
  
  # Tune model if required
  if(!model_type %in% c("lm")){
    #Tune grid using a latin hypercube design
    model_tune <- tune_grid(model_workflow, resamples = cv, grid = 40)
    model_workflow <- finalize_workflow(model_workflow, model_tune %>% select_best("rmse"))
    tune_res <- collect_metrics(model_tune)
    
  }else {
    tune_res <- NULL
  }
  
  # Train on training set
  model_obj <- fit(model_workflow, data=train)
  
  # Predict training set
  pred_train <- safe_predict(model_obj, train)
  
  # Predict test set
  pred_test <- safe_predict(model_obj, test)
  
  # Get CV error
  if(val_err){
    cv_error <- fit_resamples(model_workflow,  resamples = cv) %>% 
      unnest(.metrics) %>% 
      filter(.metric == "rmse") %>%
      dplyr::select(fold= id, RMSE = .estimate)
  } else (cv_error <- NULL)
  
  out <- tibble(model_obj = list(model_obj),
                pred_train = list(pred_train),
                pred_test = list(pred_test),
                cv_error = list(cv_error),
                tune_res = list(tune_res))
  return(out)
}

# Load data ---------------------------------------------------------------
joint <- readRDS("joint.rds")

# Make training / testing splits  ------------------------------------------------------------

set.seed(666)

# Get 100 different random seeds
seeds <- sample(1:10000, 100, replace = FALSE)
s=1

split_list <- vector("list", length=length(seeds))
for (s in 1:length(seeds)){
  set.seed(seeds[s])
  
  joint_split <- joint  %>%
    ungroup() %>%
    dplyr::select(sample_name, material_type) %>%
    distinct() %>%
    initial_split(strata = material_type, prop = .8)
  
  # Training set
  joint_train <- joint %>%
    dplyr::select(-error) %>%
    filter(sample_name %in% training(joint_split)$sample_name) %>%
    ungroup()%>%
    pivot_longer(starts_with("error"),
                 names_to="error_type",
                 values_to="error") %>%
    mutate(observed.prop =case_when(
      error_type == "error_prop_logit" ~ logit(observed.prop),
      error_type == "error_prop_log" ~ log(observed.prop),
      error_type == "error_abs_log" ~ log(observed_abs),
      TRUE ~ observed.prop
    ),
    expected.prop = case_when(
      error_type == "error_prop_logit" ~ logit(expected.prop),
      error_type == "error_prop_log" ~ log(expected.prop),
      error_type == "error_abs_log" ~ log(expected),
      TRUE ~ expected.prop
    ))%>%
    filter(expected > 0) %>%
    group_by(error_type) %>%
    nest()%>%
    dplyr::rename(train=data)
  
  
  # Testing set
  joint_test <- joint %>%
    dplyr::select(-error) %>%
    filter(sample_name %in% testing(joint_split)$sample_name) %>%
    ungroup()%>%
    pivot_longer(starts_with("error"),
                 names_to="error_type",
                 values_to="error")%>%
    mutate(observed.prop =case_when(
      error_type == "error_prop_logit" ~ logit(observed.prop),
      error_type == "error_prop_log" ~ log(observed.prop),
      error_type == "error_abs_log" ~ log(observed_abs),
      TRUE ~ observed.prop
    ),
    expected.prop = case_when(
      error_type == "error_prop_logit" ~ logit(expected.prop),
      error_type == "error_prop_log" ~ log(expected.prop),
      error_type == "error_abs_log" ~ log(expected),
      TRUE ~ expected.prop
    ))%>%
    filter(expected > 0) %>%
    group_by(error_type) %>%
    nest() %>%
    dplyr::rename(test=data)
  
  # Create Cross validation folds for model tuning 
  train_cv <- joint_train %>%
    mutate(cv = purrr::map(train, function(x){
      x %>%
        group_vfold_cv(v = 5, group = sample_id) #Ensures samples are kept together
    })) %>%
    dplyr::select(-train)
  
  
  split_list[[s]] <- joint_train %>%
    left_join(joint_test, by="error_type") %>%
    left_join(train_cv, by="error_type")%>%
    mutate(seed = seeds[s])
  
}

all_equal(split_list[[1]]$test[[1]], split_list[[2]]$test[[1]])

#
models_to_fit <- c("lm", "lasso", "rf", "xgboost", "svm", "metacal", "uncorrected")

splits <- split_list %>%
  purrr::map(function(x){
    x %>% dplyr::slice(rep(1:n(), each = length(models_to_fit))) %>%
      mutate(model_type = models_to_fit)
  }) %>%
  bind_rows() %>%
  #filter(error_type %in% c("error_prop", "error_clr", "error_alr", "error_abs")) %>%
  ungroup()


# Preprocessing  ------------------------------------------------------------

# Recipe for fitting model including fcid + material_type
bias_recipe <- 
  recipe(formula = error ~ 0 + taxon + fcid + material_type, data = splits$train[[1]]) %>% 
  step_string2factor(fcid, material_type, taxon) %>% 
  step_novel(all_nominal(), -all_outcomes()) %>% 
  step_dummy(all_nominal(), -all_outcomes(), one_hot=TRUE) %>% 
  step_zv(all_predictors())
#step_normalize(-all_nominal()) 

#Recipe for fitting model to proportions
prop_recipe <- recipe(formula =  expected.prop ~ 0 + observed.prop + taxon + fcid + material_type, data = splits$train[[1]]) %>% 
  #step_string2factor(one_of(fcid, material_type)) %>% 
  step_string2factor(taxon) %>% 
  step_novel(all_nominal(), -all_outcomes()) %>% 
  step_dummy(all_nominal(), -all_outcomes(), one_hot=TRUE) %>% 
  step_zv(all_predictors())  %>%
  step_naomit(all_predictors())


# Fit modelling function  ------------------------------------------------------------
cores <- 48

future::plan(future::multiprocess, workers = cores)
fits_all <- splits %>%
  #filter(error_type %in% c("error_clr", "error_alr", "error_abs", "error_prop")) %>%
  mutate(
    fits = furrr::future_pmap(list(train, test, cv, error_type, model_type), modelling_func)
  ) %>% 
  unnest(fits)

# write out fits_all
saveRDS(fits_all, "fits_all.rds")

nrow(fits_all)
message("Job Complete")
