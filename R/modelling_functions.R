
# Metacal Spec ------------------------------------------------------------

# This is a modified version of the metacal::center function using the "proj" method
# This has been modified to auto detect the input scale and work better with the tidymodels workflow
# Original code from https://github.com/mikemc/metacal/blob/master/R/compositional-mean.R
# Method from vandenBoogaart2006; also described in vandenBoogaart2013 and at
# https://core.ac.uk/download/pdf/132548286.pdf (Bren2008)

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


## Tests to make sure that the pre transformed errors from earlier are equivalent
## They should be because everything we are looking at is a ratio
tests <- joint %>%
  dplyr::rename(error_rat = error) %>%
  filter(sample_id %in% training(joint_split)$sample_id) %>%
  ungroup()%>%
  pivot_longer(starts_with("error"),
               names_to="error_type",
               values_to="error") %>%
  filter(expected > 0) %>%
  group_by(error_type) %>%
  nest()

# Test that ALR on original errors is equivalent to the pre-transformed ALR error
test1 <- metacal_spec((tests %>% filter(error_type == "error_rat"))$data[[1]],
                      denom="Carpophilus_hemipterus")
test2 <- metacal_spec((tests %>% filter(error_type == "error_alr"))$data[[1]],
                      denom="Carpophilus_hemipterus")
all(round(test1$.pred,4) == round(test2$.pred,4))


#Test CLR on on original erorrs is equivalent to the pre-transformed CLR error
test3 <- metacal_spec((tests %>% filter(error_type == "error_rat"))$data[[1]])
test4 <- metacal_spec((tests %>% filter(error_type == "error_clr"))$data[[1]])
all(round(test3$.pred,4) == round(test4$.pred,4))