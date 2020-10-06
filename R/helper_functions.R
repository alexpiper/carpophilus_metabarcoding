#Create safe predict function
safe_predict <- function(model,data){
  if(is.null(data) | is.null(model) | is(model, "data.frame") ){
    return(NULL)
  }
  out <- tryCatch(
    {
      predict(model, data)
    },
    error=function(e) {
      # Choose a return value in case of error
      return(NULL)
    }
  )    
  return(out)
}

# Safe fit
safe_fit <- function(model, data) {
  out <- tryCatch(
    {
     fit(model, data)
    },
    error=function(e) {
      # Choose a return value in case of error
      return(NULL)
    }
  )    
  return(out)
}


#Logit transform
logit <- function(x){
  log(x) - log(1 - x)
}

# Inverse logit transform
expit <- function(x){
  exp(x) / (1+exp(x))
}