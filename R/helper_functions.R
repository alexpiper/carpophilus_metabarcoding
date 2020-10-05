#Create safe predict function
safe_predict <- function(model,data){
  if(is.null(data) | is(model, "data.frame") ){
    return(NULL)
  }
  predict(model, data)
}


#Logit transform
logit <- function(x){
  log(x) - log(1 - x)
}

# Inverse logit transform
expit <- function(x){
  exp(x) / (1+exp(x))
}