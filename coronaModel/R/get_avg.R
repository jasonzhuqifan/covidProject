get_avg <- function(param){
  #' Calculate the average for each columns
  #'
  #' @param param parameter matrix to calculate avarege
  #' @return calcualte average with 10 burn-in values
  #' @export 
  #' 
  # Discard the first 10 iterations of each chain as burin-in
  mean(param[10:length(param)])
}