library("deSolve")

#' Maximize the log-likelihood function and compute initial estimates of beta and gamma
#'
#' @param init_b number of days we want to predict for.
#' @param init_g value of alpha1, the infectious rate by exposed people.
#' @param lower_g value of alpha1, the infectious rate by exposed people.
#' @param upper_g value of alpha1, the infectious rate by exposed people.
#' @param upper_b value of alpha1, the infectious rate by exposed people.
#' @return The prediction i.e. confirmed, death, recovered based on the params for pred_days `(I, D, C)`
#' @export
sir_mle_param <- function(init_b, init_g, lower_g, upper_g, upper_b){
  ll.pois <- function(theta, yi) {
    logl = 0
    for (i in length(confirmed$Date)) {
      yi <- rpois(1,lambda = confirmed$Total[i])
      logl <- logl + (yi * log(confirmed$Total[i]) - confirmed$Total[i])
    }
    logl
  }
  # Define the log-likelihood function assuming the observations follow poisson process
  mle <- optim(c(init_b, init_g), fn = ll.pois, control = list(fnscale = -1))
  beta2 <- mle$par[1]
  gamma2 <-mle$par[2]
  opt_param <- c(beta2,gamma2)
  while (gamma2 < lower_g | upper_g < gamma2 | beta2 <= gamma2 | upper_b < beta2) {
    mle <- optim(c(init_b, init_g), fn = ll.pois, control = list(fnscale = -1))
    beta2 <- mle$par[1]
    gamma2 <-mle$par[2]
    opt_param <- c(beta2,gamma2)
  }
  return(opt_param)
}

#' Define a helper function to enforce additional constraints on beta and gamma so that
#' the optimal parameters will be computed until they meet the requirements
#'
#' @param opt_param number of days we want to predict for.
#' @param N value of alpha1, the infectious rate by exposed people.
#' @param init_S value of alpha1, the infectious rate by exposed people.
#' @param init_I value of alpha1, the infectious rate by exposed people.
#' @param init_R value of alpha1, the infectious rate by exposed people.
#' @return The prediction i.e. confirmed, death, recovered based on the params for pred_days `(I, D, C)`
#' @export
sir_ode_fit <- function(opt_param, times, N, init_S, init_I, init_R){
  # define the ODE for SIR model
  sir <- function(t, y, parms) {
    beta <- parms[1]
    gamma <- parms[2]
    S <- y[1]
    I <- y[2]
    R <- y[3]
    return(list(c(S = -beta * S * I / N, I = beta * S * I / N - gamma * I, R = gamma * I)))
  }
  # Plug in the MLE to solve the ODE for SIR model
  res <- ode(y = c(init_S, init_I, init_R), times=times, func = sir, parms = opt_param)
  res
}

#' Define a helper function to enforce additional constraints on beta and gamma so that
#' the optimal parameters will be computed until they meet the requirements
#'
#' @param pred_days number of days we want to predict for.
#' @param beta value of alpha1, the infectious rate by exposed people.
#' @param gamma value of alpha1, the infectious rate by exposed people.
#' @param S value of alpha1, the infectious rate by exposed people.
#' @param I value of alpha1, the infectious rate by exposed people.
#' @param R value of alpha1, the infectious rate by exposed people.
#' @param N value of alpha1, the infectious rate by exposed people.
#' @return The prediction i.e. confirmed, death, recovered based on the params for pred_days `(I, D, C)`
#' @export
sir_fit <- function(pred_days, beta, gamma, S, I, R, N){
  pred.S <- c(S)
  pred.I <- c(I)
  pred.R <- c(R)
  
  for (i in 2:pred_days) {     
    S_ <- S - beta * S * I / N
    R_ <- R + gamma * I
    I_ <- N - S_ - R_ 
    S <- S_
    I <- I_
    R <- R_
    pred.S[i] <- S
    pred.I[i] <- I
    pred.R[i] <- R
  }
  list(S=pred.S, I=pred.I, R=pred.R)
}

#' Calculate quantile regerssion loss
#'
#' @param pred contact rate
#' @param data removed rate
#' @export
sir_plot <- function(pred, data){
  if (!is.null(pred$I) && !is.null(data$I)) {
    plot(pred$I, type="l", main="Infected", col="red", ylab = "Number of People", xlab = "Days", ylim = c(0, max(pred$I,data$I)))
    points(data$I, col="black")
    
    plot(pred_ode[,3],type="l", main="Infected - Mu", col="red", ylab = "Number of People", xlab = "Days", ylim = c(0, max(pred$I,data$I)))
    points(data$I, col="black")}
  
  if (!is.null(pred$R) && !is.null(data$R)) {
    plot(pred$R, type="l", main="Removed", col="red", ylab = "Number of People", xlab = "Days", ylim = c(0, max(pred$R,data$R)))
    points(data$R, col="black")
    
    plot(pred_ode[,3],type="l", main="Removed - Mu", col="red", ylab = "Number of People", xlab = "Days", ylim = c(0, max(pred$R,data$R)))
    points(data$R, col="black")}
}

