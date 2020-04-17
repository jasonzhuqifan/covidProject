#' Maximize the log-likelihood function and compute initial estimates of beta and gamma
#'
#' @param init_b Inital guess value of beta.
#' @param init_g Inital guess value of gamma.
#' @param lower_g The lower bound of gamma when using mle to optimize gamma.
#' @param upper_g The upper bound of gamma when using mle to optimize gamma.
#' @param upper_b The upper bound of beta when using mle to optimize beta.
#' @param confirmed A list of number of confirmed cases of real data `(day_1, ..., day_N)`.
#' @details Generate Yi that follows a poission distribution with a rate of I(ti).
#' Then derive a log-likelihood estimator and compute the initial beta and gamma by MLE.
#' @return The estimated parameter beta and gamma within the limit of bound values `(beta, gamma)`.
#' @export
sir_mle_param <- function(init_b, init_g, lower_g, upper_g, upper_b, confirmed){
  ll.pois <- function(theta, yi) {
    logl = 0
    for (i in length(confirmed)) {
      yi <- rpois(1,lambda = confirmed[i])
      logl <- logl + (yi * log(confirmed[i]) - confirmed[i])
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

#' Fit the SIR model using ODE to find the prediction values for susceptible, infected and removed cases
#' for a period of time, given gamma, beta and initial values of populations.
#'
#' @param opt_param The estimated parameter beta and gamma within the limit of bound values `(beta, gamma)`.
#' @param times A vector from 1 to the number of day we want to fit `(1, 2, ..., N)`.
#' @param N The value of otal population.
#' @param init_S The value of initial susceptible people.
#' @param init_I The value of initial infected people.
#' @param init_R The value of initial removed people.
#' @details The reparametrized SIR model is
#' ```
#' S_new = S_old - beta * S_old * I_old / N
#' I_new = I_old + beta * S_old * I_old / N - gamma * I_old
#' R_new = R_old + gamma * I_old
#' ```
#' @return The prediction i.e. susceptible, infected, removed `(S, I, C)`.
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

#' Fit SIR model using an iteration approach
#'
#' @param pred_days Number of days we want to predict for.
#' @param beta Value of beta, the rate of transition from susceptible to infected.
#' @param gamma Value of gamma, the rate of transition from infected to removed.
#' @param S Initial number of susceptible people.
#' @param I Initial number of confirmed infected people.
#' @param R Initial number of removed cases.
#' @param N Total population.
#' @details The reparametrized SEIRD modified model is
#' ```
#' S_new = S_old - beta * S_old * I_old / N
#' I_new = I_old + beta * S_old * I_old / N - gamma * I_old
#' R_new = R_old + gamma * I_old
#' ```
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

#' Plot SIR models on confirmed cases and removed cases.
#'
#' @param pred Prediction values of confirmed cases and removed cases. `(I, R)`.
#' @param data True values of confirmed cases and removed cases. `(I, R)`.
#' @details Plot the models of infection and removed cases seperately with prediction and real data.
#' @return Two plots of infection cases and cremoved cases respectively.
#' @export
sir_plot <- function(pred, pred_ode, data){
  if (!is.null(pred$I) && !is.null(data$I)) {
    plot(pred$I, type="l", main="SIR Infected", col="red", ylab = "Number of People", xlab = "Days", ylim = c(0, max(pred$I,data$I)))
    points(data$I, col="black")
    legend("topleft",c("Real data","Prediction"), cex=.9, col=c("black","red"), pch=c('o','--'))
    }
  if (!is.null(pred$R) && !is.null(data$R)) {
    plot(pred$R, type="l", main="SIR Removed", col="red", ylab = "Number of People", xlab = "Days", ylim = c(0, max(pred$R,data$R)))
    points(data$R, col="black")
    legend("topleft",c("Real data","Prediction"), cex=.9, col=c("black","red"), pch=c('o','--'))
  }
}

