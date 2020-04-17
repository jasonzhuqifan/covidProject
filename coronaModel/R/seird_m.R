#' Fit SEIRD modified model
#'
#' @param pred_days Number of days we want to predict for.
#' @param a1 Value of alpha1, the infectious rate by exposed people.
#' @param a2 Value of alpha2, the infectious rate by infected people.
#' @param b Value of beta, the rate of transition from exposed to infected.
#' @param s Value of sigma, the rate of transition from infected to recovered.
#' @param g Value of gmma, the rate of transition from infected to death.
#' @param E Initial number of infected people with no symptom.
#' @param I Initial number of confirmed infected people with symptom.
#' @param D Initial number of death.
#' @param C Initial number of recovered people.
#' @param a1_dec_rate The rate of decreasing alpha1, default is 1.
#' @param a2_dec_rate The rate of decreasing alpha2, default is 1.
#' @details The reparametrized SEIRD modified model is
#' ```
#' E_new = E_old + alpha1 * E_old + alpha2 * I_old - beta * E_old
#' I_new = I_old + beta * E_old - sigma * I_old - gamma * I_old
#' D_old = D_old + gamma * I_old
#' C_old = C_old + sigma * I_old
#' ```
#' where `alpha1` and `alpha2` can be given a constant decaying rate by a1_dec_rate and a2_dec_rate respectively.
#' @return The prediction i.e. confirmed, death, recovered based on the params for pred_days `(I, D, C)`.
#' @export
seird_m_fit <- function(pred_days, a1, a2, b, s, g, E, I, D, C, a1_dec_rate=1, a2_dec_rate=1){
  pred.E <- c()
  pred.I <- c()
  pred.D <- c()
  pred.C <- c()

  for (i in 1:pred_days) {

    E_ <- E + a1 * E + a2 * I - b * E
    I_ <- I + b * E - s * I - g * I
    D_ <- D + g * I
    C_ <- C + s * I
    E <- E_
    I <- I_
    D <- D_
    C <- C_
    pred.E[i] <- E
    pred.I[i] <- I
    pred.D[i] <- D
    pred.C[i] <- C
    a1 <- a1_dec_rate * a1
    a2 <- a2_dec_rate * a2
  }

  list(E=pred.E, I=pred.I, D=pred.D, C=pred.C)
}


#' Plot SIERD modified models on confirmed cases, death cases and recovered cases.
#'
#' @param pred Prediction values of confirmed cases, death cases and recovered cases. `(I, D, C)`.
#' @param data True values of confirmed cases, death cases and recovered cases. `(I, D, C)`.
#' @details Plot the given data seperately with prediction and real data.
#' @return The plots of confirmed cases, death cases and recovered cases.
#' @export
seird_m_plot <- function(pred, data){
  if (!is.null(pred$I) && !is.null(data$I)) {
    plot(pred$I, type="l", main="SEIRD Confirmed", col="red", ylab = "Number of People", xlab = "Iterations", ylim = c(0, max(pred$I + data$I)))
    points(data$I, col="black")
    legend("topleft",c("Real data","Prediction"), cex=.9, col=c("black","red"), pch=c('o','--'))
  }
  if (!is.null(pred$D) && !is.null(data$D)) {
    plot(pred$D, type="l", main="SEIRD Death", col="red", ylab = "Number of People", xlab = "Iterations", ylim = c(0, max(pred$D + data$D)))
    points(data$D, col="black")
    legend("topleft",c("Real data","Prediction"), cex=.9, col=c("black","red"), pch=c('o','--'))
  }
  if (!is.null(pred$C) && !is.null(data$C)) {
    plot(pred$C, type="l", main="SEIRD Recovered", col="red", ylab = "Number of People", xlab = "Iterations", ylim = c(0, max(pred$C + data$C)))
    points(data$C, col="black")
    legend("topleft",c("Real data","Prediction"), cex=.9, col=c("black","red"), pch=c('o','--'))
  }
}

