#' Fit SEIR model
#'
#' @param pred_days number of days we want to predict for.
#' @param a1 value of alpha1, the infectious rate by exposed people.
#' @param a2 value of alpha2, the infectious rate by infected people.
#' @param b value of beta, the rate of transition from exposed to infected.
#' @param s value of sigma, the rate of transition from infected to recovered.
#' @param g value of gmma, the rate of transition from infected to death.
#' @param E initial number of infected people with no symptom.
#' @param I initial number of confirmed infected people with symptom.
#' @param D initial number of death.
#' @param C initial number of recovered people.
#' @param a1_decay_rate the rate of decreasing alpha1, default is 1.
#' @param a2_decay_rate the rate of decreasing alpha2, default is 1.
#' @param s_decay_rate the rate of increasing sigma, default is 1.
#' @return The prediction i.e. confirmed, death, recovered based on the params for pred_days `(I, D, C)`
#' @export
seir_fit <- function(pred_days, a1, a2, b, s, g, E, I, D, C, a1_dec_rate=1, a2_dec_rate=1, s_inc_rate=1){
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
    s <- s_inc_rate * s
  }

  list(E=pred.E, I=pred.I, D=pred.D, C=pred.C)
}


#' Plot SIER models on confirmed cases, death cases and recovered cases.
#'
#' @param pred prediction values of confirmed cases, death cases and recovered cases. `(I, D, C)`
#' @param data true values of confirmed cases, death cases and recovered cases. `(I, D, C)`
#' @return The plots of confirmed cases, death cases and recovered cases.
#' @export
seir_plot <- function(pred, data){
  if (!is.null(pred$I) && !is.null(data$I)) {
    plot(pred$I, type="l", main="Confirmed", col="red", ylab = "Number of People", xlab = "Iterations", ylim = c(0, max(pred$I + data$I)))
    points(data$I, col="black")
  }
  if (!is.null(pred$D) && !is.null(data$D)) {
    plot(pred$D, type="l", main="Death", col="red", ylab = "Number of People", xlab = "Iterations", ylim = c(0, max(pred$D + data$D)))
    points(data$D, col="black")
  }
  if (!is.null(pred$C) && !is.null(data$C)) {
    plot(pred$C, type="l", main="Recovered", col="red", ylab = "Number of People", xlab = "Iterations", ylim = c(0, max(pred$C + data$C)))
    points(data$C, col="black")
  }
}

