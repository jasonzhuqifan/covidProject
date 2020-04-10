#' Calculate quantile regerssion loss
#'
#'@param pred_days number of days we want to predict for
#' @param a1 value of alpha1
#' @param a2 value of alpha2
#' @param b value of beta
#' @param s value of sigma
#' @param g value of gmma
#' @param E initial number of infected people with no symptom
#' @param I initial number of confirmed infected people with symptom
#' @param D initial number of death
#' @param C initial number of recovered people
#' @return The prediction i.e. confirmed, death, recovered based on the params for pred_days `(I, D, C)`
#' @export
sier_predict <- function(pred_days, a1, a2, b, s, g, E, I, D, C){
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
  }

  list(E=pred.E, I=pred.I, D=pred.D, C=pred.C)
}


#' Calculate quantile regerssion loss
#'
#' @param a1 value of y-X*beta
#' @param a2 value of probability
#' @return z x ( tau - 0)$ if z < 0 or $ z x (tau - 1)  otherwise
#' @export
sier_plot <- function(pred, data){
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

