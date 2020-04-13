library("deSolve")
require(coronaModel)

# Import data
data.infected <- read.csv("time_series_covid19_confirmed_global.csv", sep=",", header = TRUE)
data.infected <- data.infected[data.infected[,"Country.Region"]=="Italy",]
infected <- colSums(data.infected[,5:ncol(data.infected)])

data.recovered <- read.csv("time_series_covid19_recovered_global.csv", sep=",", header = TRUE)
data.recovered<- data.recovered[data.recovered[,"Country.Region"]=="Italy",]
recovered <- colSums(data.recovered[,5:ncol(data.recovered)])

data.death <- read.csv("time_series_covid19_deaths_global.csv", sep=",", header = TRUE)
data.death <- data.death[data.death[,"Country.Region"]=="Italy",]
death <- colSums(data.death[,5:ncol(data.death)])

removed <- recovered + death

# Assume a fixed population of 60.36 billion in Italy
N <- 60360000
init_b = 0.28
init_g = 1/12
lower_g = 0.08
upper_g = 0.09
upper_b = 0.3
min_init_I = 20
max_init_R = 0

ii = 1
while (infected[ii] < min_init_I && removed[ii] <= max_init_R) {
  ii = ii + 1
}
infected = infected[ii:length(infected)]
removed = removed[ii:length(removed)]

init_I = infected[1]
init_R = removed[1]

# Output the optimal parameters beta and gamma respectively
opt_param <- sir_mle_param(init_b, init_g, lower_g, upper_g, upper_b, infected) # init beta, init gamma, lower_g, upper_g, upper_b

# Calculate the R0 number
R0.hat <- opt_param[1]/opt_param[2]

# Use ode to find the predictions for the SIR model
pred_ode <- sir_ode_fit(opt_param, c(1:length(infected)), N, N-init_I-init_R, init_I, init_R) # remove hardcode later TBD..

# Use iterations to find the predictions for the SIR model
pred <- sir_fit(length(infected), opt_param[1], opt_param[2], N - init_I, init_I, init_R, N)

data <- list(I=infected,R=removed)

sir_plot(pred, data)
