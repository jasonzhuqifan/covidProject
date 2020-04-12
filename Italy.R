library("deSolve")
require(coronaModel)
# import Italy data
confirmed <- read.csv("Italy.csv", sep=",", header = TRUE)
confirmed.test <- confirmed[1:nrow(confirmed),]
confirmed.train <- confirmed[(nrow(confirmed)-20):nrow(confirmed)-28,]
removed <- read.csv("Italy-removed.csv", sep=",", header = TRUE)
removed.test <- removed[1:nrow(removed),]
removed.train <- removed[(nrow(removed)-20):nrow(removed)-28,]

# Assume a fixed population of 60.36 billion in Italy
N <- 60360000

# Output the optimal parameters beta and gamma respectively
opt_param <- sir_mle_param(0.28, 1/12, 0.08, 0.09, 0.3) # init beta, init gamma, lower_g, upper_g, upper_b

# Calculate the R0 number
R0.hat <- opt_param[1]/opt_param[2]

# Use ode to find the predictions for the SIR model
pred_ode <- sir_ode_fit(opt_param, confirmed$Date, N, N-20, 20, 0) # remove hardcode later TBD..

# Use iterations to find the predictions for the SIR model
pred <- sir_fit(49, opt_param[1], opt_param[2], 60359980, 20, 0, 60360000)

data <- list(I=confirmed$Total,R=removed$Total)

sir_plot(pred, data)
