## Course: Acma 490
## Title:  Homework 1
## Author: Nathan Esau
## Date:   Jan 19, 2017

# Read in the weekly time series data
tsdata <- read.csv("D:/acma490/Homework/Homework1/hw1_data.csv")
tsdata$Date <- as.Date(tsdata$Date, format="%d-%b-%y")
tsdata <- tsdata[order(tsdata$Date),] # sort by date

# the time series to analyze
y <- tsdata$Series
x <- y - mean(y)

# 1 - fit AR(1) model by implementing least square method
fit.ar1 <- function(x) {
  
  N <- length(x)
  
  # based on (3.17) and (3.18)
  phi1 = sum(x[2:N] * x[1:(N-1)]) / sum(x[1:(N-1)]^2)
  sigma2 = sum((x[2:N] - phi1 * x[1:(N-1)])^2) / (N - 1)
  
  names(phi1) <- "ar1"
  
  list(coef = phi1, sigma2 = sigma2)
}

model1 <- fit.ar1(x)

# 2 - fit AR(1) model using a package
model2 <- arima(x, c(1,0,0), include.mean = FALSE, method = 'CSS')
model3 <- arima(x, c(1,0,0), include.mean = FALSE, method = 'CSS-ML')

# 3 - compare the fitted models
models <- data.frame(
  model = c("model1", "model2", "model3"),
  ar1 = c(model1$coef, model2$coef, model3$coef),
  sigma2 = c(model1$sigma2, model2$sigma2, model3$sigma2),
  call = c("fit.ar1", "arima(method = 'CSS')",
           "arima(method = 'CSS-ML')")
)

## Additional Code

# green's function for AR(1)
green.ar1 <- function(model, nterm = 10) {
  
  phi1 <- model1$coef
  lambda1 <- phi1
  g1 <- 1
  
  lambda1^(1:nterm)
}

green.ar1(model1)
ARMAtoMA(model1$coef, lag.max = 10)
