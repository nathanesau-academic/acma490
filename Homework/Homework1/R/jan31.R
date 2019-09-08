## Course:      Acma 490
## Description: Case Study 1
## Author:      Nathan Esau

library(forecast)

tsdata <- read.csv("D:/acma490/hw1/hw1_shortdata.csv")

y <- tsdata$Series

plot(y, type = 'l')

fit.ar1 <- function(x) {
  
  N <- length(x)
  
  # based on (3.17) and (3.18)
  phi1 = sum(x[2:N] * x[1:(N-1)]) / sum(x[1:(N-1)]^2)
  sigma2 = sum((x[2:N] - phi1 * x[1:(N-1)])^2) / (N - 1)
  
  names(phi1) <- "ar1"
  
  list(coef = phi1, sigma2 = sigma2)
}

ar1model <- arima(y - mean(y), order=c(1,0,0), method='CSS', include.mean = FALSE)

x <- y - mean(y)
fit.ar1(x)

plot(forecast(ar1model, h=25))
