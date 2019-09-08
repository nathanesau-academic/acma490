tsdata <- read.csv("D:/acma490/Notes/Chapter3//hw1_data.csv")
tsdata$Date <- as.Date(tsdata$Date, format="%d-%b-%y")
tsdata <- tsdata[order(tsdata$Date),] # sort by date

# the time series to analyze
x <- tsdata$Series

# Estimate parameters for AR model

fit.ar <- function(x, p = 1) { # similar to arima()
  
  N <- length(x)
  
  # don't estimate an intercept
  xx <- embed(x, p+1)
  phi <- lm(xx[,1] ~ xx[,-1] - 1)$coef
  sigma2 <- sum((xx[,1] - as.matrix(xx[,2:(p+1)]) %*% phi)^2) / (N - 1)
      
  names(phi) <- sapply(1:p, function(x) paste0("ar",x))
  
  list(coef = phi, sigma2 = sigma2)
}

model <- fit.ar(x, p = 1)

# Green's function for an AR model

green.ar <- function(model, nterm = 10) { # similar to ARMAtoMA()
  phi <- model$coef[grepl("ar", names(model$coef))]
  
  lambda <- polyroot(c(-rev(phi), 1))
  
  n <- length(phi)
  g <- numeric(n)

  if(n == 1)
  {
    g[n] = 1
  } else {
    for(j in 1:n) {
      
      g[j] = 1
  
      for(i in 1:n) {
        if(i != j) {
          g[j] = g[j] / (lambda[j] - lambda[i])
        }
      }
    } 
  }
  
  G <- numeric(nterm)
  offset = n-1
  for(j in 1:nterm + offset) {
    G[j-offset] = sum(g * lambda^j)
  }
  
  Re(G)
}

model <- fit.ar(x, p = 20)
green.ar(model, 10)
ARMAtoMA(model$coef, lag.max=10)
