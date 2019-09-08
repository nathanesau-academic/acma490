# Read in the weekly time series data
tsdata <- read.csv("D:/acma490/Notes/Chapter5/hw1_data.csv")
tsdata$Date <- as.Date(tsdata$Date, format="%d-%b-%y")
tsdata <- tsdata[order(tsdata$Date),] # sort by date

# the time series to analyze
y <- tsdata$Series
x <- y - mean(y)

ARtoOU <- function(armodel, Delta = 1) {
  phi1 <- armodel$coef
  sigma2a <- armodel$sigma2
  
  alpha <- -log(phi1) / Delta
  sigma2 <- 2*alpha*sigma2a / (1 - phi1^2)
  
  list(alpha = alpha, sigma2 = sigma2)
}

OUtoAR <- function(oumodel, Delta = 1) {
  alpha = oumodel$alpha
  sigma2 <- oumodel$sigma2
  
  phi1 <- exp(-alpha*Delta)
  sigma2a <- sigma2 / (2*alpha) * (1 - phi1^2)
  
  list(coef = phi1, sigma2 = sigma2a)
}

# EstimateSecondOrderSDE <- function(x) {
# 
#   N <- length(x)
#   
#   f <- function(alpha, Delta = 1) {
#   
#     alpha0 <- alpha[1]
#     alpha1 <- alpha[2]
#       
#     mu <- polyroot(c(-alpha0,-alpha1,1))
#     mu1 <- mu[1]
#     mu2 <- mu[2]
#     
#     phi1 <- exp(mu1*Delta) + exp(mu2*Delta)
#     phi2 <- -exp((mu1 + mu2)*Delta)
#     
#     a <- alpha1/2
#     b <- sqrt(as.complex(alpha1^2 + 4*alpha0))/2
# 
#     P <- ifelse(alpha1^2 >= 4*alpha0,
#                   (b*sinh(2*a*Delta) - sinh(2*b*Delta))/
#                   (2*a*sinh(b*Delta)*cosh(a*Delta) - 2*b*sinh(a*Delta)*cosh(b*Delta)),
#                 (b*sinh(2*a*Delta) - a*sin(2*b*Delta))/
#                 (2*a*sin(b*Delta)*cosh(a*Delta) - 2*b*sinh(a*Delta)*cos(b*Delta)))
#     
#     lambda1 <- exp(mu1*Delta)
#     lambda2 <- exp(mu2*Delta)
# 
#     P1 <- (-mu1 * (1 + lambda1^2) * (1 - lambda2^2) + mu2 * (1 + lambda2^2) * (1 - lambda1^2)) /
#       (2*(mu1 * lambda1 * (1 - lambda2^2) - mu2 * lambda2 * (1 - lambda1^2)))
#    
#     print(alpha) 
#     print((P - P1))
#     
#     theta1 <- c(-P + sqrt(P^2 - 1), -P - sqrt(P^2 - 1))
#     theta1 <- theta1[which(abs(theta1) == min(abs(theta1)))]
#     
#     residuals <- numeric(N)
#     
#     for(i in 3:N) {
#       pred <- phi1 * x[i-1] + phi2 * x[i-2] + theta1 * residuals[i-1]
#       residuals[i] <- x[i] - pred
#     }
#   
#     Re(sum(residuals^2))
#   }
# 
#   alpha <- optim(par = c(-0.04, -0.5), fn = f,
#                method = 'L-BFGS-B', upper = c(-0.000001, -0.000001))$par
# 
#   print(alpha)
# }

ARMAtoSDE <- function(armamodel, Delta = 1) {
  phi1 <- armamodel$coef[1]
  phi2 <- armamodel$coef[2]
  
  lambda1 <- 0.5 * (phi1 + sqrt(phi1^2 + 4*phi2))
  lambda2 <- 0.5 * (phi1 - sqrt(phi1^2 + 4*phi2))
  
  mu1 <- log(lambda1) / Delta
  mu2 <- log(lambda2) / Delta
  
  alpha1 = mu1 + mu2
  alpha0 = -mu1 * mu2
  
  list(alpha0 = alpha0, alpha1 = alpha1)
}

SDEtoARMA <- function(sdemodel, Delta = 1) {
  alpha0 <- sdemodel$alpha0
  alpha1 <- sdemodel$alpha1
  sigma2 <- sdemodel$sigma2
  
  # eigenvalues
  mu1 <- 0.5 * (alpha1 + sqrt(alpha1^2 + 4*alpha0))
  mu2 <- 0.5 * (alpha1 - sqrt(alpha1^2 + 4*alpha0))
  
  lambda1 <- exp(mu1 * Delta)
  lambda2 <- exp(mu2 * Delta)
  
  phi1 <- lambda1 + lambda2
  phi2 <- -lambda1 * lambda2
  
  P = (-mu1*(1+lambda1^2)*(1-lambda2^2) + mu2*(1+lambda2^2)*(1-lambda1^2)) /
       (mu1*lambda1*(1-lambda2^2) - mu2*lambda2*(1-lambda1^2))
  
  theta1 <- Re(polyroot(c(1,P,1)))
  theta1 <- theta1[which(abs(theta1) == min(abs(theta1)))]
  
  sigma2a = sigma2 * (lambda1 - lambda2)^2 / (2*mu1*(mu1^2 - mu2^2) * (lambda1 - theta1)) /
    ((lambda1 - theta1)/(1 - lambda1^2) - (lambda2 - theta1)/(1 - lambda1*lambda2))
  
  coef <- c(phi1, phi2, theta1)
  names(coef) <- c("ar1", "ar2", "ma1")
  
  list(coef = coef, sigma2 = sigma2a)
}

armodel <- arima(x, order = c(1, 0, 0),
                 method = 'CSS',
                 include.mean = FALSE)

oumodel <- ARtoOU(armodel)
armodel <- OUtoAR(oumodel)

sdemodel <- list(
  alpha0 = -0.04,
  alpha1 = -0.5,
  sigma2 = 0.01
)

# armamodel <- arima(x, order = c(2, 0, 1),
#                    method = 'CSS',
#                    include.mean = FALSE)

armamodel <- list(coef = c(1.05, -0.095, -0.05))

sdemodel <- ARMAtoSDE(armamodel, Delta = 5)
sdemodel$sigma2 = 0.01 # arbitrary here
SDEtoARMA(sdemodel, Delta = 5)
