## Exercise 4.4, 4.5, 4.6

ygbm.ev <- function(t, delta = 0.07)
{
   delta * t
}

ygbm.cov <- function(s, t, sigma = 0.1)
{
  sigma^2 * (min(s,t)^2 * max(s, t) / 2 - min(s,t)^3 / 6)
}

ygbm.var <- function(t, sigma = 0.1)
{
  sigma^2 * t^3 / 3
}

y.ev <- function(t, delta = 0.05, delta0 = 0.08, alpha = 0.10)
{
  delta*t + (delta0 - delta) * (1 - exp(-alpha*t)) / alpha
}

y.var <- function(t, alpha = 0.10, sigma = 0.01)
{
  sigma^2/alpha^2 * t + sigma^2 / (2*alpha^3) * (-3 + 4*exp(-alpha*t) - exp(-2*alpha*t))
}

y.cov <- function(s, t, alpha = 0.10, sigma = 0.01)
{
  sigma^2 / alpha^2 * min(s, t) +
  sigma^2 / (2*alpha^3) * (-2 + 2*exp(-alpha*s) + 2*exp(-alpha*t) 
                           - exp(-alpha*(abs(t-s))) - exp(-alpha*(t+s)))
}

## Present Value function

pvgbm.moment <- function(t, moment, mu = -ygbm.ev(t), sigma2 = ygbm.var(t))
{
  exp(moment*mu + 0.5*sigma2*moment^2)
}

pvgbm.ev <- function(t, mu = -ygbm.ev(t), sigma2 = ygbm.var(t))
{
  pvgbm.moment(t, moment = 1, mu = mu, sigma2 = sigma2)
}

pvgbm.cov <- function(s, t)
{
  mu <- -(ygbm.ev(s) + ygbm.ev(t))
  sigma2 <- ygbm.var(s) + ygbm.var(t) + 2*ygbm.cov(s,t)
  EXY <- exp(mu + 0.5*sigma2)
  EXEY <- pvgbm.ev(s) * pvgbm.ev(t)
  
  EXY - EXEY
}

pvgbm.var <- function(t)
{
  pvgbm.moment(t, moment = 2) - pvgbm.moment(t, moment = 1)^2
}

pv.moment <- function(t, moment, mu = -y.ev(t), sigma2 = y.var(t))
{
  exp(moment*mu + 0.5*sigma2*moment^2) 
}

pv.ev <- function(t, mu = -y.ev(t), sigma2 = y.var(t))
{
  pv.moment(t, moment = 1, mu = mu, sigma2 = sigma2)
}

pv.cov <- function(s, t)
{
  mu <- -(y.ev(s) + y.ev(t))
  sigma2 <- y.var(s) + y.var(t) + 2*y.cov(s,t)
  EXY <- exp(mu + 0.5*sigma2)
  EXEY <- pv.ev(s) * pv.ev(t)
  
  EXY - EXEY
}

pv.var <- function(t)
{
  pv.moment(t, moment = 2) - pv.moment(t, moment = 1)^2
}

## Annuity

ann.ev <- function(n)
{
  total = 0
  
  for(t in 1:n)
  {
    total = total + pv.ev(t)
  }
  
  total
}

ann.var <- function(n)
{
  total = 0
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      total = total + pv.cov(i, j)
    }
  }
  total
}

pv.ev(t = 10) # 0.50599342
pv.var(t = 10) # 0.0043000
y.cov(s = 5, t = 10) # 0.005957969
pv.cov(s = 5, t = 10) # 0.00209574

ann.ev(n = 10) # 6.88418889
ann.var(n = 10) # 0.13503997

pvgbm.ev(t = 10) # 0.61672421
pvgbm.var(t = 10) # 0.01289196
ygbm.cov(s = 5, t = 10) # 0.010416667
pvgbm.cov(s = 5, t = 10) # 0.0050398647