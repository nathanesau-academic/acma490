## Exercise 3.3

ar1.ev <- function(t, delta = 0.05, delta0 = 0.08, phi1 = 0.90)
{
  phi1^t * (delta0 - delta) + delta
}

ar1.cov <- function(s, t, phi1 = 0.90, sigma = 0.01)
{
  phi1^(abs(s-t)) * (1 - phi1^(2*min(s,t))) / (1 - phi1^2) * sigma^2
}

y.ev <- function(t)
{
  sum(ar1.ev(0:(t-1)))
}

y.var <- function(t)
{
  total <- 0
  for(s in 0:(t-1))
  {
    for(r in 0:(t-1))
    {
      total <- total + ar1.cov(s,r)
    }
  }
  total
}

y.cov <- function(s, t)
{
  total <- 0
  for(m in 0:(s-1))
  {
    for(r in 0:(t-1))
    {
      total <- total + ar1.cov(m,r)
    }
  }
  total
}

## Present Value function

pv.moment <- function(t, moment, mu = -y.ev(t), sigma2 = y.var(t))
{
  exp(mu*moment + 0.5*sigma2*moment^2) 
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

ann.ev(5) # 4.000048
sqrt(ann.var(5)) # 0.081630
