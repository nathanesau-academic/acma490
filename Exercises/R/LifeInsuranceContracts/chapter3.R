## Exericse 3.2

setwd("D:/acma490/Exercises/R/LifeInsuranceContracts/")

pdf("exercise32_1.pdf")

SIGMA <- 0.02

oumodel = iratemodel(params = list(delta = 0.05, delta0 = 0.03, alpha = 0.10, 
                                   sigma = 0.02), "ou")
pv.var <- function(t) pv.moment(t,2,oumodel) - pv.moment(t,1,oumodel)^2
plot(pv.var, 0, 100, ylab = "Var[pv(t)]", xlab = "t")
maxt = optim(par=c(20), fn=function(t)-pv.var(t),
             method="Brent", lower=0, upper=50)$par

maxt
dev.off()

library(lattice)
library(latex2exp)

oumodel = iratemodel(params = list(delta = 0.05, delta0 = 0.03, alpha = 0.10, 
                                   sigma = 0.02), "ou")

g <- expand.grid(y = seq(1,70,5), x = seq(20,70,5))
g$z <- numeric(nrow(g))

counter = 1

for(i in 1:nrow(g))
{
  term = insurance(params = list(n = g$y[counter], d = 1), "isingle", "term")
  mort = mortassumptions(params = list(x = g$x[counter], table = "MaleMort82"))
  
  g$z[counter] = z.sd(term, mort, oumodel)
  counter = counter + 1
}

pdf("exercise32.pdf")

wireframe(z ~ y * x, data = g,
          drape = TRUE,
          col = 'black',
          col.regions = 'white',
          aspect = c(1.0, 0.8),
          colorkey = FALSE,
          xlab = "n",
          ylab = "issue age",
          zlab = "",
          screen = list(z = 340, x = -70),
          scales = list(arrows = FALSE, col="black", font = 10,
                        cex= 1.0),
          par.settings = list(regions=list(alpha = 0.3),
                              axis.line = list(col = "transparent")),
          zoom = 0.90)

dev.off()
