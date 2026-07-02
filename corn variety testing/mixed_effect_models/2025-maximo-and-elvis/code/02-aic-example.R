library(nlraa)
library(nlme)

## This is the script we wrote on Feb 13th looking at
## why comparing the AIC of the same model fit to
## different data is not a good idea. The AIC will vary
## greatly with the units of the data, just like it does
## in this example.

x <- 0:300
y <- linp(x, a = 7, b = 0.03, xs = 180)
y <- y + rnorm(length(x), mean = 0, sd = 0.5)


plot(x, y)

## data in t/ha
fit1 <- nls(y ~ SSlinp(x, a, b, xs))
fit1 |> AIC()



## data in kg/ha
y2 <- y * 1000
plot(x, y2)

fit2 <- nls(y2 ~ SSlinp(x, a, b, xs))
fit2 |> AIC()


