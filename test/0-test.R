#' test 1
library(spatstat)
#source("R/logi.engine.R")
# fn <- "logi.engine"
# env <- as.environment("package:spatstat")
# assignInNamespace(fn, logi.engine_with_vb, envir=env)
#library(vblogistic)

#source("R/formula.R")

library(ppmvb)


x <- rStrauss(100, 0.1, 0.06)

set.seed(1)
f0 <- ppm(x, trend = ~x, interaction=Strauss(0.06), method="logi")
set.seed(1)
f1 <- ppmvb(x, trend = ~x, interaction=Strauss(0.06), verbose=T, xi=4, eps=0.1)


print(exp(rbind(coef(f0),coef(f1))))

