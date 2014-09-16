#' test 1
library(devtools)
load_all(".")

library(spatstat)
library(splines)
f0 <- ppm(nztrees, ~ bs(x,df=3), method="logi")
f1 <- ppmvb(nztrees, ~ bs(x,df=3))


print(exp(rbind(coef(f0),coef(f1))))

