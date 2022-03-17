#!/usr/bin/env R 


# https://cran.r-project.org/web/packages/randomForestSRC/index.html
library(randomForestSRC)
library(survival)

data(pbc, package = "randomForestSRC")


pbc.obj <- rfsrc(Surv(days,status) ~ ., pbc, importance = TRUE)
#find.interaction(pbc.obj, method = "vimp", nvar = 8)

plot(pbc.obj)


# Last res to death ----

# Van alle resecties tijd tot dood of laatste event nodig






