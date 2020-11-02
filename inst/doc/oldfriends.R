## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(brokenstick)

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("growthcharts/brokenstick@V0.62.1")

## ----fit-old, eval=FALSE------------------------------------------------------
#  data <- brokenstick::smocc_200
#  fit <- brokenstick(y = data$hgt.z, x = data$age, subjid = data$subjid)

## ----fit-new------------------------------------------------------------------
data <- brokenstick::smocc_200

# formula interface
fit1 <- brokenstick(hgt.z ~ age | id, data)

# XY interface - numeric vector
fit2 <- with(data, brokenstick(age, hgt.z, id))

# XY interface - data.frame
fit3 <- with(data, brokenstick(data.frame(age), hgt.z, id))

# XY interface - matrix
tt <- as.matrix(data[, c(1, 2, 7)])
fit4 <- brokenstick(tt[, "age", drop = FALSE],
                    tt[, "hgt.z", drop = FALSE],
                    tt[, "id", drop = FALSE])

## ----predict-old, eval = FALSE------------------------------------------------
#  # predict at observed data
#  p1 <- predict(fit)
#  
#  # predict at knots
#  p2 <- predict(fit, at = "knots")
#  
#  # predict at both observed data and knots
#  p3 <- predict(fit, at = "both")
#  
#  # predict knots, broad version
#  p4 <- predict(fit, at = "knots", output = "broad")

## ----predict-new--------------------------------------------------------------
# predict at observed data
p1 <- predict(fit1, data)

# predict at knots
p2 <- predict(fit1, data, x = "knots")

# predict at both observed data and knots
p3 <- predict(fit1, data, x = "knots", strip_data = FALSE)

# predict knots, broad matrix
p4 <- predict(fit1, data, x = "knots", shape = "wide")

## ----plot-old, eval = FALSE---------------------------------------------------
#  ids <- c(10001, 10005, 10022)
#  plot(fit, ids = ids)

## ----plot-new, fig.height=3, fig.width=7--------------------------------------
ids <- c(10001, 10005, 10022)
plot(fit1, data, group = ids, what = "all")

## ----explain-old, eval=FALSE--------------------------------------------------
#  get_pev(fit)

## ----explain-new--------------------------------------------------------------
get_r2(fit1, data)

