## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(brokenstick)

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("growthcharts/brokenstick@9b969af")

## ----fit-v1, eval = FALSE-----------------------------------------------------
#  data <- brokenstick::smocc_200
#  
#  # formula interface
#  fit1 <- brokenstick(hgt_z ~ age | id, data)
#  
#  # XY interface - numeric vector
#  # Deprecated in v2.0.0
#  fit2 <- with(data, brokenstick(age, hgt_z, id))
#  
#  # XY interface - data.frame
#  # Deprecated in v2.0.0
#  fit3 <- with(data, brokenstick(data.frame(age), hgt_z, id))
#  
#  # XY interface - matrix
#  # Deprecated in v2.0.0
#  tt <- as.matrix(data[, c(1, 2, 7)])
#  fit4 <- brokenstick(tt[, "age", drop = FALSE],
#                      tt[, "hgt_z", drop = FALSE],
#                      tt[, "id", drop = FALSE])

## ----fit-v2-------------------------------------------------------------------
data <- brokenstick::smocc_200

# formula interface
fit1 <- brokenstick(hgt_z ~ age | id, data)

## ----predict-v1---------------------------------------------------------------
# predict at observed data
p1 <- predict(fit1, data)

# predict at knots
p2 <- predict(fit1, data, x = "knots")

# predict at both observed data and knots
p3 <- predict(fit1, data, x = "knots", strip_data = FALSE)

# predict knots, broad matrix
p4 <- predict(fit1, data, x = "knots", shape = "wide")

## ----predict-v2---------------------------------------------------------------
# predict at observed data
p1 <- predict(fit1)

# predict at knots
p2 <- predict(fit1, x = "knots")

# predict at both observed data and knots
p3 <- predict(fit1, x = "knots", strip_data = FALSE)

# predict knots, broad matrix
p4 <- predict(fit1, x = "knots", shape = "wide")

## ----plot-v1, fig.height=3, fig.width=7---------------------------------------
ids <- c(10001, 10005, 10022)
plot(fit1, data, group = ids, what = "all")

## ----plot-v2, fig.height=3, fig.width=7---------------------------------------
ids <- c(10001, 10005, 10022)
plot(fit1, group = ids, what = "all")

