## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      fig.width = 7, fig.height = 3,
                      dev = "png")

## ----message = FALSE----------------------------------------------------------
require(brokenstick)
require(dplyr)
library(ggplot2)

## ----smocc_200----------------------------------------------------------------
data <- brokenstick::smocc_200
head(data, 3)

## ----plotcm, echo = FALSE, fig.height=7, warning=FALSE------------------------
ggplot(data[1:500, ], aes(x = age, y = hgt, group = id, color = as.factor(id))) +
  geom_line(size = 0.1) + 
  geom_point(size = 0.7) +
  scale_colour_viridis_d(option = "cividis") +
  xlab("Age (years)") +
  ylab("Length (cm)") +
  theme_light() +
  theme(legend.position = "none")

## ----plotsds, fig.height=5, warning=FALSE-------------------------------------
ggplot(data[1:500, ], aes(x = age, y = hgt_z, group = id, color = as.factor(id))) +
  geom_line(size = 0.1) + 
  geom_point(size = 0.7) +
  scale_colour_viridis_d(option = "cividis") +
  xlab("Age (years)") +
  ylab("Length SDS") +
  theme_light() +
  theme(legend.position = "none")

## ----figure1, warning = FALSE-------------------------------------------------
set.seed(123)
fit <- brokenstick(hgt ~ age | id, data)
ids <- c(10001, 10005, 10022)
plot(fit, group = ids, what = "all",
     xlab = "Age (years)", ylab = "Length (cm)")

## ----zscore, warning = FALSE--------------------------------------------------
fit0 <- brokenstick(hgt_z ~ age | id, data)
plot(fit0, group = ids, what = "all", 
     xlab = "Age (years)", ylab = "Length (SDS)")

## ----plotfit2, cache = TRUE, warning=FALSE------------------------------------
fit2 <- brokenstick(hgt_z ~ age | id, data = data, knots = 0:3)
plot(fit2, group = ids, xlab = "Age (years)", ylab = "Length (SDS)")

## -----------------------------------------------------------------------------
fit2

## ----fit9, cache = TRUE, warning = FALSE--------------------------------------
knots <- round(c(0, 1, 2, 3, 6, 9, 12, 15, 18, 24)/12, 4)
fit9 <- brokenstick(hgt_z ~ age | id, data = data, 
                    knots = knots, boundary = c(0, 3))

## ----plotfit9, echo = FALSE, warning = FALSE----------------------------------
plot(fit9, group = ids, xlab = "Age (years)", ylab = "Length (SDS)")

## -----------------------------------------------------------------------------
p1 <- predict(fit2)
head(p1)
identical(nrow(data), nrow(p1))

## -----------------------------------------------------------------------------
p2 <- predict(fit2, x = "knots")
head(p2)
nrow(p2)

## -----------------------------------------------------------------------------
p3 <- predict(fit2, x = "knots", strip_data = FALSE)
table(p3$.source)

## -----------------------------------------------------------------------------
get_r2(fit2)

## -----------------------------------------------------------------------------
get_r2(fit9)

## ----subjleveldata------------------------------------------------------------
subj <- data %>%
  select(id, sex, ga, bw) %>% 
  group_by(id) %>% 
  slice(1)
head(subj, 3)

## ----predictatknots-----------------------------------------------------------
bs <- predict(fit9, x = "knots", shape = "wide")
data <- bind_cols(subj, select(bs, -id))
head(data, 3)

## ----lm1----------------------------------------------------------------------
fit1_lm <- lm(`2` ~ sex + ga + I(bw / 1000), data = data)
summary(fit1_lm)

## ----lm2----------------------------------------------------------------------
fit2_lm <- lm(`2` ~ sex + ga + I(bw / 1000) + `0`, data = data)
summary(fit2_lm)

