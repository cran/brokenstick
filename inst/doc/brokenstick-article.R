## ---- include=FALSE-----------------------------------------------------------
options(tinytex.verbose = TRUE)

## ----setup, include=FALSE-----------------------------------------------------
require("brokenstick")
require("dplyr")
require("tidyr")
require("ggplot2")
require("lattice")

old <- options(digits = 2)
knitr::opts_chunk$set(echo = TRUE, cache = FALSE, fig.height = 7, fig.width = 7, out.width="100%", dev = "png")

colors <- mice::mdc(1:6)
knitr::knit_hooks$set(set.theme = function(before, options, envir) {
  if (before) 
    lattice::trellis.par.set(list(
      strip.background = list(col = "grey95"),
      superpose.symbol = list(col = colors[1:2]),
      superpose.line = list(col = colors[1:2]),
      plot.symbol = list(col = colors[1:2])))
})
# knitr::knit_hooks$set(set.palette = function(before, options, envir) {
#   if (before) {
#     opal <- palette(colors)
#     opar <- par(cex.main = 1.0)
#   }
# })

## ----visits, fig.height = 4, fig.cap = "Abacus plot of observation times for the first 20 children of the SMOCC data.", echo=FALSE, set.theme = TRUE----
data <- smocc_200
data$month <- data$age * 12
schedule <- c(0, 0.92, 1.84, 3, 6, 9, 12, 15, 18, 24)
schedule_labels <- c("b", "4w", "8w", "3m", "6m", "9m", 
                     "12m", "15m", "18m", "24m")
lattice::dotplot(
  reorder(id, rev(id)) ~ month, data = data[1:200, ],  
  pch = 19,  col = colors[4], xlab = "Age (in months)",
  cex = 0.7, 
  scales = list(x = list(tck = 0, at = schedule,
                         labels = schedule_labels)),
  panel = function(...) {
    panel.refline(v = schedule)
    panel.dotplot(...)
  }
)

## ----plotfit2echo, echo = TRUE, eval = FALSE----------------------------------
#  ids <- c(10001, 10005, 10022)
#  fit2 <- brokenstick(hgt.z ~ age | id, smocc_200, knots = 0:3)
#  knots <- c(0, 0.0833, 0.1667, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)
#  fit9 <- brokenstick(hgt.z ~ age | id, smocc_200,
#                      knots = knots, boundary = c(0, 3))
#  m2 <- plot(fit2, smocc_200, group = ids,
#             xlab = "Age (years)", ylab = "Length (SDS)")
#  m9 <- plot(fit9, smocc_200, group = ids,
#             xlab = "Age (years)", ylab = "Length (SDS)")
#  gridExtra::grid.arrange(m2, m9, nrow = 2)

## ----plotfit2x, echo = FALSE, fig.height = 6, fig.cap = "Broken stick model with two (top) and nine (bottom) line segments for three children. Blue = observed data, Red = Fitted broken stick curves.", set.theme = TRUE----
ids <- c(10001, 10005, 10022)
fit2 <- brokenstick(hgt.z ~ age | id, smocc_200, knots = 0:3)
m2 <- plot(fit2, smocc_200, group = ids,
           xlab = "Age (years)", ylab = "Length (SDS)")
m9 <- plot(brokenstick::fit_200, smocc_200, group = ids,
           xlab = "Age (years)", ylab = "Length (SDS)")
gridExtra::grid.arrange(m2, m9, nrow = 2)

## -----------------------------------------------------------------------------
library(splines)
data <- brokenstick::smocc_200
brk <- c(0, 0.5, 1, 2)
X <- bs(data$age, knots = brk, Boundary.knots = c(0, 3), degree = 1)
colnames(X) <- paste("age", c(brk, 3), sep = "_")
data <- cbind(data[, c("id", "age", "hgt.z")], X)
head(data)

## ----warning=FALSE, echo=FALSE------------------------------------------------
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## -----------------------------------------------------------------------------
library(lme4)
ctl <- .makeCC("warning", tol = 4e-3)
f <- hgt.z ~ 0 + age_0 + age_0.5 + age_1 + age_2 + age_3 +
  (0 + age_0 + age_0.5 + age_1 + age_2 + age_3 | id)
fit <- lmer(f, data, 
            control = lmerControl(check.conv.grad = ctl))

### fitted trajectories for all persons
bse <- t(t(ranef(fit)$id) + fixef(fit))
head(round(bse, 3), 3)

## -----------------------------------------------------------------------------
ctl <- control_brokenstick(
  lmer = lmerControl(check.conv.grad = 
                       .makeCC("warning", tol = 4e-3)))
mod <- brokenstick(hgt.z ~ age | id, data = smocc_200, 
                   knots = c(0, 0.5, 1, 2), boundary = c(0, 3),
                   control = ctl)
head(predict(mod, smocc_200, x = "knots", shape = "wide"), 3)

## ----dataprep-----------------------------------------------------------------
library(brokenstick)
head(smocc_200, 3)

## ----y2z----------------------------------------------------------------------
library(AGD)
z <- with(smocc_200, y2z(y = hgt, 
                         x = age, 
                         sex = ifelse(sex == "male", "M", "F"), 
                         ref = nl4.hgt))
identical(z, smocc_200$hgt.z)

## ----zscores, echo = FALSE, fig.height=3, fig.cap = "Distribution of height SDS for 200 Dutch children.", set.theme = TRUE, warning=FALSE----
ggplot(smocc_200, aes(hgt.z)) +
  xlim(-6, 6) + xlab("Height SDS") +
  theme_light() +
  geom_histogram(binwidth = 0.1, colour = "white", 
                 fill = colors[4], size = 0.1) +
  stat_function(fun = function(x) dnorm(x, 0, 1) * 194,
                color = colors[5], size = 0.7)

## ----z2y----------------------------------------------------------------------
y <- with(smocc_200, z2y(z = hgt.z, 
                         x = age, 
                         sex = ifelse(sex == "male", "M", "F"), 
                         ref = nl4.hgt))
all.equal(y, smocc_200$hgt, tol = 0.0001)

## ----plotsds, echo=FALSE, fig.height=4, fig.cap = "Length growth of 52 infants expressed in the Z-score scale.", set.theme = TRUE, warning=FALSE----
ggplot(smocc_200[1:500, ], aes(x = age, y = hgt.z, group = id, color = as.factor(id))) +
  geom_line(size = 0.2) + 
  geom_point(size = 0.7) +
  scale_colour_viridis_d(option = "viridis") +
  xlab("Age (years)") +
  ylab("Length SDS") +
  theme_light() +
  theme(legend.position = "none")

## ----line1, fig.height=3, fig.cap = "Simple linear model with one line anchored at the extremes.", set.theme = TRUE, warning=FALSE----
fit <- brokenstick(hgt.z ~ age | id, smocc_200)
ids <- c(10001, 10005, 10022)
plot(fit, new_data = data, group = ids, what = "all",
     xlab = "Age (years)", ylab = "Length (cm)")

## ----plotfit2z, fig.height=3, fig.cap = "Broken stick model with two lines.", set.theme = TRUE, warning=FALSE----
fit2 <- brokenstick(hgt.z ~ age | id, smocc_200, knots = 0:2)
plot(fit2, data, group = ids, 
     xlab = "Age (years)", ylab = "Length (SDS)")

## -----------------------------------------------------------------------------
fit2

## ----fit9, cache = TRUE, warning = FALSE, eval = FALSE------------------------
#  knots <- round(c(0, 1, 2, 3, 6, 9, 12, 15, 18, 24)/12, 4)
#  fit9 <- brokenstick(hgt.z ~ age | id, data = smocc_200,
#                      knots = knots, boundary = c(0, 3))

## ----echo = FALSE-------------------------------------------------------------
fit9 <- brokenstick::fit_200

## ----plotfit9, echo = FALSE,fig.height=3, fig.cap = "Broken stick model with nine lines.", set.theme = TRUE, warning=FALSE----
plot(fit9, data, group = ids, xlab = "Age (years)", ylab = "Length (SDS)")

## -----------------------------------------------------------------------------
p1 <- predict(fit9, smocc_200)
head(p1, 3)

## -----------------------------------------------------------------------------
p2 <- predict(fit9, smocc_200, x = "knots")
head(p2, 3)

## -----------------------------------------------------------------------------
p3 <- predict(fit9, smocc_200, x = "knots", strip_data = FALSE)
head(p3, 3)

## ----predx--------------------------------------------------------------------
head(predict(fit9, smocc_200, x = c(0.42, 1.33, 4), shape = "wide"), 3)

## ----pred10001----------------------------------------------------------------
predict(fit9, smocc_200, group = 10001, shape = "vector")

## ----pred10001x---------------------------------------------------------------
tail(predict(fit9, smocc_200, x = c(0.42, 1.33), y = c(-0.5, -1), 
             group = c(10001, 10001), strip_data = FALSE), 3)

## -----------------------------------------------------------------------------
data <- data.frame(
  age = c(0, 0.12, 0.32, 0.62, 1.1, 0.25, 0.46),
  hgt.z = c(-1.2, -1.8, -1.7, -1.9, -2.1, -1.9, -1.5),
  id = c(rep("Fred", 5), rep("Alice", 2)))
p <- predict(fit9, data, x = "knots", strip_data = FALSE)

## ----echo = TRUE, fig.height=3, fig.width=5.2, out.width = "67%", fig.cap="Alice and Fred - observed (blue) and fitted (red) trajectory.", set.theme=TRUE----
plot(fit9, data, ylim = c(-2.5, 0), xlab = "Age (years)", ylab = "Length (SDS)")

## ----echo = FALSE, fig.height=4, fig.cap="Predicted versus observed values."----
pred <- predict(fit_200, smocc_200, shape = "vector")
oldpar <- par(mfrow = c(1, 2))

# Z-score scale
MASS::eqscplot(x = smocc_200$hgt.z, xlab = "Height (SDS)", 
               y = pred, ylab = "Predicted Z-score", 
               pch = ".", xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = "grey")

# cm scale
y_cm <- with(smocc_200, 
             z2y(z = hgt.z, 
                 x = age, 
                 sex = ifelse(sex == "male", "M", "F"), 
                 ref = nl4.hgt))
yhat_cm <- with(smocc_200, 
                z2y(z = pred, 
                    x = age, 
                    sex = ifelse(sex == "male", "M", "F"), 
                    ref = nl4.hgt))
MASS::eqscplot(x = y_cm, xlab = "Height (cm)",
               y = yhat_cm, ylab = "Predicted (cm)", pch = ".")
abline(0, 1, col = "grey")

par <- par(oldpar)

## -----------------------------------------------------------------------------
knots <- round(c(0, 1/3, 1, 2, 4, 6, 10, 14, 24, 29), 3)
labels <- c("birth", "4m", "1y", "2y", "4y", "6y", "10y", "14y", "24y", "")

## ----tbc1, echo=FALSE, fig.height=4, fig.cap = "Body Mass Index (BMI) SDS by log(age + 0.2) (Terneuzen cohort)", set.theme = TRUE, warning=FALSE----
ggplot(mice::tbc, aes(age + 0.2, bmi.z)) +
  theme_light() +
  scale_x_log10(name = "Age", breaks = knots + 0.2, 
                labels = labels, minor_breaks = NULL) +
  scale_y_continuous(name = "BMI SDS", breaks = -4:4, 
                     limits = c(-4, 4)) +
  geom_point(size = 0.5, col = colors[1], shape = 20, na.rm = TRUE)

## ----tbc-lmer, cache=TRUE-----------------------------------------------------
ctl <- control_brokenstick(
  lmer = lmerControl(check.conv.grad = .makeCC("warning", 0.02, NULL),
                     check.conv.singular = .makeCC("ignore", 0.001)))
fit_lmer <- brokenstick(bmi.z ~ age | id, data = mice::tbc,
                        knots = knots, boundary = c(0, 29),
                        control = ctl)

## ----tbc-trajectories-lmer, fig.height=5.4, fig.cap = "Body Mass Index (BMI) SDS trajectories of six subjects, observed (blue) and fitted (red). lmer method.", set.theme = TRUE, warning=FALSE----
ids <- c(8, 1259, 2447, 7019, 7460, 7646)
plot(fit_lmer, mice::tbc, group = ids,
     ylab = "BMI SDS", xlab = "Age (years)")

## ----tbc-kr, cache=TRUE-------------------------------------------------------
fit_kr <- brokenstick(bmi.z ~ age | id, data = mice::tbc,
                      knots = knots, boundary = c(0, 29),
                      method = "kr", seed = 41441)

## ----tbc-trajectories-kr, echo=FALSE, fig.height=5.4, fig.cap = "Body Mass Index (BMI) SDS trajectories of six subjects, observed (blue) and fitted (red). kr method.", set.theme = TRUE, warning=FALSE----
plot(fit_kr, mice::tbc, group = ids, 
     ylab = "BMI SDS", xlab = "Age (years)")

## ----tbc-data-----------------------------------------------------------------
tbc1 <- mice::tbc %>% 
  filter(!is.na(ao) & first) %>% 
  select(id, nocc, sex)
tbc2 <- mice::tbc.target %>% 
  filter(id %in% tbc1$id)
prd <- predict(fit_kr, mice::tbc, x = "knots", 
               shape = "wide", group = tbc1$id)
data <- bind_cols(prd, 
                  select(tbc1, -id), 
                  select(tbc2, -id))
head(data, 3)

## ----tbc-trajectories-all, echo=FALSE, fig.height=4, fig.cap = "Body Mass Index (BMI) SDS trajectories for 92 subjects, coloured by adult overweight status.", set.theme = TRUE, warning=FALSE----
pd <- data %>% 
  select(id, ao, `0`:`24`) %>% 
  tidyr::pivot_longer(`0`:`24`, names_to = "age", values_to = "bmi") %>% 
  mutate(age = as.numeric(age), 
         logage = log10(age + 0.2))
ggplot(pd, aes(x = age + 0.2, y = bmi, group = id, colour = factor(ao))) +
  theme_light() +
  theme(legend.position = "none") +
  geom_line(data = subset(pd, ao == 0)) +
  geom_line(data = subset(pd, ao == 1)) +
  scale_x_log10(name = "Age (log10(year + 0.2))", breaks = knots + 0.2, 
                labels = labels, minor_breaks = NULL) +
  scale_y_continuous(name = "BMI SDS", breaks = -4:4, 
                     limits = c(-4, 3)) +
  scale_colour_manual(values = c("grey", colors[5]))

## -----------------------------------------------------------------------------
m1 <- lm(bmi.z.jv ~ `6`, data)
m2 <- lm(bmi.z.jv ~ `6` + I(`6`-`4`), data)
anova(m1, m2)

## ----t2t, cache=TRUE----------------------------------------------------------
fit <- brokenstick(hgt.z ~ age | id, data = smocc_200, 
                   knots = 1:4/2, boundary = c(0, 3))
t2t <- fit$omega + diag(fit$sigma2, ncol(fit$omega))
round(cov2cor(t2t), 2)

## ----t2tm, cache=TRUE---------------------------------------------------------
fit <- brokenstick(hgt.z ~ age | id, data = smocc_200, 
                   knots = seq(0, 2, 0.1), boundary = c(0, 3),
                   method = "kr")
t2t <- fit$omega + diag(fit$sigma2, ncol(fit$omega))
dim(t2t)

## ----weightloss-data, fig.height=4.5, echo=FALSE, fig.cap = "Daily body weight (KG) for 12 subjects under three conditions.", set.theme = TRUE, warning=FALSE----
data <- brokenstick::weightloss
ggplot(data, aes(day, body_weight, group = subject, 
                 colour = condition)) +
  scale_x_continuous(name = "Day", breaks = c(0, 21, 42, 63), minor_breaks = c(7, 14, 28, 35, 49, 56)) +
  ylab("Body weight (KG)") +
  geom_line() + geom_point(size = 0.7) +
  theme_light() +
  theme(legend.position = "bottom")

## ----fig5, echo=FALSE, eval=FALSE---------------------------------------------
#  ggplot(data, aes(factor(week), body_weight)) +
#    geom_boxplot() +
#    facet_wrap(vars(subject), scales = "free") +
#    theme_light()

## ----weightloss-constant, fig.height=5, fig.cap = "Constant model. Observed and fitted trajectories for a model that summarises each experimental period by a constant.", set.theme = TRUE, warning=FALSE----
fit0 <- brokenstick(body_weight ~ day | subject, data, 
                    knots = c(0, 21, 42, 63), degree = 0)
plot(fit0, data, size_y = 0, color_y = rep("grey", 2), what = "all",
     scales = "free_y", xlab = "Day", ylab = "Body weight (KG)", 
     n_plot = 12, ncol = 4)

## -----------------------------------------------------------------------------
prd <- data.frame(predict(fit0, data, x = "knots", shape = "wide"))
control <- prd[, 2]
diet <- prd[, 3]
diet[c(4, 12)] <- prd[c(4, 12), 4]
activity <- prd[, 4]
activity[c(4, 12)] <- prd[c(4, 12), 3]
round(data.frame(diet_control = diet - control, 
                 activity_control = activity - control, 
                 activity_diet = activity - diet), 1)

## ----weightloss-slope, fig.height=5, fig.cap = "Broken stick model. Observed and fitted trajectories for a model that summarises each experimental period by a line.", set.theme = TRUE, warning=FALSE----
fit1 <- brokenstick(body_weight ~ day | subject, data, 
                    knots = c(0, 21, 42, 63))
plot(fit1, data, size_y = 0, color_y = rep("grey", 2), what = "all",
     size_yhat = 1.5, scales = "free_y", , xlab = "Day", ylab = "Body weight (KG)",
     n_plot = 12, ncol = 4)

## -----------------------------------------------------------------------------
prd <- data.frame(predict(fit1, data, x = "knots", shape = "wide"))
control <- prd[, 3] - prd[, 2]
diet <- prd[, 4] - prd[, 3]
diet[c(4, 12)] <- prd[c(4, 12), 5] - prd[c(4, 12), 4]
activity <- prd[, 5] - prd[, 4]
activity[c(4, 12)] <- prd[c(4, 12), 4] - prd[c(4, 12), 3]
round(data.frame(control = control, 
                 diet = diet, 
                 activity = activity), 1)

## -----------------------------------------------------------------------------
df <- data.frame(y = c(diet, activity),
                 act = rep(c(0, 1), each = 12),
                 per2 = rep(c(rep(1, 3), 0, rep(1, 7), 0), 2))
coef(lm(y ~ act, data = df))
coef(lm(y ~ act + per2, data = df))

## ----zscores1, echo=TRUE------------------------------------------------------
boy <- data.frame(x = c(1/12, 14/12), y = c(52.6, 81.7))
ref <- AGD::nl4.hgt
boy$z <- AGD::y2z(y = boy$y, x = boy$x, sex = "M", ref = ref)
boy$z

## ----inter1, echo = FALSE, fig.height = 4, fig.cap = "Linear interpolation in the cm scale results in an unrealistic trajectory at intermediate ages.", set.palette = TRUE----
sds <- c(-2, -1, 0, 1, 2)
age <- round(seq(0, 1.25, 1/48), 3)
z <- rep(sds, times = length(age))
x <- rep(age, each = length(sds))
w <- z2y(z = z, x = x, sex = 'M', ref = ref)
w <- matrix(w, ncol = length(sds), byrow = TRUE)
dimnames(w) <- list(age, sds)

oldpar <- par(mfrow = c(1, 2))
matplot(x = as.numeric(rownames(w)), y = w, type = "l", lty = 1, 
        col = "grey", lwd = c(1, 1, 1.5, 1, 1),
        xlab = "Age (years)", ylab = "Length (cm)")
matpoints(x = boy$x, y = boy$y, pch = 20, 
          cex = 1.5, col = colors[1], type = "o")
v <- matrix(c(c(0, 1.25), rep(-2:2, each = 2)), 
            ncol = 6, byrow = FALSE)
matplot(x = v[,1], y = v[,2:6], 
        type = "l", lty = 1,
        col = "grey", lwd = c(1, 1, 1.5, 1, 1),
        ylim = c(-2.5, 2.5),
        xlab = "Age (years)", ylab = "Length (SDS)")
yout <- approx(x = boy$x, y = boy$y, xout = as.numeric(rownames(w)))$y
zout <- y2z(x = as.numeric(rownames(w)), y = yout, ref = ref)
matpoints(x = boy$x, y = boy$z, pch = 20, 
          cex = 1.5, col = colors[1], type = "p")
matpoints(x = as.numeric(rownames(w)), y = zout, pch = 20, 
          cex = 1.5, col = colors[1], type = "l")

## ----inter2, echo = FALSE, fig.height = 4, fig.cap = "Linear interpolation in the Z-score scale results in a more realistic trajectory at intermediate ages.", set.palette = TRUE----
zout <- approx(x = boy$x, y = boy$z, xout = as.numeric(rownames(w)))$y
yout <- z2y(x = as.numeric(rownames(w)), z = zout, ref = ref)

oldpar <- par(mfrow = c(1, 2))
matplot(x = as.numeric(rownames(w)), y = w, type = "l", lty = 1, 
        col = "grey", lwd = c(1, 1, 1.5, 1, 1),
        xlab = "Age (years)", ylab = "Length (cm)")
matpoints(x = boy$x, y = boy$y, pch = 20, 
          cex = 1.5, col = colors[1], type = "p")
matpoints(x = as.numeric(rownames(w)), y = yout, pch = 20, 
          cex = 1.5, col = colors[1], type = "l")

v <- matrix(c(c(0, 1.25), rep(-2:2, each = 2)), 
            ncol = 6, byrow = FALSE)
matplot(x = v[,1], y = v[,2:6], 
        type = "l", lty = 1,
        col = "grey", lwd = c(1, 1, 1.5, 1, 1),
        ylim = c(-2.5, 2.5),
        xlab = "Age (years)", ylab = "Length (SDS)")
matpoints(x = boy$x, y = boy$z, pch = 20, 
          cex = 1.5, col = colors[1], type = "p")
matpoints(x = as.numeric(rownames(w)), y = zout, pch = 20, 
          cex = 1.5, col = colors[1], type = "l")

## ----brokenstickinterpolation, echo = TRUE------------------------------------
# prepare data input
age <- round(seq(2/24, 28/24, 1/24), 3)
z <- rep(NA, length(age))
z[1] <- boy$z[1]; z[length(z)] <- boy$z[2]

# predict with broken stick model
zout <- predict(fit_200, x = age, y = z, shape = "vector")

# convert predicted values to Y-scale
yout <- z2y(x = age, z = zout, ref = ref)

## ----inter3, echo = FALSE, fig.height = 4, fig.cap = "Broken stick model fitted in the SDS scale results in a most realistic expected trajectory at intermediate ages.", set.palette = TRUE----
# create the plots
oldpar <- par(mfrow = c(1, 2))
matplot(x = as.numeric(rownames(w)), y = w, type = "l", lty = 1, 
        col = "grey", lwd = c(1, 1, 1.5, 1, 1),
        xlab = "Age (years)", ylab = "Length (cm)")
matpoints(x = age, y = yout, pch = 20, 
          cex = 1.5, col = colors[1], type = "l")
matpoints(x = boy$x, y = boy$y, pch = 20, 
          cex = 1.5, col = colors[1], type = "p")

v <- matrix(c(c(0, 1.25), rep(-2:2, each = 2)), 
            ncol = 6, byrow = FALSE)
matplot(x = v[,1], y = v[,2:6], 
        type = "l", lty = 1,
        col = "grey", lwd = c(1, 1, 1.5, 1, 1),
        ylim = c(-2.5, 2.5),
        xlab = "Age (years)", ylab = "Length (SDS)")
matpoints(x = age, y = zout, col = colors[2], type = "l")
matpoints(x = age, y = zout, pch = 20, 
          cex = 0.3, col = colors[2], type = "p")
matpoints(x = boy$x, y = boy$z, pch = 20, 
          cex = 1.5, col = colors[1], type = "p")
par(oldpar)

## ----mi-code, cache=TRUE, fig.height=3, fig.cap = "Observed data plotted on top of 20 imputed trajectories.", set.theme = TRUE, warning=FALSE----
ctl <- control_brokenstick(kr = control_kr(imp_skip = 10))
knots <- round(c(0, 1, 2, 3, 6, 9, 12, 15, 18, 24)/12, 4)
data <- bind_rows(smocc_200[!is.na(smocc_200$hgt.z), ],
                  expand.grid(id = unique(smocc_200$id), age = knots))
fit_kr <- brokenstick(hgt.z ~ age | id, data = data,
                      knots = knots, boundary = c(0, 3),
                      method = "kr", control = ctl, seed = 15244)
plot(fit_kr, new_data = data, show = c(TRUE, FALSE, TRUE),
     group = c(10001, 10005, 10022),
     xlab = "Age (years)", ylab = "Length (SDS)")

## ----t2t-imp------------------------------------------------------------------
expand.grid(id = unique(smocc_200$id), age = knots) %>%
  bind_cols(as.data.frame(fit_kr$draws)) %>%
  pivot_longer(cols = starts_with("V"), names_to = "imp") %>%
  pivot_wider(id_cols = c("id", "imp"), names_from = "age") %>%
  select(-id, -imp) %>%
  cor() %>%
  round(2)

## ----echo=FALSE, fig.height=5, fig.cap = "Curve matching. Predict infant length at 14 months given length data up to 6 months using 10 matches.", set.theme = TRUE, warning=FALSE----
knitr::include_graphics("figures/JAMES_tryout.png")

## ----cm1, cache=TRUE----------------------------------------------------------
donor_data <- smocc_200 %>% 
  filter(id != "10001")
target_data <- smocc_200 %>% 
  filter(id == "10001" & age < 0.51)

# fit brokenstick model at time level
knots <- round(c(0, 1, 2, 3, 6, 9, 12, 15, 18, 24)/12, 4)
fit <- brokenstick(hgt.z ~ age | id, data = donor_data,
                   knots = knots, boundary = c(0, 3),
                   method = "kr", seed = 15244)

## ----cm2----------------------------------------------------------------------
# predict with matching model at child level
covariates <- donor_data %>% 
  group_by(id) %>% 
  slice(1)
bse <- predict(fit, donor_data, x = "knots", shape = "wide")
donors <- bind_cols(covariates, select(bse, -id))
model <- lm(`1.25` ~ `0` + `0.0833` + `0.1667` + `0.25` + `0.5`
            + sex + ga + bw, data = donors)
summary(model)

## ----cm3----------------------------------------------------------------------
donors_pred <- predict(model)
names(donors_pred) <- donors$id

target <- bind_cols(
  slice(target_data, 1), 
  select(predict(fit, target_data, x = "knots", shape = "wide"), -id))
target_pred <- predict(model, newdata = target)

matches <- sort(abs(donors_pred - target_pred))[1:10]
matches

## ----cm4, fig.height=4, fig.cap = "Curve matching. Observed and fitted trajectories of 10 matches for subject 10001.", set.theme = TRUE, warning=FALSE----
ids <- as.numeric(names(matches))
plot(fit, new_data = donor_data, group = ids, 
     xlim = c(0, 1.4), size_y = 1, size_yhat = 0,
     xlab = "Age (years)", ylab = "Length (SDS)", 
     ncol = 5)

