# Lee-Carter Style Forecast of Age-Specific Annual Deaths

# Init ------------------------------------------------------------

library(yaml)
library(readxl)
library(tidyverse)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  global = './src/00-global_objects.R',
  artsakh = './dat/Artsakh.xlsx'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  data = './dat/artsakh_predictions.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# global objects and functions
global <- source(paths$input$global)

# constants specific to this analysis
cnst <- within(list(), {
  n_age = length(config$skeleton$age)
  n_year = length(config$skeleton$year)
  n_sex = 2
  fitting_period = c('2015', '2016', '2017', '2018', '2019')
  cv_training_period = c('2015', '2016', '2017', '2018')
  cv_test_period = c('2019')
})

# list containers for analysis artifacts
artsakh <- list()

# Import data -----------------------------------------------------

artsakh$raw <- read_excel(paths$input$artsakh)

# Transform to array ----------------------------------------------

artsakh$array <-
  array(
    NA, dim = c(cnst$n_age, cnst$n_year, cnst$n_sex),
    dimnames = list(
      'age' = config$skeleton$age,
      'year' = config$skeleton$year,
      'sex' = config$skeleton$sex
    )
  )

artsakh$array[,,'Male'] <-
  unlist(artsakh$raw[,grep('^male', names(artsakh$raw))])
artsakh$array[,,'Female'] <-
  unlist(artsakh$raw[,grep('^female', names(artsakh$raw))])

# Lee-carter style model ------------------------------------------

# Predict annual age-specific death counts n steps ahead
# X: matrix of age specific death counts by year
# npred: number of years to predict
# nsim: number of simulation draws for prediction intervals
# p: coverage of prediction intervals
LeeCarterPrediction <- function (X, npred = 1, nsim = 100, p = 0.95) {
  
  # we fit the model
  # log(D(x,t) + 1) = a(x) + b(x)*k(t) + Norm(0, sigma)
  # via singular value decomposition of age (x) by period (t) matrix
  # of logged death counts D(x,t).
  # we use log(x+1) to avoid issues with 0's.
  log1p_Dxt <- log1p(X)
  ax <- rowMeans(log1p_Dxt)
  log1p_Dxt_normalized <- sweep(log1p_Dxt, 1, ax, FUN = '-')
  SVD <- svd(log1p_Dxt_normalized)
  
  # highest singular value
  s <- SVD$d[1]
  # Lee-Carter k(t) index
  # annual index of log death count level
  kt <- SVD$v[,1]
  # Lee-Carter b(x) coefficients
  # rates of change in age specific log death counts
  # wrt. changes in overall death counts level k(t)
  bx <- SVD$u[,1]
  
  # random walk with drift forecast of k(t) death count index
  kt_arima <- arima(kt, order = c(0, 1, 0), include.mean = TRUE)
  kt_forecast <- predict(kt_arima, n.ahead = npred, se.fit = TRUE)
  kt_avg <- c(kt, kt_forecast$pred)
  kt_se <- c(rep(0, length.out = length(kt)), kt_forecast$se)
  
  # simulate draws from kt predictive distribution
  N = nsim
  kt_sim <- rbind(
    kt_avg,
    mapply(
      FUN = function (kt_avg, kt_se) {
        rnorm(N, mean = kt_avg, sd = kt_se)
      }, kt_avg, kt_se
    ))
  
  # the outer product of the b(x) and k(t) vectors
  # gives us the expected log death counts after
  # adding back the previously subtracted row-means.
  # we remove the log transform.
  Dxt_expected_sim <-
    array(NA, dim = c(length(bx), length(kt_avg), N+1))
  for (i in 1:(N+1)) {
    Dxt_expected_sim[,,i] <-
      expm1(s*outer(bx, kt_sim[i,])+ax)
    Dxt_expected_sim[,,i][Dxt_expected_sim[,,i]<0] <- 0
  }
  
  # add Poisson variation around expected death counts
  Dxt_predicted_sim <-
    apply(Dxt_expected_sim[,,-1],
          1:3, function(x) rpois(n = 1, lambda = x))
  
  # quantiles of prediction intervals
  Dxt_predicted_qlo <-
    apply(
      Dxt_predicted_sim, MARGIN = 1:2, quantile,
      probs = (1-p)/2, type = 1, na.rm = TRUE
    )
  Dxt_predicted_qhi <-
    apply(
      Dxt_predicted_sim, MARGIN = 1:2, quantile,
      probs = (1+p)/2, type = 1, na.rm = TRUE
    )
  
  predictions <- list(
    expected = Dxt_expected_sim[,,1],
    lo = Dxt_predicted_qlo,
    hi = Dxt_predicted_qhi
  )
  
  return(predictions)
  
}

artsakh$lc_male <- LeeCarterPrediction(
  artsakh$array[,cnst$fitting_period,'Male'],
  nsim = 100, p = 0.95
)

artsakh$lc_female <- LeeCarterPrediction(
  artsakh$array[,cnst$fitting_period,'Female'],
  nsim = 100, p = 0.95
)

# Predicted vs. observed ------------------------------------------

artsakh$skeleton <- 
  expand.grid(
    age = unlist(config$skeleton$age),
    year = config$skeleton$year,
    sex = unlist(config$skeleton$sex)
  ) %>%
  mutate(
    observed = c(artsakh$array),
    predicted = c(artsakh$lc_male$expected, artsakh$lc_female$expected),
    lo = c(artsakh$lc_male$lo, artsakh$lc_female$lo),
    hi = c(artsakh$lc_male$hi, artsakh$lc_female$hi)
  )

artsakh$skeleton %>%
  ggplot() +
  aes(x = year, color = sex, group = sex) +
  geom_point(
    aes(y = observed, shape = sex),
  ) +
  geom_crossbar(
    aes(ymin = lo, y = predicted, ymax = hi),
    filter(artsakh$skeleton, year == 2020)
  ) +
  scale_color_manual(
    values = c(`Male` = config$figspec$colors$sex$male,
               `Female` = config$figspec$colors$sex$female)
  ) +
  scale_shape_manual(values = c(`Male` = 1, `Female` = 4)) +
  facet_wrap(~age, scales = 'free_y') +
  labs(y = 'Deaths', x = 'Year', color = NULL) +
  MyGGplotTheme(panel_border = TRUE, grid = '')

# Cross validation ------------------------------------------------

artsakh$skeleton_cv <- 
  expand.grid(
    age = unlist(config$skeleton$age),
    year = config$skeleton$year,
    sex = unlist(config$skeleton$sex)
  ) %>%
  mutate(
    observed = c(artsakh$array),
    sample = ifelse(year < 2019, 'training', 'test')
  ) %>%
  filter(year < 2020)

# Poisson GLM prediction for 2019
artsakh$poisson_glm <-
  CountGAM(
    artsakh$skeleton_cv,
    observed ~ year*age, poisson(link = "log"),
    col_year = year, nsim = 500, col_stratum = 'sex',
    col_sample = 'sample')

# Lee-Carter prediction for 2019
artsakh$lc_cv_male <- LeeCarterPrediction(
  artsakh$array[,cnst$cv_training_period,'Male'],
  nsim = 500, p = 0.95, npred = 1
)
artsakh$lc_cv_female <- LeeCarterPrediction(
  artsakh$array[,cnst$cv_training_period,'Female'],
  nsim = 500, p = 0.95, npred = 1
)
 
artsakh$cv_predictions <-
  artsakh$skeleton_cv %>%
  mutate(
    predicted_lc = c(artsakh$lc_cv_male$expected, artsakh$lc_cv_female$expected),
    lo_lc = c(artsakh$lc_cv_male$lo, artsakh$lc_cv_female$lo),
    hi_lc = c(artsakh$lc_cv_male$hi, artsakh$lc_cv_female$hi),
    predicted_poisson = c(artsakh$poisson_glm$predicted),
    lo_poisson =
      apply(
        artsakh$poisson_glm[,grep('^simulated', names(artsakh$poisson_glm))],
        1, quantile, p = 0.025
      ),
    hi_poisson = apply(
      artsakh$poisson_glm[,grep('^simulated', names(artsakh$poisson_glm))],
      1, quantile, p = 0.975
    )
  ) %>%
  as_tibble()

# cross validation statistics
artsakh$cv_predictions %>%
  filter(year == 2019) %>%
  pivot_longer(cols = c(predicted_lc, lo_lc, hi_lc,
                        predicted_poisson, lo_poisson, hi_poisson)) %>%
  separate(name, into = c('measure', 'model')) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  group_by(model) %>%
  summarise(
    # weighted mean percentage error
    wmpe = sum(observed - predicted) / sum(observed),
    # weighted mean absolute percentage error
    wmape = sum(abs(observed - predicted)) / sum(observed),
    # coverage of prediction interval
    coverage = sum(observed > lo & observed < hi) / n()
  )

# Lee-Carter does not suffer from the exploding prediction
# interval in age groups with many 0's as does Poisson
artsakh$cv_predictions %>%
  filter(year == 2019, sex == 'Female', age == '20 - 24') %>%
  glimpse()
  
# Export ----------------------------------------------------------

# export results of analysis
