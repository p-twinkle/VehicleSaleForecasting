## ------------------------------------------
# Project :  Time Series Forecasting for RPM
# Author  :  Twinkle Panda
# Date.   :  24th April 2025
## ------------------------------------------

## 0. Load required packages ----------------------------------------------------
library(readr)       # for reading CSV files
library(tsibble)     # for time series tibbles
library(lubridate)   # for date handling
library(dplyr)       # for data manipulation
library(ggplot2)     # for plotting
library(fable)       # for modeling and forecasting
library(feasts)      # for decomposition and ACF/PACF
library(scales)      # for comma_format()
library(tidyr)
library(nortest)
library(stlplus)
library(fabletools)

rm(list=ls()) # clear all objects from the current workspace/environment
options(tibble.width = Inf) # print all columns of a tibble

# setwd("/Users/twinklepanda/Desktop/Me/MSBA/4.3 Spring/TimeSeries/Project")

## 1. Read data and convert to tsibble -----------------------------------------
rpm_df <- read_csv("RPM.csv") %>%
  mutate(observation_date = yearmonth(observation_date)) %>%
  as_tsibble(index = observation_date)

## 2. RPM Over Time ------------------------------------------------------------
autoplot(rpm_df, RPM)

# Pretty Plot
rpm_df %>%
  autoplot(RPM, colour = "darkblue", size = 0.55) +
  # highlight the COVID intervention period
  annotate("rect",
           xmin = yearmonth("2020 Jan"),
           xmax = yearmonth("2021 Jan"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2,
           fill = "red") +
  labs(
    title    = "Monthly Revenue Passenger Miles (RPM)",
    subtitle = "Jan 2000 – Dec 2024",
    x        = NULL,
    y        = "RPM (billions)",
    # caption  = "Data source: US Bureau of Transportation Statistics"
  ) +
  scale_y_continuous(
    labels = comma_format(scale = 1e-6, suffix = "M"), # show in millions
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  theme_light(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(size = 12, margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

## 3. LogRPM -----------------------------------------------------------------------
rpm_df <- rpm_df %>%
  mutate(LogRPM = log(RPM))

rpm_df %>%
  autoplot(LogRPM, colour = "darkblue", size = 0.55) +
  # highlight the COVID intervention period
  annotate("rect",
           xmin = yearmonth("2020 Jan"),
           xmax = yearmonth("2021 Jan"),
           ymin = -Inf, ymax = Inf,
           alpha = 0.2,
           fill = "red") +
  labs(
    title    = "Log of Monthly Revenue Passenger Miles (LogRPM)",
    subtitle = "Jan 2000 – Dec 2024",
    x        = NULL,
    y        = "RPM (billions)",
    # caption  = "Data source: US Bureau of Transportation Statistics"
  ) +
  scale_y_continuous(
    labels = comma_format(scale = 1e-6, suffix = "M"), # show in millions
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  theme_light(base_size = 14) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(size = 12, margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

## 4. STL decomposition on log(RPM) --------------------------------------------
rpm_stl <- rpm_df %>%
  model(
    # 7-month window LOESS for seasonality, robust to outliers
    STL(LogRPM ~ season(window = 7), robust = TRUE)
  ) %>%
  components()

# Plot the four STL components
rpm_stl_long <- rpm_stl %>%
  select(observation_date, LogRPM, trend, season_year, remainder) %>%
  pivot_longer(
    -observation_date,
    names_to  = "component",
    values_to = "value"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c("LogRPM", "trend", "season_year", "remainder"),
      labels = c(
        "Data (log RPM)",
        "Trend",
        "Seasonal",
        "Remainder"
      )
    )
  )

ggplot(rpm_stl_long, aes(x = observation_date, y = value)) +
  geom_line(color = "darkblue") +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  scale_x_yearmonth(
    date_breaks = "5 years",
    date_labels = "%Y",
    expand      = expansion(add = c(0, 0))
  ) +
  labs(
    title    = "STL Decomposition of log(RPM)",
    subtitle = "Jan 2000 – Dec 2024",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text        = element_text(face = "bold", size = 12),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray90")
  )

## 5. Stationarity Checks ----------------------------------------------------------
rpm_df <- rpm_df %>%
  # compute log and bring in season_year from the STL components
  left_join(
    rpm_stl %>% select(observation_date, season_year),
    by = "observation_date"
  ) %>%
  mutate(
    LogRPM        = log(RPM),               # log scale
    LogRPM_sa     = LogRPM - season_year,   # seasonally adjusted
    diff1         = difference(LogRPM_sa),  # first difference
    diff12        = difference(LogRPM_sa, 12)# seasonal difference
  )

# Plot first difference of LogRPM_sa
mean_diff1 <- mean(rpm_df$diff1, na.rm=TRUE)

rpm_df %>% autoplot(diff1) +
  geom_point(aes(y=diff1)) +
  geom_hline(aes(yintercept=mean_diff1), lty=2) +
  theme_classic() +
  ggtitle("First differences of LogRPM_sa vs. Time") + xlab("Time") + ylab("First differences of LogRPM_sa")

# Plot LogRPM_sa
autoplot(rpm_df, LogRPM_sa)

# KPSS tests
unitroot_kpss(rpm_df$LogRPM_sa) # kpss_pvalue: 0.0787
unitroot_kpss(rpm_df$diff1) # kpss_pvalue: 0.1

# ACF plots
rpm_df %>% ACF(LogRPM_sa) %>% 
  autoplot() + ggtitle("ACF: seasonally adjusted log(RPM)")

rpm_df %>% ACF(diff1) %>% 
  autoplot() + ggtitle("ACF: LogRPM_sa First Difference")


## 6. Forecast LogRPM using model suggested by ARIMA(0,0,0) --------------------

# Let ARIMA select a model to forecast LogRPM
result_ARIMA_LogRPM <- rpm_df %>%
  model(ARIMA(LogRPM ~ PDQ(0,0,0),
              stepwise=FALSE, approximation=FALSE, trace=TRUE))
report(result_ARIMA_LogRPM) 

# ACF of residuals
result_ARIMA_LogRPM %>% augment() %>% ACF(.resid) %>% autoplot()

# Augment
result_ARIMA_Y_augment <- augment(result_ARIMA_LogRPM)

# Plot the actual data (LogRPM) and the 4-periods ahead forecast
result_ARIMA_Y_forecast <- result_ARIMA_LogRPM %>% forecast(h=4)

# Plot in-sample and out-of-sample forecasts
result_ARIMA_Y_forecast %>% autoplot(rpm_df) +
  geom_line(aes(y = result_ARIMA_Y_augment$.fitted), color = "red", lty = 2) +
  ggtitle("Black represents LogRPM; Red represents 4-steps ahead in-sample forecasts)") +
  xlab("Time") + ylab("LogRPM") +
  theme_classic()

# Standard error
forecast_errors <- result_ARIMA_Y_augment$.resid
sd_forecast_errors <- sd(forecast_errors, na.rm = TRUE)
print(sd_forecast_errors) # standard error: 0.196965

# Residual Plots
gg_tsresiduals(result_ARIMA_LogRPM)

# AD test
ad.test(result_ARIMA_Y_augment$.resid)

## 7. Forecast LogRPM_sa using model suggested by ARIMA(0,0,0) -----------------

# Let ARIMA select a model to forecast LogRPM_sa
result_ARIMA_LogA <- rpm_df %>%
  model(ARIMA(LogRPM_sa ~ PDQ(0,0,0),
              stepwise=FALSE, approximation=FALSE, trace=TRUE))
report(result_ARIMA_LogA) 

# ACF of residuals
result_ARIMA_LogA %>% augment() %>% ACF(.resid) %>% autoplot()

# Augment
result_ARIMA_Y_augment <- augment(result_ARIMA_LogA)

# Plot the actual data (LogRPM_SA) and the 4-periods ahead forecast
result_ARIMA_Y_forecast <- result_ARIMA_LogA %>% forecast(h=4)

# Plot in-sample forecasts
result_ARIMA_Y_forecast %>% autoplot(rpm_df) +
  geom_line(aes(y = result_ARIMA_Y_augment$.fitted), color = "red", lty = 2) +
  ggtitle("Black represents LogRPM_SA; Red represents 4-steps ahead in-sample forecasts)") +
  xlab("Time") + ylab("LogRPM_SA") +
  theme_classic()

# Standard error
forecast_errors <- result_ARIMA_Y_augment$.resid
sd_forecast_errors <- sd(forecast_errors, na.rm = TRUE)
print(sd_forecast_errors) # standard error: 0.1709829

# Residual Plots
gg_tsresiduals(result_ARIMA_LogA)

# AD test
ad.test(result_ARIMA_Y_augment$.resid)

## 8. Train Test Split ---------------------------------------------------------

rm(list=ls()) # clear all objects from the current workspace/environment
options(tibble.width = Inf) # print all columns of a tibble

rpm_raw <- read_csv("RPM.csv") %>%
  mutate(observation_date = yearmonth(observation_date)) %>%
  as_tsibble(index = observation_date)

# Recompute STL on the full series to get the seasonal pattern for EVERY month
full_stl <- rpm_raw %>%
  mutate(LogRPM = log(RPM)) %>%
  model(
    STL(LogRPM ~ season(window = 7), robust = TRUE)
  ) %>%
  components() %>%
  select(observation_date, season_full = season_year)

# Attach that full‐series seasonal component back to rpm_raw
rpm_full <- rpm_raw %>%
  mutate(LogRPM = log(RPM)) %>%
  left_join(full_stl, by = "observation_date") %>%
  mutate(
    LogRPM_sa_full = LogRPM - season_full
  )

# Now split rpm_full INTO train & test, but keep LogRPM_sa_full
total_rows <- nrow(rpm_full)
n_test     <- 12

rpm_train  <- rpm_full %>% slice_head(n = total_rows - n_test)
rpm_test   <- rpm_full %>% slice_tail(n = n_test)

## 9. Covid Intervention -------------------------------------------------------

# Create a working copy of rpm_train and blank out COVID months
rpm_train_int <- rpm_train %>%
  mutate(
    RPM_i      = if_else(
      observation_date >= yearmonth("2020 Mar") &
        observation_date <= yearmonth("2022 Feb"),
      NA_real_, 
      RPM
    ),
    LogRPM_i   = log(RPM_i)
  )

# Run STL+ to get components across the NA gap (s.window=7 LOESS)
ts_raw    <- ts(rpm_train_int$LogRPM_i, frequency = 12)
comps_int <- stlplus(ts_raw, s.window = 7, robust = TRUE)

# Attach the interpolated seasonal component and form LogRPM_sa_int
rpm_train_int <- rpm_train_int %>%
  mutate(
    season_i      = comps_int$data$seasonal,
    LogRPM_sa_int = comps_int$data$raw - comps_int$data$seasonal
  )

# Check: Plot the seasonally-adjusted series, including the stitched gap
ggplot(rpm_train_int, aes(observation_date)) +
  geom_line(aes(y = LogRPM_i), color = "grey70") +
  geom_line(aes(y = season_i),    color = "blue",    lty = 2) +
  geom_line(aes(y = LogRPM_sa_int),color = "darkgreen") +
  labs(
    title = "Pandemic adjusted STL Components",
    subtitle = "Grey=raw log, Blue=seasonal, Green=seasonally adjusted",
    y     = "log(RPM)"
  ) +
  theme_minimal()

## 10. Fit 2 competing models --------------------------------------------------

# Fit two models on LogRPM_sa_int
# model 1
mod1 <- rpm_train_int %>%
  model(ARIMA(LogRPM_sa_int ~ PDQ(0,0,0),
              stepwise=FALSE, approximation=FALSE, trace=TRUE))
report(mod1) 

mod1 %>% gg_tsresiduals()

# model 2
mod2 <- rpm_train_int %>%
  model(
    SeasonalMA = ARIMA(
      LogRPM_sa_int ~ pdq(0,1,1)   # non‐seasonal: d=1, q=1
      + PDQ(0,0,1),  # seasonal: D=0, Q=1 at lag 12
      stepwise      = FALSE,
      approximation = FALSE
    )
  )
report(mod2)
mod2 %>% gg_tsresiduals()

## 11. Bootstrap 12-month forecasts ------------------------------------------

n_boot <- 100000

# Forecast from models
fc1_raw <- mod1 %>%
  forecast(h = 12, bootstrap = TRUE, times = n_boot)
fc2_raw <- mod2 %>%
  forecast(h = 12, bootstrap = TRUE, times = n_boot)

## 12. Prepare hold-out with seasonally-adjusted log -------------------------

# rpm_test comes from rpm_full and already has LogRPM_sa_full
rpm_test_int <- rpm_test %>%
  transmute(
    observation_date,
    LogRPM_sa_int = LogRPM_sa_full
  )

## 13. Point-forecast accuracy ----------------------------------------------

acc1 <- accuracy(fc1_raw, rpm_test_int)
acc2 <- accuracy(fc2_raw, rpm_test_int)

bind_rows(
  Base       = acc1  %>% select(.model, ME, RMSE, MAE),
  SeasonalMA = acc2  %>% select(.model, ME, RMSE, MAE)
) %>%
  print()

## 14. 95% interval coverage -----------------------------------------------

pi1 <- hilo(fc1_raw, level = 95)
pi2 <- hilo(fc2_raw, level = 95)

lower1 <- pi1$`95%`$lower
upper1 <- pi1$`95%`$upper

lower2 <- pi2$`95%`$lower
upper2 <- pi2$`95%`$upper

# 3) Grab the actual test values in the same order
obs <- rpm_test_int$LogRPM_sa_int

# 4) Compute empirical coverage
cov1 <- mean(obs >= lower1 & obs <= upper1, na.rm = TRUE)
cov2 <- mean(obs >= lower2 & obs <= upper2, na.rm = TRUE)

tibble(
  model    = c("Base", "SeasonalMA"),
  coverage = c(cov1, cov2)
) %>% print()


## ---- END ----

