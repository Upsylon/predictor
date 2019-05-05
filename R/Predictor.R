#' @title Generates a baseline based on historical data
#' @description This function enables to generate a baseline based on historical
#' sales data and to plot the graphs linked to the different step of the
#' process to have a deep view on a given series. It eliminates the role of external impactor that have a significant
#' effect on the sales to obtain a smooth baseline.
#' @param sales_data A vector containing historical sales data
#' @param promo_done A logical variable specifying if promotions are done for the product.
#' @param sizeroll An odd integer that determine the width of the rolling median
#' and rolling average to take.
#' @param smoother The smoother that should be considered. It can be "mean", "median" or "loess".
#' @param showgraph A logical value. If TRUE, returns the different
#' graphs obtained after each step of the process.
#' @return A vector containing the baseline.
#' @return A vector composed of binary variables that indicate whether an external
#' impactor was detected in a given period.
#' @return In option, the different graphs obtained in the process if showgraph is TRUE.
#' @author Grandadam Patrik
#' @importFrom magrittr %>%
#' @export
#' @examples
#' baseline() ###################### A COMPLETER ###############

baseline <- function(sales_data,
                     promo_done = FALSE,
                     sizeroll = 11,
                     smoother = "mean",
                     showgraph = FALSE) {
  
  sales_data <- data.frame(sales_data)
  week <- c(1:nrow(sales_data))
  
  if (promo_done == TRUE) {
    
    my_baseline <-  vector()
    proba_baseline <- vector()
    limit_baseline <- vector()
    
    mydata_baseline <- sales_data
    index_baseline <- !is.na(sales_data)
    mydata_baseline[is.na(mydata_baseline)] <- 0
    
    rollmed_baseline <- zoo::rollapply(mydata_baseline,
                                       width = sizeroll,
                                       FUN = median,
                                       fill = NA,
                                       partial = TRUE,
                                       coredata = FALSE) %>%
      as.data.frame
    
    myres_baseline <- mydata_baseline - rollmed_baseline
    
    mix_baseline <-
      mixtools::normalmixEM(myres_baseline[index_baseline[, 1], 1],
                            maxit = 10 ^ 8,
                            maxrestarts = 10 ^ 4)
    
    promo <- matrix(nrow = nrow(sales_data), ncol = ncol(sales_data))
    
    ## returning a matrix that will contain 1 when the product is on promotion and 0 otherwise
    for (j in 1:length(mix_baseline$posterior[, 1])) {
      if (mix_baseline$mu[1] > mix_baseline$mu[2]) {
        if (mix_baseline$posterior[j, 1] < 0.5) {
          promo[index_baseline[, 1], 1][j] <- 0
        } else if ((mix_baseline$posterior[j, 1] > 0.5) && mix_baseline$x[j] > 0) {
          promo[index_baseline[, 1], 1][j] <- 1
          mydata_baseline[index_baseline[, 1], 1][j] <-
            rollmed_baseline[index_baseline[, 1], 1][j]
        } else if ((mix_baseline$posterior[j, 1] > 0.5) && mix_baseline$x[j] <= 0) {
          promo[index_baseline[, 1], 1][j] <- 0
        } else {
          NA
        }
      } else {
        if (mix_baseline$posterior[j, 1] > 0.5) {
          promo[index_baseline[, 1], 1][j] <- 0
        } else if ((mix_baseline$posterior[j, 1] < 0.5) && mix_baseline$x[j] > 0) {
          promo[index_baseline[, 1], 1][j] <- 1
          mydata_baseline[index_baseline[, 1], 1][j] <-
            rollmed_baseline[index_baseline[, 1], 1][j]
        } else if ((mix_baseline$posterior[j, 1] < 0.5) && mix_baseline$x[j] <= 0) {
          promo[index_baseline[, 1], 1][j] <- 0
        } else {
          NA
        }
      }
    }
  } else if (promo_done == FALSE) {
    mydata_baseline <- sales_data
    promo <- rep(0, nrow(sales_data))
  }
  
  
  if (smoother == "mean") {
    smoothed_baseline <- zoo::rollapply(mydata_baseline,
                                        sizeroll ,
                                        FUN = mean,
                                        fill = NA,
                                        partial = TRUE,
                                        coredata = FALSE) %>%
      as.data.frame
  } else if (smoother == "median") {
    smoothed_baseline <- zoo::rollapply(mydata_baseline,
                                        sizeroll ,
                                        FUN = median,
                                        fill = NA,
                                        partial = TRUE,
                                        coredata = FALSE) %>%
      as.data.frame
  } else if (smoother == "loess") {
    smoothed_baseline <- stats::predict(
      loess(mydata_baseline[, 1] ~ c(1:length(mydata_baseline[, 1])), span = 0.4)
    ) %>% as.data.frame
  }
  
  a <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(week, smoothed_baseline[, 1]), col = 'red') +
    ggplot2::geom_line(ggplot2::aes(week, sales_data[, 1])) +
    ggplot2::ggtitle(stringr::str_c("The final baseline")) +
    ggplot2::labs(x = "Weeks", y = "Sales Quantity") +
    my_theme()
  
  
  
  if (showgraph == TRUE) {
    return(list(x = sales_data,
                baseline = smoothed_baseline[, 1],
                promotions = promo,
                plot = a))
  } else {
    return(list(x = sales_data, baseline = smoothed_baseline[, 1], promotions = promo))
  }
}


#' @title Creates the forecast of the baseline of sales based on historical data
#' @description This function enables to create the baseline sales forecast based on historical
#' sales data. It generates an historical baseline, considers a multitude of models
#' and selects the model that had the best performance on the testing set.
#' @param sales_data A vector containing historical sales data
#' @param frequency A numerical value specifying the frequency of the seasonality
#' @param start A vector of length 2 with the date of the first observation.
#' It contains first the year and then the day/week/month according to your data.
#' @param end A vector of length 2 with the date of the last observation.
#' It contains first the year and then the day/week/month according to your data.
#' @param forecast_horizon An integer value specifying the number of observations to forecast.
#' @param promo_done A logical variable specifying if promotions are done for the product.
#' @return A list containing the select model, the associated graphs, the predictions and the
#' confidence intervals, the accuracy measures and the same elements for all other considered models.
#' @author Grandadam Patrik
#' @importFrom magrittr %>%
#' @export
#' @examples
#' predict_baseline ###################### A COMPLETER ###############

predict_baseline <- function(sales_data,
                             frequency = 52,
                             start = c(2014, 1),
                             end = c(2019, 12),
                             forecast_horizon = 52,
                             promo_done = FALSE
) {
  
  # the original data, train and test set converted into a time series
  ts_actuals <- stats::ts(
    sales_data, start = start, end = end, frequency = forecast_horizon)
  ts_actuals_train <- base::subset(
    ts_actuals, end = length(ts_actuals) - forecast_horizon) # train set
  ts_actuals_test <- base::subset(
    ts_actuals, start = length(ts_actuals) - forecast_horizon + 1) # test set
  
  
  # Creation of training and testing set (not time series format)
  size.tr.set <- length(sales_data) - forecast_horizon
  size.te.set <- forecast_horizon
  tr.set <- sales_data[1:size.tr.set]
  te.set <- sales_data[(size.tr.set+1):length(sales_data)]
  
  # Creation of the baseline and training/testing set of the baseline
  original_baseline <- baseline(sales_data, promo_done = promo_done)$baseline
  ts_baseline <- stats::ts(
    original_baseline, start = start, end = end, frequency = forecast_horizon)
  ts_baseline_train <- base::subset(
    ts_baseline, end = length(ts_actuals) - forecast_horizon) # train set
  ts_baseline_test <- base::subset(
    ts_baseline, start = length(ts_actuals) - forecast_horizon + 1) # test set
  
  
  ## auto.arima on the baseline
  arima_baseline_train <- ts_baseline_train %>% forecast::auto.arima() # arima on train set
  fcst_arima_baseline <- arima_baseline_train %>%
    forecast::forecast(h = forecast_horizon, level = 0.95) # forecasts baseline test
  fcst_arima_future <- ts_baseline %>% forecast::auto.arima() %>%
    forecast::forecast(h = forecast_horizon, level = 0.95) # arima on everything
  
  plot_arima_baseline <- ggplot2::autoplot(ts_baseline_test) +
    ggplot2::autolayer(fcst_arima_future,
                       series = "Forecasted future baseline", level = FALSE) +
    ggplot2::autolayer(ts_baseline, series = "Baseline") +
    ggplot2::autolayer(fcst_arima_baseline,
                       series = "Forecasted baseline on testing set", level = FALSE) +
    ggplot2::autolayer(ts_actuals, series = "Actuals") +
    ggplot2::geom_ribbon(data = ts_baseline_test,
                         ggplot2::aes(ymin = fcst_arima_baseline$lower, ymax = fcst_arima_baseline$upper),
                         fill = "blue", alpha = "0.3") +
    ggplot2::geom_ribbon(data = fcst_arima_future$mean,
                         ggplot2::aes(ymin = fcst_arima_future$lower, ymax = fcst_arima_future$upper),
                         fill = "red", alpha = "0.3") +
    ggplot2::scale_color_manual(values = c("black", "grey", "blue", "red")) +
    ggplot2::ggtitle("Actuals, baseline and forecasted baseline using ARIMA") +
    ggplot2::xlab("Years") +
    ggplot2::ylab("Sales") +
    my_theme()
  acc_arima_baseline <- fcst_arima_baseline %>% forecast::accuracy(ts_baseline) # accuracy
  mape_arima_baseline <- acc_arima_baseline[2, "MAPE"] # MAPE
  arima_model_baseline <- list(NAME = "Arima model",
                               TEST_SET = fcst_arima_baseline$mean,
                               FORECAST = fcst_arima_future$mean,
                               lower_95 = fcst_arima_future$lower,
                               Upper_95 = fcst_arima_future$upper,
                               PLOT = plot_arima_baseline,
                               ACCURACIES = acc_arima_baseline,
                               MAPE = mape_arima_baseline)
  
  ## stlf on the baseline
  stlf_baseline_train <- ts_baseline_train %>% forecast::stlf(level = 0.95) # stlf on train set
  fcst_stlf_baseline <- stlf_baseline_train %>%
    forecast::forecast(h = forecast_horizon) # forecasts baseline test
  fcst_stlf_baseline_future <- ts_baseline %>% forecast::stlf(level = 0.95) %>%
    forecast::forecast(h = forecast_horizon) # stlf on everything
  
  plot_stlf_baseline <- ggplot2::autoplot(ts_baseline_test) +
    ggplot2::autolayer(fcst_stlf_baseline_future,
                       series = "Forecasted future baseline", level = FALSE) +
    ggplot2::autolayer(ts_baseline, series = "Baseline") +
    ggplot2::autolayer(fcst_stlf_baseline,
                       series = "Forecasted baseline on testing set", level = FALSE) +
    ggplot2::autolayer(ts_actuals, series = "Actuals") +
    ggplot2::geom_ribbon(data = ts_baseline_test,
                         ggplot2::aes(ymin = fcst_stlf_baseline$lower,
                                      ymax = fcst_stlf_baseline$upper),
                         fill = "blue", alpha = "0.3") +
    ggplot2::geom_ribbon(data = fcst_stlf_baseline_future$mean,
                         ggplot2::aes(ymin = fcst_stlf_baseline_future$lower,
                                      ymax = fcst_stlf_baseline_future$upper),
                         fill = "red", alpha = "0.3") +
    ggplot2::scale_color_manual(values = c("black", "grey", "blue", "red")) +
    ggplot2::ggtitle("Actuals, baseline and forecasted baseline using stlf") +
    ggplot2::xlab("Years") +
    ggplot2::ylab("Sales") +
    my_theme()
  acc_stlf_baseline <- fcst_stlf_baseline %>% forecast::accuracy(ts_baseline) # accuracy
  mape_stlf_baseline <- acc_stlf_baseline[2, "MAPE"] # MAPE
  stlf_model_baseline <- list(NAME = "stlf model",
                              TEST_SET = fcst_stlf_baseline$mean,
                              FORECAST = fcst_stlf_baseline_future$mean,
                              lower_95 = fcst_stlf_baseline_future$lower,
                              Upper_95 = fcst_stlf_baseline_future$upper,
                              PLOT = plot_stlf_baseline,
                              ACCURACIES = acc_stlf_baseline,
                              MAPE = mape_stlf_baseline)
  
  
  ## nnet on the baseline
  nnet_baseline_train <- ts_baseline_train %>% forecast::nnetar() # nnet on train set
  fcst_nnet_baseline <- nnet_baseline_train %>%
    forecast::forecast(h = forecast_horizon) # forecasts baseline test
  fcst_nnet_baseline_future <- ts_baseline %>% forecast::nnetar() %>%
    forecast::forecast(h = forecast_horizon) # nnet on everything
  
  plot_nnet_baseline <- ggplot2::autoplot(ts_baseline_test) +
    ggplot2::autolayer(fcst_nnet_baseline_future,
                       series = "Forecasted future baseline", level = FALSE) +
    ggplot2::autolayer(ts_baseline, series = "Baseline") +
    ggplot2::autolayer(fcst_nnet_baseline,
                       series = "Forecasted baseline on testing set", level = FALSE) +
    ggplot2::autolayer(ts_actuals, series = "Actuals") +
    ggplot2::scale_color_manual(values = c("black", "grey", "blue", "red")) +
    ggplot2::ggtitle("Actuals, baseline and forecasted baseline using Neural Network") +
    ggplot2::xlab("Years") +
    ggplot2::ylab("Sales") +
    ggplot2::scale_x_continuous(breaks = seq(2014, 2020, 1)) +
    my_theme()
  acc_nnet_baseline <- fcst_nnet_baseline %>% forecast::accuracy(ts_baseline) # accuracy
  mape_nnet_baseline <- acc_nnet_baseline[2, "MAPE"] # MAPE
  nnet_model_baseline <- list(NAME = "nnet model",
                              TEST_SET = fcst_nnet_baseline$mean,
                              FORECAST = fcst_nnet_baseline_future$mean,
                              lower_95 = fcst_nnet_baseline_future$lower,
                              Upper_95 = fcst_nnet_baseline_future$upper,
                              PLOT = plot_nnet_baseline,
                              ACCURACIES = acc_nnet_baseline,
                              MAPE = mape_nnet_baseline)
  
  
  ### bootstrapping on historical data
  n.boot <- 200L # number of simulations
  boot_train_baseline <- ts_baseline_train %>%
    forecast::bld.mbb.bootstrap(num = n.boot)
  boot_baseline <- ts_baseline %>%
    forecast::bld.mbb.bootstrap(num = n.boot)
  
  ### using stlf on the bootstrapped series
  fcst_boot_stlf <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
  fcst_boot_stlf_future <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
  for(i in seq(n.boot)) {
    fcst_boot_stlf[, i] <-
      forecast::stlf(boot_train_baseline[[i]], h = forecast_horizon, level = 95)$mean
    fcst_boot_stlf_future[, i] <-
      forecast::stlf(boot_baseline[[i]], h = forecast_horizon, level = 95)$mean
  }
  
  fcst_boot_stlf_baseline <- list(
    mean = stats::ts(rowMeans(fcst_boot_stlf),
                     start = stats::time(ts_baseline_test)[[1]], frequency = frequency),
    lower = stats::ts(apply(fcst_boot_stlf, 1, quantile, prob = 0.025),
                      start = stats::time(ts_baseline_test)[[1]], frequency = frequency),
    upper = stats::ts(apply(fcst_boot_stlf, 1, quantile, prob = 0.975),
                      start = stats::time(ts_baseline_test)[[1]], frequency = frequency)
  )
  
  
  fcst_boot_stlf_future <- list(
    mean = stats::ts(rowMeans(fcst_boot_stlf_future),
                     start = stats::time(ts_baseline_test)[[1]] + 1, frequency = frequency),
    lower = stats::ts(apply(fcst_boot_stlf_future, 1, quantile, prob = 0.025),
                      start = stats::time(ts_baseline_test)[[1]] + 1, frequency = frequency),
    upper = stats::ts(apply(fcst_boot_stlf_future, 1, quantile, prob = 0.975),
                      start = stats::time(ts_baseline_test)[[1]] + 1, frequency = frequency)
  )
  
  
  
  
  plot_boot_stlf_baseline <- ggplot2::autoplot(ts_baseline_test) +
    ggplot2::autolayer(fcst_boot_stlf_future$mean,
                       series = "Forecast of future sales") +
    ggplot2::autolayer(ts_baseline, series = "baseline") +
    ggplot2::autolayer(fcst_boot_stlf_baseline$mean,
                       series = "Forecast on testing set") +
    ggplot2::geom_ribbon(data = ts_baseline_test,
                         ggplot2::aes(ymin = fcst_boot_stlf_baseline$lower, ymax = fcst_boot_stlf_baseline$upper),
                         fill = "blue", alpha = "0.3") +
    ggplot2::geom_ribbon(data = fcst_boot_stlf_future$mean,
                         ggplot2::aes(ymin = fcst_boot_stlf_future$lower, ymax = fcst_boot_stlf_future$upper),
                         fill = "red", alpha = "0.3") +
    ggplot2::autolayer(ts_actuals, series = "Actuals") +
    ggplot2::scale_color_manual(values = c("black", "grey", "red", "blue")) +
    ggplot2::ggtitle("Actuals, baseline and forecasted baseline using bootstrap and STLF") +
    ggplot2::xlab("Years") +
    ggplot2::ylab("Sales") +
    ggplot2::scale_x_continuous(breaks = seq(2014, 2020, 1)) +
    my_theme()
  
  acc_boot_stlf_baseline <- forecast::accuracy(fcst_boot_stlf_baseline$mean, ts_baseline)
  mape_boot_stlf_baseline <- acc_boot_stlf_baseline[, "MAPE"] # MAPE
  boot_stlf_model_baseline <- list(NAME = "boot_stlf model",
                                   TEST_SET = fcst_boot_stlf_baseline$mean,
                                   FORECAST = fcst_boot_stlf_future$mean,
                                   lower_95 = fcst_boot_stlf_future$lower,
                                   Upper_95 = fcst_boot_stlf_future$upper,
                                   PLOT = plot_boot_stlf_baseline,
                                   ACCURACIES = acc_boot_stlf_baseline,
                                   MAPE = mape_boot_stlf_baseline)
  
  
  
  ### using nnet on the bootstrapped series
  fcst_boot_nnet <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
  fcst_boot_nnet_future <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
  model_nnet <- list()
  model_nnet_future <- list()
  for(i in seq(n.boot)) {
    model_nnet[[i]] <- forecast::nnetar(boot_train_baseline[[i]])
    fitted_values <- model_nnet[[i]] %>% forecast::forecast(h = forecast_horizon)
    fitted_values <- fitted_values$mean
    fcst_boot_nnet[, i] <- fitted_values # forecasts on test set
    
    model_nnet_future[[i]] <- forecast::nnetar(boot_baseline[[i]])
    fitted_values_future <- model_nnet_future[[i]] %>% forecast::forecast(h = forecast_horizon)
    fitted_values_future <- fitted_values_future$mean
    fcst_boot_nnet_future[, i] <- fitted_values_future # forecasts of future values of baseline
  }
  
  fcst_boot_nnet_baseline <- list(
    mean = stats::ts(rowMeans(fcst_boot_nnet),
                     start = stats::time(ts_baseline_test)[[1]], frequency = frequency),
    lower = stats::ts(apply(fcst_boot_nnet, 1, quantile, prob = 0.025),
                      start = stats::time(ts_baseline_test)[[1]], frequency = frequency),
    upper = stats::ts(apply(fcst_boot_nnet, 1, quantile, prob = 0.975),
                      start = stats::time(ts_baseline_test)[[1]], frequency = frequency)
  )
  
  
  fcst_boot_nnet_future <- list(
    mean = stats::ts(rowMeans(fcst_boot_nnet_future),
                     start = stats::time(ts_baseline_test)[[1]] + 1, frequency = frequency),
    lower = stats::ts(apply(fcst_boot_nnet_future, 1, quantile, prob = 0.025),
                      start = stats::time(ts_baseline_test)[[1]] + 1, frequency = frequency),
    upper = stats::ts(apply(fcst_boot_nnet_future, 1, quantile, prob = 0.975),
                      start = stats::time(ts_baseline_test)[[1]] + 1, frequency = frequency)
  )
  
  
  
  plot_boot_nnet_baseline <- ggplot2::autoplot(ts_baseline_test) +
    ggplot2::autolayer(fcst_boot_nnet_future$mean,
                       series = "Forecast of future sales") +
    ggplot2::autolayer(ts_baseline, series = "baseline") +
    ggplot2::autolayer(fcst_boot_nnet_baseline$mean,
                       series = "Forecast on testing set") +
    ggplot2::geom_ribbon(data = ts_baseline_test,
                         ggplot2::aes(ymin = fcst_boot_nnet_baseline$lower, ymax = fcst_boot_nnet_baseline$upper),
                         fill = "blue", alpha = "0.3") +
    ggplot2::geom_ribbon(data = fcst_boot_nnet_future$mean,
                         ggplot2::aes(ymin = fcst_boot_nnet_future$lower, ymax = fcst_boot_nnet_future$upper),
                         fill = "red", alpha = "0.3") +
    ggplot2::autolayer(ts_actuals, series = "Actuals") +
    ggplot2::scale_color_manual(values = c("black", "grey", "red", "blue")) +
    ggplot2::ggtitle("Actuals, baseline and forecasted baseline using bootstrap and nnet") +
    ggplot2::xlab("Years") +
    ggplot2::ylab("Sales") +
    ggplot2::scale_x_continuous(breaks = seq(2014, 2020, 1)) +
    my_theme()
  
  acc_boot_nnet_baseline <- forecast::accuracy(fcst_boot_nnet_baseline$mean, ts_baseline)
  mape_boot_nnet_baseline <- acc_boot_nnet_baseline[, "MAPE"] # MAPE
  
  boot_nnet_model_baseline <- list(NAME = "boot_nnet model",
                                   TEST_SET = fcst_boot_nnet_baseline$mean,
                                   FORECAST = fcst_boot_nnet_future$mean,
                                   lower_95 = fcst_boot_nnet_future$lower,
                                   Upper_95 = fcst_boot_nnet_future$upper,
                                   PLOT = plot_boot_nnet_baseline,
                                   ACCURACIES = acc_boot_nnet_baseline,
                                   MAPE = mape_boot_nnet_baseline)
  
  
  
  mapes <- c(mape_arima_baseline = mape_arima_baseline,
             mape_boot_nnet_baseline = mape_boot_nnet_baseline,
             mape_boot_stlf_baseline = mape_boot_stlf_baseline,
             mape_nnet_baseline = mape_nnet_baseline,
             mape_stlf_baseline = mape_stlf_baseline)
  
  if(names(which.min(mapes)) ==  "mape_arima_baseline") {
    retained_model <- arima_model_baseline
  } else if(names(which.min(mapes)) ==  "mape_boot_nnet_baseline") {
    retained_model <- boot_nnet_model_baseline
  } else if(names(which.min(mapes)) ==  "mape_boot_stlf_baseline") {
    retained_model <- boot_stlf_model_baseline
  } else if(names(which.min(mapes)) ==  "mape_nnet_baseline") {
    retained_model <- nnet_model_baseline
  } else if(names(which.min(mapes)) ==  "mape_stlf_baseline") {
    retained_model <- stlf_model_baseline
  }
  
  
  return(list(selected_model = retained_model,
              all_models = list(arima_baseline = arima_model_baseline,
                                bootstrapped_nnet_baseline = boot_nnet_model_baseline,
                                bootstrapped_stlf_baseline = boot_stlf_model_baseline,
                                neural_network = nnet_model_baseline,
                                stlf_baseline =   stlf_model_baseline)))
  
  
  
  
  
}


#' @title Creates sales forecast based on historical data
#' @description This function enables to create sales forecast based on historical
#' sales data. It considers a multitude of models and selects the model that had
#' the best performance on the testing set.
#' @param sales_data A vector containing historical sales data
#' @param frequency A numerical value specifying the frequency of the seasonality
#' @param start A vector of length 2 with the date of the first observation.
#' It contains first the year and then the day/week/month according to your data.
#' @param end A vector of length 2 with the date of the last observation.
#' It contains first the year and then the day/week/month according to your data.
#' @param forecast_horizon An integer value specifying the number of observations to forecast.
#' @param promo_done A logical variable specifying if promotions are done for the product.
#' @param future_promotions Optionnal: a vector composed on binary variables which are
#' equal to 1 when there will be a promotion in the forecasting horizon and to 0 otherwise.
#' @return A list containing the select model, the associated graphs, the predictions and the
#' confidence intervales, the accuracy measures and the same elements for all other considered models.
#' @author Grandadam Patrik
#' @importFrom magrittr %>%
#' @export
#' @examples
#' predict_sales ###################### A COMPLETER ###############


predict_sales <- function(sales_data,
                          frequency = 52,
                          start = c(2014, 1),
                          end = c(2019, 12),
                          forecast_horizon = 52,
                          promo_done = FALSE,
                          future_promotions = NA) {
  
  if (length(future_promotions) != forecast_horizon && !is.na(future_promotions)) {
    warning("future_promotions should be the same length as the forecasting horizon")
  } else {
    
    # changing NA to 0
    sales_data[is.na(sales_data)] <- 0
    
    # the original data, train and test set converted into a time series
    ts_actuals <- stats::ts(
      sales_data, start = start, end = end, frequency = forecast_horizon)
    ts_actuals_train <- base::subset(
      ts_actuals, end = length(ts_actuals) - forecast_horizon) # train set
    ts_actuals_test <- base::subset(
      ts_actuals, start = length(ts_actuals) - forecast_horizon + 1) # test set
    
    
    # Creation of training and testing set (not time series format)
    size.tr.set <- length(sales_data) - forecast_horizon
    size.te.set <- forecast_horizon
    tr.set <- sales_data[1:size.tr.set]
    te.set <- sales_data[(size.tr.set + 1):length(sales_data)]
    
    
    
    
    
    ## auto.arima on historic data
    arima_actuals_train <- ts_actuals_train %>% forecast::auto.arima() # arima on train set
    fcst_arima_actuals <- arima_actuals_train %>%
      forecast::forecast(h = forecast_horizon, level = 0.95) # forecasts actuals test
    arima_future_model <- ts_actuals %>% forecast::auto.arima()
    fcst_arima_future <- arima_future_model %>%
      forecast::forecast(h = forecast_horizon, level = 0.95) # arima on everything
    
    plot_arima_actuals <- ggplot2::autoplot(ts_actuals_test) +
      ggplot2::autolayer(fcst_arima_future, series = "Forecast of future sales", level = FALSE) +
      ggplot2::autolayer(ts_actuals, series = "Actuals") +
      ggplot2::autolayer(fcst_arima_actuals, series = "Forecast on testing set", level = FALSE) +
      ggplot2::geom_ribbon(data = ts_actuals_test,
                           ggplot2::aes(ymin = fcst_arima_actuals$lower, ymax = fcst_arima_actuals$upper),
                           fill = "blue", alpha = "0.3") +
      ggplot2::geom_ribbon(data = fcst_arima_future$mean,
                           ggplot2::aes(ymin = fcst_arima_future$lower, ymax = fcst_arima_future$upper),
                           fill = "red", alpha = "0.3") +
      ggplot2::scale_color_manual(values = c("black", "red", "blue")) +
      ggplot2::ggtitle("Actuals and forecasted total sales using ARIMA") +
      ggplot2::xlab("Years") +
      ggplot2::ylab("Sales") +
      my_theme()
    acc_arima_actuals <- fcst_arima_actuals %>% forecast::accuracy(ts_actuals) # accuracy
    mape_arima_actuals <- acc_arima_actuals[2, "MAPE"] # MAPE
    arima_model <- list(NAME = "Arima model",
                        MODEL = arima_future_model,
                        TEST_SET = fcst_arima_actuals$mean,
                        FORECAST = fcst_arima_future$mean,
                        lower_95 = fcst_arima_future$lower,
                        Upper_95 = fcst_arima_future$upper,
                        PLOT = plot_arima_actuals,
                        ACCURACIES = acc_arima_actuals,
                        MAPE = mape_arima_actuals)
    
    ## stlf on historic data
    stlf_actuals_train <- ts_actuals_train %>% forecast::stlf(level = 0.95) # stlf on train set
    fcst_stlf_actuals <- stlf_actuals_train %>%
      forecast::forecast(h = forecast_horizon) # forecasts actuals test
    stlf_future_model <- ts_actuals %>% forecast::stlf(level = 0.95)
    fcst_stlf_future <- stlf_future_model %>%
      forecast::forecast(h = forecast_horizon) # stlf on everything
    
    plot_stlf_actuals <- ggplot2::autoplot(ts_actuals_test) +
      ggplot2::autolayer(fcst_stlf_future, series = "Forecast of future sales", level = FALSE) +
      ggplot2::autolayer(ts_actuals, series = "Actuals") +
      ggplot2::autolayer(fcst_stlf_actuals, series = "Forecast on testing set", level = FALSE) +
      ggplot2::geom_ribbon(data = ts_actuals_test,
                           ggplot2::aes(ymin = fcst_stlf_actuals$lower, ymax = fcst_stlf_actuals$upper),
                           fill = "blue", alpha = "0.3") +
      ggplot2::geom_ribbon(data = fcst_stlf_future$mean,
                           ggplot2::aes(ymin = fcst_stlf_future$lower, ymax = fcst_stlf_future$upper),
                           fill = "red", alpha = "0.3") +
      ggplot2::scale_color_manual(values = c("black", "red", "blue")) +
      ggplot2::ggtitle("Actuals and forecasted total sales using STLF") +
      ggplot2::xlab("Years") +
      ggplot2::ylab("Sales") +
      my_theme()
    acc_stlf_actuals <- fcst_stlf_actuals %>% forecast::accuracy(ts_actuals) # accuracy
    mape_stlf_actuals <- acc_stlf_actuals[2, "MAPE"] # MAPE
    stlf_model <- list(NAME = "stlf model",
                       MODEL = stlf_future_model,
                       TEST_SET = fcst_stlf_actuals$mean,
                       FORECAST = fcst_stlf_future$mean,
                       lower_95 = fcst_stlf_future$lower,
                       Upper_95 = fcst_stlf_future$upper,
                       PLOT = plot_stlf_actuals,
                       ACCURACIES = acc_stlf_actuals,
                       MAPE = mape_stlf_actuals)
    
    
    ## nnet on historic data
    nnet_actuals_train <- ts_actuals_train %>% forecast::nnetar(level = 0.95) # nnet on train set
    fcst_nnet_actuals <- nnet_actuals_train %>%
      forecast::forecast(h = forecast_horizon) # forecasts actuals test
    nnet_future_model <- ts_actuals %>% forecast::nnetar(level = 0.95)
    fcst_nnet_future <- nnet_future_model %>%
      forecast::forecast(h = forecast_horizon) # nnet on everything
    
    plot_nnet_actuals <- ggplot2::autoplot(ts_actuals_test) +
      ggplot2::autolayer(fcst_nnet_future, series = "Forecast of future sales", level = FALSE) +
      ggplot2::autolayer(ts_actuals, series = "Actuals") +
      ggplot2::autolayer(fcst_nnet_actuals, series = "Forecast on testing set", level = FALSE) +
      ggplot2::scale_color_manual(values = c("black", "red", "blue")) +
      ggplot2::ggtitle("Actuals and forecasted total sales using Neural Network") +
      ggplot2::xlab("Years") +
      ggplot2::ylab("Sales") +
      
      my_theme()
    acc_nnet_actuals <- fcst_nnet_actuals %>% forecast::accuracy(ts_actuals) # accuracy
    mape_nnet_actuals <- acc_nnet_actuals[2, "MAPE"] # MAPE
    nnet_model <- list(NAME = "nnet model",
                       TEST_SET = fcst_nnet_actuals$mean,
                       MODEL = nnet_future_model,
                       FORECAST = fcst_nnet_future$mean,
                       lower_95 = fcst_nnet_future$lower,
                       Upper_95 = fcst_nnet_future$upper,
                       PLOT = plot_nnet_actuals,
                       ACCURACIES = acc_nnet_actuals,
                       MAPE = mape_nnet_actuals)
    
    
    
    
    
    
    
    
    ### bootstrapping on historical data
    n.boot <- 200L # number of simulations
    boot_train <- ts_actuals_train %>%
      forecast::bld.mbb.bootstrap(num = n.boot)
    boot_actuals <- ts_actuals %>%
      forecast::bld.mbb.bootstrap(num = n.boot)
    
    ### using stlf on the bootstrapped series
    fcst_boot_stlf <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
    fcst_boot_stlf_future <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
    for(i in seq(n.boot)) {
      fcst_boot_stlf[, i] <-
        forecast::stlf(boot_train[[i]], h = forecast_horizon, level = 95)$mean
      fcst_boot_stlf_future[, i] <-
        forecast::stlf(boot_actuals[[i]], h = forecast_horizon, level = 95)$mean
    }
    
    fcst_boot_stlf_actuals <- structure(list(
      mean = stats::ts(rowMeans(fcst_boot_stlf),
                       start = stats::time(ts_actuals_test)[[1]], frequency = frequency),
      lower = stats::ts(apply(fcst_boot_stlf, 1, quantile, prob = 0.025),
                        start = stats::time(ts_actuals_test)[[1]], frequency = frequency),
      upper = stats::ts(apply(fcst_boot_stlf, 1, quantile, prob = 0.975),
                        start = stats::time(ts_actuals_test)[[1]], frequency = frequency),
      class = "forecast"
    )
    )
    
    fcst_boot_stlf_future <- structure(list(
      mean = stats::ts(rowMeans(fcst_boot_stlf_future),
                       start = stats::time(ts_actuals_test)[[1]] + 1, frequency = frequency),
      lower = stats::ts(apply(fcst_boot_stlf_future, 1, quantile, prob = 0.025),
                        start = stats::time(ts_actuals_test)[[1]] + 1, frequency = frequency),
      upper = stats::ts(apply(fcst_boot_stlf_future, 1, quantile, prob = 0.975),
                        start = stats::time(ts_actuals_test)[[1]] + 1, frequency = frequency),
      class = "forecast"
    )
    )
    
    
    
    plot_boot_stlf_actuals <- ggplot2::autoplot(ts_actuals_test) +
      ggplot2::autolayer(fcst_boot_stlf_future$mean,
                         series = "Forecast of future sales") +
      ggplot2::autolayer(ts_actuals, series = "Actuals") +
      ggplot2::autolayer(fcst_boot_stlf_actuals$mean,
                         series = "Forecast on testing set") +
      ggplot2::geom_ribbon(data = ts_actuals_test,
                           ggplot2::aes(ymin = fcst_boot_stlf_actuals$lower, ymax = fcst_boot_stlf_actuals$upper),
                           fill = "blue", alpha = "0.3") +
      ggplot2::geom_ribbon(data = fcst_stlf_future$mean,
                           ggplot2::aes(ymin = fcst_boot_stlf_future$lower, ymax = fcst_boot_stlf_future$upper),
                           fill = "red", alpha = "0.3") +
      ggplot2::scale_color_manual(values = c("black", "red", "blue")) +
      ggplot2::ggtitle("Actuals and forecasted total sales using bootstrap and STLF") +
      ggplot2::xlab("Years") +
      ggplot2::ylab("Sales") +
      
      my_theme()
    
    acc_boot_stlf_actuals <- forecast::accuracy(fcst_boot_stlf_actuals$mean, ts_actuals)
    mape_boot_stlf_actuals <- acc_boot_stlf_actuals[, "MAPE"] # MAPE
    
    
    boot_stlf_model <- list(NAME = "boot_stlf model",
                            TEST_SET = fcst_boot_stlf_actuals$mean,
                            FORECAST = fcst_boot_stlf_future$mean,
                            lower_95 = fcst_boot_stlf_future$lower,
                            Upper_95 = fcst_boot_stlf_future$upper,
                            PLOT = plot_boot_stlf_actuals,
                            ACCURACIES = acc_boot_stlf_actuals,
                            MAPE = mape_boot_stlf_actuals)
    
    ### using nnet on the bootstrapped series
    fcst_boot_nnet <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
    fcst_boot_nnet_future <- matrix(0, ncol = n.boot, nrow = forecast_horizon)
    model_nnet <- list()
    model_nnet_future <- list()
    for(i in seq(n.boot)) {
      model_nnet[[i]] <- forecast::nnetar(boot_train[[i]])
      fitted_values <- model_nnet[[i]] %>% forecast::forecast(h = forecast_horizon)
      fitted_values <- fitted_values$mean
      fcst_boot_nnet[, i] <- fitted_values # forecasts on test set
      
      model_nnet_future[[i]] <- forecast::nnetar(boot_actuals[[i]])
      fitted_values_future <- model_nnet_future[[i]] %>% forecast::forecast(h = forecast_horizon)
      fitted_values_future <- fitted_values_future$mean
      fcst_boot_nnet_future[, i] <- fitted_values_future # forecasts of future values
    }
    
    
    fcst_boot_nnet_actuals <- structure(list(
      mean = stats::ts(rowMeans(fcst_boot_nnet),
                       start = stats::time(ts_actuals_test)[[1]], frequency = frequency),
      lower = stats::ts(apply(fcst_boot_nnet, 1, quantile, prob = 0.025),
                        start = stats::time(ts_actuals_test)[[1]], frequency = frequency),
      upper = stats::ts(apply(fcst_boot_nnet, 1, quantile, prob = 0.975),
                        start = stats::time(ts_actuals_test)[[1]], frequency = frequency),
      class = "forecast"
    )
    )
    
    fcst_boot_nnet_future <- structure(list(
      mean = stats::ts(rowMeans(fcst_boot_nnet_future),
                       start = stats::time(ts_actuals_test)[[1]] + 1, frequency = frequency),
      lower = stats::ts(apply(fcst_boot_nnet_future, 1, quantile, prob = 0.025),
                        start = stats::time(ts_actuals_test)[[1]] + 1, frequency = frequency),
      upper = stats::ts(apply(fcst_boot_nnet_future, 1, quantile, prob = 0.975),
                        start = stats::time(ts_actuals_test)[[1]] + 1, frequency = frequency),
      class = "forecast"
    )
    )
    
    plot_boot_nnet_actuals <- ggplot2::autoplot(ts_actuals_test) +
      ggplot2::autolayer(fcst_boot_nnet_future$mean,
                         series = "Forecast of future sales") +
      ggplot2::autolayer(ts_actuals, series = "Actuals") +
      ggplot2::autolayer(fcst_boot_nnet_actuals$mean,
                         series = "Forecast on testing set") +
      ggplot2::geom_ribbon(data = ts_actuals_test,
                           ggplot2::aes(ymin = fcst_boot_nnet_actuals$lower, ymax = fcst_boot_nnet_actuals$upper),
                           fill = "blue", alpha = "0.3") +
      ggplot2::geom_ribbon(data = fcst_nnet_future$mean,
                           ggplot2::aes(ymin = fcst_boot_nnet_future$lower, ymax = fcst_boot_nnet_future$upper),
                           fill = "red", alpha = "0.3") +
      ggplot2::scale_color_manual(values = c("black", "red", "blue")) +
      ggplot2::ggtitle("Actuals and forecasted total sales using bootstrap and Neural Network") +
      ggplot2::xlab("Years") +
      ggplot2::ylab("Sales") +
      
      my_theme()
    
    acc_boot_nnet_actuals <- forecast::accuracy(fcst_boot_nnet_actuals$mean, ts_actuals)
    mape_boot_nnet_actuals <- acc_boot_nnet_actuals[, "MAPE"] # MAPE
    
    
    boot_nnet_model <- list(NAME = "boot_nnet model",
                            TEST_SET = fcst_boot_nnet_actuals$mean,
                            FORECAST = fcst_boot_nnet_future$mean,
                            lower_95 = fcst_boot_nnet_future$lower,
                            Upper_95 = fcst_boot_nnet_future$upper,
                            PLOT = plot_boot_nnet_actuals,
                            ACCURACIES = acc_boot_nnet_actuals,
                            MAPE = mape_boot_nnet_actuals)
    
    if(promo_done == TRUE) {
      promos <- baseline(sales_data, promo_done = TRUE)$promotions %>% as.vector()
      promos_train <- promos[1:size.tr.set]
      promos_test <- promos[(size.tr.set + 1):length(sales_data)]
      
      dyna_actuals_train <- ts_actuals_train %>%
        forecast::auto.arima(xreg = promos_train) # dyna on train set
      fcst_dyna_actuals <- dyna_actuals_train %>%
        forecast::forecast(xreg = promos_test, h = forecast_horizon, level = 0.95) # forecasts
      fcst_dyna_future <- ts_actuals %>%
        forecast::auto.arima(xreg = promos) %>%
        forecast::forecast(h = forecast_horizon, xreg = future_promotions, level = 0.95) # forecasts
      
      
      plot_dyna_actuals <- ggplot2::autoplot(ts_actuals_test) +
        ggplot2::autolayer(fcst_dyna_future, series = "Forecast of future sales", level = FALSE) +
        ggplot2::autolayer(ts_actuals, series = "Actuals") +
        ggplot2::autolayer(fcst_dyna_actuals, series = "Forecast on testing set", level = FALSE) +
        ggplot2::geom_ribbon(data = ts_actuals_test,
                             ggplot2::aes(ymin = fcst_dyna_actuals$lower, ymax = fcst_dyna_actuals$upper),
                             fill = "blue", alpha = "0.3") +
        ggplot2::geom_ribbon(data = fcst_dyna_future$mean,
                             ggplot2::aes(ymin = fcst_dyna_future$lower, ymax = fcst_dyna_future$upper),
                             fill = "red", alpha = "0.3") +
        ggplot2::scale_color_manual(values = c("black", "red", "blue")) +
        ggplot2::ggtitle("Actuals and forecasted total sales using a dynamic model") +
        ggplot2::xlab("Years") +
        ggplot2::ylab("Sales") +
        my_theme()
      acc_dyna_actuals <- fcst_dyna_actuals %>% forecast::accuracy(ts_actuals) # accuracy
      mape_dyna_actuals <- acc_dyna_actuals[2, "MAPE"] # MAPE
      
    } else {
      dyna_actuals_train <- NULL
      fcst_dyna_actuals <- NULL
      plot_dyna_actuals <- NULL
      acc_dyna_actuals <- NULL
      mape_dyna_actuals <- 10^8
    }
    
    dyna_model <- list(NAME = "Dynamic model",
                       FORECAST = fcst_dyna_actuals$mean,
                       lower_95 = fcst_dyna_actuals$lower,
                       Upper_95 = fcst_dyna_actuals$upper,
                       PLOT = plot_dyna_actuals,
                       ACCURACIES = acc_dyna_actuals,
                       MAPE = mape_dyna_actuals)
    
    # performance of each model
    mapes <- c(mape_arima_actuals = mape_arima_actuals,
               mape_boot_nnet_actuals = mape_boot_nnet_actuals,
               mape_boot_stlf_actuals = mape_boot_stlf_actuals,
               mape_dyna_actuals = mape_dyna_actuals,
               mape_nnet_actuals = mape_nnet_actuals,
               mape_stlf_actuals = mape_stlf_actuals)
    
    if(names(which.min(mapes)) ==  "mape_arima_actuals") {
      retained_model <- arima_model
    } else if(names(which.min(mapes)) ==  "mape_boot_nnet_actuals") {
      retained_model <- boot_nnet_model
    } else if(names(which.min(mapes)) ==  "mape_boot_stlf_actuals") {
      retained_model <- boot_stlf_model
    } else if(names(which.min(mapes)) ==  "mape_dyna_actuals") {
      retained_model <- dyna_model
    } else if(names(which.min(mapes)) ==  "mape_nnet_actuals") {
      retained_model <- nnet_model
    } else if(names(which.min(mapes)) ==  "mape_stlf_actuals") {
      retained_model <- stlf_model
    }
    
    
    return(list(selected_model = retained_model,
                all_models = list(arima = arima_model,
                                  bootstrapped_nnet = boot_nnet_model,
                                  bootstrapped_stlf= boot_stlf_model,
                                  dynamic_model = dyna_model,
                                  neural_network = nnet_model,
                                  stlf_baseline =   stlf_model)))
    
    
  }
}











#' @title A great theme for graphs issued from the ggplot2 plackage
#' @description This function enables to have a nice looking theme for ggplot graphs.
#' @param base_size The base size, no need to change it.
#' @param base_family The base_family, no need to change it
#' @author Grandadam Patrik
#' @export
my_theme <- function(base_size = 10, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 0.5),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.major = ggplot2::element_line(color = "grey"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "aliceblue"),
      strip.background = ggplot2::element_rect(
        color = "grey", size = 1, fill = "lightgrey"),
      strip.text = ggplot2::element_text(color = "black", size = 10, face = "bold"),
      legend.position = "bottom",
      legend.justification = "top",
      legend.box = "horizontal",
      legend.box.background = ggplot2::element_rect(colour = "grey50"),
      legend.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "grey", fill = NA, size = 0.5)
    )
}






