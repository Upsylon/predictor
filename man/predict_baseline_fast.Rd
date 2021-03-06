% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Predictor.R
\name{predict_baseline_fast}
\alias{predict_baseline_fast}
\title{Creates the forecast of the baseline of sales based on historical data using only fast models}
\usage{
predict_baseline_fast(sales_data, frequency = 52, start = c(2014, 1),
  forecast_horizon = 52, size.te.set = 52, promo_done = FALSE,
  criterion = "accuracy", sizeroll = 11)
}
\arguments{
\item{sales_data}{A vector containing historical sales data.}

\item{frequency}{A numerical value specifying the frequency of the seasonality.}

\item{start}{A vector of length 2 with the date of the first observation.
It contains first the year and then the day/week/month according to your data.}

\item{forecast_horizon}{An integer value specifying the number of observations to forecast.}

\item{size.te.set}{An integer value specifying the size of the testing set.}

\item{promo_done}{A logical variable specifying if promotions are done for the product.}

\item{criterion}{A string variable specifying the selection criterion that should be used to
select the model ("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1", "Theil's U"). "accuracy"
can also be used to reflect the needs of the company.}

\item{sizeroll}{The window of the moving average or moving median when using 
the baseline() function.}

\item{smoother}{The smoother that should be considered when using the baseline() function. 
It can be "mean", "median" or "loess".}
}
\value{
A list containing the select model, the associated graphs, the predictions and the
confidence intervals, the accuracy measures and the same elements for all other considered models.
}
\description{
This function has the same purpose as the predict_baseline() function but uses 
only the models which computation time are the smallest.
}
\examples{
data("mydata")
my_baseline <- predict_baseline_baseline(mydata, promo_done = TRUE, criterion = "MAPE")
my_baseline$selected_model$PLOT # the plot of the selected model
my_baseline$selected_model$FORECAST # the forecast of the selected model
my_baseline$selected_model$ACCURACIES # the accuracies of the selected model
}
\author{
Grandadam Patrik
}
