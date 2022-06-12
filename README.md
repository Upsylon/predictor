## Context :books:

This package comes with my Master Thesis in *Business Analytics* orientation of the Master in Management at HEC Lausanne.  
I obtained the grade of 6/6.  
The slides of the presentation are available [here](https://github.com/Upsylon/predictor/blob/master/Thesis%20presentation.pdf).

## Objectives of the package :dart:

1. Automate the forecasting process of time series  
--> Based on historical sales data, forecasts are computed automatically.  
2. Ensure a reproducible methodology.   
--> The library provides a variety of easily usable functions
3. Present the results visually and interactively  
--> This package comes with a [ **Shiny App**](https://upsy.shinyapps.io/Predictor/)

The methodology includes:  
- An outlier detection algorithm.  
- Predictive models:  
     - Time-Series: ARIMA, STLF, ETS  
     - Machine Learning: Neural-Networks  

## Downloading the package :cd:

The package can be downloaded using the following functions:

```{r, eval = FALSE, echo = TRUE}
if (!require("devtools"))
  install.packages("devtools")
  
devtools::install_github("Upsylon/predictor") # install the package
library(predictor) # load the package
```

## Utilisation :computer:

The package is composed of different functions. Their objectives and utilisation is detailed in their description.

