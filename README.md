## Objectives of the package :dart:

The objective package is designed to **automate the forecasting process** of sales.   
It implements a complete procedure that leads to accurate results using a **reproducible methodology**.   
Finally, it makes the life of the user simpler by **minimizing the code** needed.     

Based on historical sales data, the forecasts are computed automatically.

It includes an outlier detection algorithm and the use of multiple predictive models (ARIMA, STLF, ETS, Neural-Networks). The performance of each model is then evaluated. The best model is selected and the predictions are computed.

This package comes with a **Shiny App**:  
https://upsy.shinyapps.io/Predictor_Baseline/

## Downloading the package :cd:

The package can be downloaded using the following functions:

```{r, eval = FALSE, echo = TRUE}
if (!require("devtools"))
  install.packages("devtools")
  
devtools::install_github("Upsylon/predictor") # install the package
library(predictor) # load the package
```

## Utilisation :computer:

The package is composed of six functions:  
    - baseline()  
    - predict_baseline()  
    - predict_sales()  
    - my_theme()  
    - my_theme()  
    - oos_simulation()
    
Their objectives and utilisation is detailed in their description.



