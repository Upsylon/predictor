## Context :books:

This package has been built as part of my Master Thesis in the *Business Analytics* orientation of the Master in Management at HEC Lausanne, for which I obtained the grade of 6/6.  
The slides of the presentation are available in this repository.

## Objectives of the package :dart:

The objective package is designed to **automate the forecasting process** of sales.   
Based on historical sales data, the forecasts are computed automatically.  

It implements a complete procedure that leads to accurate results using a **reproducible methodology**.   
Finally, it makes the life of the user simpler by **minimizing the code** needed.     

It includes an outlier detection algorithm and the use of multiple predictive models (ARIMA, STLF, ETS, Neural-Networks). The performance of each model is then evaluated. The best model is selected and the predictions are computed.

This package comes with a **Shiny App**:  
https://upsy.shinyapps.io/Predictor/

With its ease of use, the Shiny App would enable a baby to make accurate predictions! :grin:  
The user can **load data**, **vizualize** the shape of the predictions and **download** the results.   
The default parameters ensure a meaningful initial setting that the user is free to change. All this **without any code**.

The "Sample.csv" file in this repository can be used for demonstration purpose. 

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



