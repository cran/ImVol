#'@title Volume Prediction of Trees Using Linear and Nonlinear Allometric Equations
#' @param data Datasets
#' @import stats tidyverse nls2 caret ggplot2 dplyr tidyr
#' @return
#' \itemize{
#'   \item results: Results
#' }
#' @export
#'
#' @examples
#' library("ImVol")
#' data <- system.file("extdata", "data_test.csv", package = "ImVol")
#' Data<- read.csv(data)
#' Vol <- ImVoL(Data)
#' @references
#' \itemize{
#'\item Sharma, A., Gupta.R.K., Mahajan, P.K and Shilp. 2017.Volume Equation for Cedrus deodara in Kullu District of Himachal Pradesh. Indian Journal of Ecology . 44(6) : 781-784.http://dx.doi.org/10.13140/RG.2.2.33786.62407
#'\item Tabassum, A., Jeelani, M.I.,Sharma,M., Rather, K R ., Rashid, I and Gul,M.2022.  Predictive Modelling of Height and Diameter Relationships of Himalayan Chir Pine .Agricultural Science Digest - A Research Journal. DOI:10.18805/ag.D-555
#'}
ImVoL <- function(data) {
  Metric<-NULL
  Value<-NULL
  Model<-NULL
  # Initialize a list to store results for different splits and cross-validation methods
  all_results <- list()

  # Define the list of splits
  splits <- list(
    "50:50" = c(0.5, 0.5),
    "70:30" = c(0.7, 0.3),
    "80:20" = c(0.8, 0.2)
  )

  # Perform modeling for each split and cross-validation method
  for (split_name in names(splits)) {
    split <- splits[[split_name]]

    # Split the data
    sample_indices <- sample(1:nrow(data), split[1] * nrow(data))
    train_data <- data[sample_indices, ]
    test_data <- data[-sample_indices, ]

    # Define the list of models
    models<-NULL
    models <- list(
      Linear.D1 = lm(Volume ~ Diameter, data = train_data),
      Logarithmic.D2 = lm(Volume ~ log(Diameter), data = train_data),
      Inverse.D3 = lm(Volume ~ 1 / Diameter, data = train_data),
      Quadratic.D4 = lm(Volume ~ Diameter + I(Diameter^2), data = train_data),
      Cubic.D5 = lm(Volume ~ Diameter + I(Diameter^2) + I(Diameter^3), data = train_data),
      Compound.D6 = nls2(Volume ~ a * b^Diameter, start = list(a = 1, b = 1), algorithm = "brute-force", data = train_data),
      Power.D7 = nls2(Volume ~ a * Diameter^b, start = list(a = 1, b = 0.01), algorithm = "brute-force", data = train_data),
      Exponential.D8 = nls2(Volume ~ a * exp(Diameter * b), start = list(a = 1, b = 0.01), algorithm = "brute-force", data = train_data),
      Linear.HD1 = lm(Volume ~ Height + Diameter, data = train_data),
      Logarithmic.HD2 = lm(Volume ~ Height + log(Diameter), data = train_data),
      Inverse.HD3 = lm(Volume ~ Height + I(1 / Diameter), data = train_data),
      Quadratic.HD4 = lm(Volume ~ Height + Diameter + I(Diameter^2), data = train_data),
      Cubic.HD5 = lm(Volume ~ Height + Diameter + I(Diameter^2) + I(Diameter^3), data = train_data),
      Compound.HD6 = nls2(Volume ~ a * b^Diameter, start = list(a = 1, b = 1), algorithm = "brute-force", data = train_data),
      Power.HD7 = nls2(Volume ~ a * Diameter^b, start = list(a = 1, b = 0.01), algorithm = "brute-force", data = train_data),
      Exponential.HD8 = nls2(Volume ~ a * exp(Diameter * b), start = list(a = 1, b = 0.01), algorithm = "brute-force", data = train_data)
    )

    # Initialize data frames to store results and model summaries
    results_train <- data.frame(Model = character(0), RMSE = numeric(0), MAE = numeric(0), AIC = numeric(0), PER = numeric(0))
    results_test <- data.frame(Model = character(0), RMSE = numeric(0), MAE = numeric(0), AIC = numeric(0), PER = numeric(0))
    model_summaries <- list()

    # Function to calculate RMSE
    calculate_RMSE <- function(model, data) {
      predictions <- predict(model, newdata = data)
      sqrt(mean((data$Volume - predictions)^2))
    }

    # Function to calculate MAE
    calculate_MAE <- function(model, data) {
      predictions <- predict(model, newdata = data)
      mean(abs(data$Volume - predictions))
    }

    # Function to calculate AIC
    calculate_AIC <- function(model) {
      AIC(model)
    }

    # Function to calculate Prediction Error Rate (PER)
    calculate_PER <- function(model, data) {
      predictions <- predict(model, newdata = data)
      sum(abs(data$Volume - predictions) > 0.1) / nrow(data)
    }

    # Fit and evaluate models
    for (model_name in names(models)) {
      model <- models[[model_name]]

      # Calculate metrics for training data
      rmse_train <- calculate_RMSE(model, train_data)
      mae_train <- calculate_MAE(model, train_data)
      aic_train <- calculate_AIC(model)
      per_train <- calculate_PER(model, train_data)

      # Calculate metrics for testing data
      rmse_test <- calculate_RMSE(model, test_data)
      mae_test <- calculate_MAE(model, test_data)
      aic_test <- calculate_AIC(model)
      per_test <- calculate_PER(model, test_data)

      # Store results
      results_train <- rbind(results_train, data.frame(Model = model_name, RMSE = rmse_train, MAE = mae_train, AIC = aic_train, PER = per_train))
      results_test <- rbind(results_test, data.frame(Model = model_name, RMSE = rmse_test, MAE = mae_test, AIC = aic_test, PER = per_test))

      # Store model summaries
      model_summaries[[model_name]] <- summary(model)
    }

    # Select the best models based on RMSE, MAE, AIC, and PER (lower is better)
    best_rmse_model <- results_test[which.min(results_test$RMSE), ]
    best_mae_model <- results_test[which.min(results_test$MAE), ]
    best_aic_model <- results_test[which.min(results_test$AIC), ]
    best_per_model <- results_test[which.min(results_test$PER), ]

    # Print the best models
    message("Best RMSE Model:", best_rmse_model$Model, "RMSE:", best_rmse_model$RMSE, "\n")
    message("Best MAE Model:", best_mae_model$Model, "MAE:", best_mae_model$MAE, "\n")
    message("Best AIC Model:", best_aic_model$Model, "AIC:", best_aic_model$AIC, "\n")
    message("Best PER Model:", best_per_model$Model, "PER:", best_per_model$PER, "\n")

    # Create a data frame for test results
    test_results <- results_test %>% gather(Metric, Value, -Model)

    # Create a bar plot with a title reflecting the type of split
    plot_title <- paste("Split Type:", split_name)
    plot <- ggplot(test_results, aes(x = Model, y = Value, fill = Metric)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      facet_wrap(~Metric, scales = "free_y") +
      labs(x = "Model", y = "Value", fill = "Metric", title = plot_title) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Store results for this split and cross-validation method
    result_name <- paste(split_name, model_name, sep = "_")
    all_results[[result_name]] <- list(
      results_train = results_train,
      results_test = results_test,
      model_summaries = model_summaries,
      plot = plot


    )
  }

  # Return all results
  return(all_results)
}
