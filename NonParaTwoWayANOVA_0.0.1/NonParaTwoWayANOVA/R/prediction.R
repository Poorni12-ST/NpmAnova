#' Internal function for the library NpmAnova
#'
#' @description Computes predictions for a given dataset using the model
#' containing main and interaction effects.
#'
#' @param df_encoded A data frame containing the dataset with encoded variables.
#' @param variables A character vector of variable names used in the model.
#' @param final_model A data frame containing the estimated effects for main and interaction variables.
#' @return A numeric vector of predicted values.
#' @keywords internal


predicted_values <- function(df_encoded, variables, final_model) {
  if (!all(variables %in% names(df_encoded))) {
    stop("Some variables are missing in the dataset.")
  }
  if (!is.data.frame(final_model) || !"global_mean" %in% colnames(final_model)) {
    stop("Invalid final_model structure.")
  }

  pred <- numeric(nrow(df_encoded))  # Initialize predictions

  for (i in 1:nrow(df_encoded)) {
    y <- 0
    for (j in seq_along(variables)) {
      var_name <- variables[j]
      if (!var_name %in% colnames(final_model)) {
        stop(paste("Variable", var_name, "not found in final_model."))
      }

      if (df_encoded[i, var_name] == 0) {
        x <- final_model["0", var_name]
      } else {
        x <- final_model["1", var_name]
      }
      y <- y + x
    }

    pred[i] <- y + final_model["0", "global_mean"]
  }

  return(pred)
}
