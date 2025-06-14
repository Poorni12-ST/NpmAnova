#----------------------compute_derivatives

#' Internal function for the library NpmAnova
#'
#' @description Computes the residuals, Hessian (second-order derivative),
#' and adjusted residuals (z-values) used for the model updates.
#'
#' @param y_true A numeric vector of true response values.
#' @param y_pred A numeric vector of predicted values.
#' @return A list containing:
#'   \item{G}{Residuals (y_pred - y_true).}
#'   \item{H}{Hessian (vector of 1s, assuming a simple squared error loss).}
#'   \item{z}{Adjusted residuals (negative residuals divided by Hessian).}
#' @keywords internal


compute_derivatives <- function(y_true, y_pred) {
  if (length(y_true) != length(y_pred)) {
    stop("y_true and y_pred must have the same length.")
  }

  residuals <- y_pred - y_true
  hessian <- rep(1, length(y_true))  # Assuming a constant Hessian of 1
  z <- -residuals / hessian

  return(list(G = residuals, H = hessian, z = z))
}

#----------------------fit_tree

#' Internal function for the library NpmAnova
#'
#' @description Fits a one-level decision tree (stump) for a given variable, computing
#' the mean residuals and sum of squared errors (SSE) for binary splits.
#'
#' @param data A data frame containing the dataset.
#' @param variable A string specifying the binary predictor variable to split on.
#' @param z A numeric vector of adjusted residuals from the model.
#' @return A list containing:
#'   \item{variable}{The variable used for splitting.}
#'   \item{split_1_mean}{Mean residual for group where variable == 1.}
#'   \item{split_0_mean}{Mean residual for group where variable == 0.}
#'   \item{sse}{Sum of squared errors for the split.}
#' @keywords internal

fit_tree <- function(data, variable, z) {
  if (!variable %in% names(data)) {
    stop("Variable not found in data.")
  }
  if (!is.numeric(z) || length(z) != nrow(data)) {
    stop("z must be a numeric vector of the same length as data rows.")
  }

  split_1 <- data[[variable]] == 1
  if (sum(split_1) == 0 || sum(!split_1) == 0) {
   stop("The variable does not contain both 0 and 1 values.")
  }

  mean_1 <- mean(z[split_1])
  mean_0 <- mean(z[!split_1])
  sse <- sum((z[split_1] - mean_1)^2) + sum((z[!split_1] - mean_0)^2)

  return(list(variable = variable, split_1_mean = mean_1, split_0_mean = mean_0, sse = sse))
}

#----------------------select_best_variable

#' Internal function for the library NpmAnova
#'
#' @description Identifies the variable that minimizes the sum of squared errors (SSE)
#' when used as a binary split in the tree.
#'
#' @param data A data frame containing the dataset.
#' @param variables A character vector of variable names to evaluate.
#' @param z A numeric vector of adjusted residuals from the model.
#' @return A string representing the variable that results in the smallest SSE.
#' @keywords internal

select_best_variable <- function(data, variables, z) {
  if (!is.character(variables) || length(variables) == 0) {
    stop("variables must be a non-empty character vector.")
  }
  if (!is.numeric(z) || length(z) != nrow(data)) {
    stop("z must be a numeric vector of the same length as data rows.")
  }

  sse_values <- sapply(variables, function(var) fit_tree(data, var, z)$sse, USE.NAMES = TRUE)
  best_var <- variables[which.min(sse_values)]

  return(best_var)
}

#----------------------update_predictions

#' Internal function for the library NpmAnova
#'
#' @description Updates the current model predictions using the learned effects
#' from the fitted tree.
#'
#' @param data A data frame containing the dataset.
#' @param current_pred A numeric vector of current predicted values.
#' @param tree A list from `fit_tree()` containing the variable, split means, and SSE.
#' @param learning_rate A numeric value controlling how much to update predictions.
#' @return A numeric vector of updated predictions.
#' @keywords internal

update_predictions <- function(data, current_pred, tree, learning_rate) {
  if (!"variable" %in% names(tree) || !"split_1_mean" %in% names(tree) || !"split_0_mean" %in% names(tree)) {
    stop("Invalid tree structure. Ensure tree is from fit_tree().")
  }
  if (!is.numeric(current_pred) || length(current_pred) != nrow(data)) {
    stop("current_pred must be a numeric vector with length equal to number of rows in data.")
  }
  if (!is.numeric(learning_rate) || learning_rate <= 0) {
    stop("learning_rate must be a positive number.")
  }

  split_1 <- data[[tree$variable]] == 1
  current_pred[split_1] <- current_pred[split_1] + learning_rate * tree$split_1_mean
  current_pred[!split_1] <- current_pred[!split_1] + learning_rate * tree$split_0_mean

  return(current_pred)
}


