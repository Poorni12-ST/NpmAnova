#' Internal function for the library NpmAnova (Perform a Permutation Test to Compute P-Values)
#'
#' @description Evaluates the significance of an effect by running a permutation test
#' and comparing the model performance with and without the variable.
#'
#' @param org_model A trained model obtained from `MainModel()`.
#' @param df_encoded A data frame containing the dataset with encoded variables.
#' @param response A string specifying the response variable name.
#' @param main_vars_Full A character vector of all main effect variables.
#' @param interaction_vars_Full A character vector of all interaction effect variables.
#' @param main_var_inp A character vector of main effect variables to be included in the reduced model.
#' @param interaction_var_inp A character vector of interaction effect variables to be included in the reduced model.
#' @param n An integer specifying the number of permutations (default = 100).
#' @param learning_rate A numeric vector for the learning rates.
#' @return A numeric p-value representing the statistical significance of the removed variable.
#' @keywords internal

perform_permutation_test <- function(org_model, df_encoded, response,
                                     main_vars_Full, interaction_vars_Full,
                                     main_var_inp, interaction_var_inp, n = 100, learning_rate) {

  # heck if df_encoded is a valid data frame
  if (!is.data.frame(df_encoded)) {
    stop("df_encoded must be a data frame.")
  }

  # heck if response variable exists and is numeric
  if (!(response %in% names(df_encoded))) {
    stop(paste("Response variable", response, "not found in df_encoded."))
  }
  if (!is.numeric(df_encoded[[response]])) {
    stop("Response variable must be numeric.")
  }

  # nsure org_model is a valid list
  if (!is.list(org_model) || length(org_model) == 0) {
    stop("org_model must be a valid model list.")
  }

  # nsure variable inputs are character vectors
  if (!is.character(main_vars_Full) || !is.character(interaction_vars_Full)) {
    stop("main_vars_Full and interaction_vars_Full must be character vectors.")
  }

  if (!is.character(main_var_inp) || !is.character(interaction_var_inp)) {
    stop("main_var_inp and interaction_var_inp must be character vectors.")
  }

  # nsure n is a positive integer
  n <- floor(n)
  if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive integer.")
  }

  # Compute original model predictions
  all_vars <- c(main_vars_Full, interaction_vars_Full)
  pp_org <- predicted_values(df_encoded, all_vars, org_model)

  # Compute residuals and R² for the original model
  residual <- df_encoded[[response]] - pp_org
  sse <- sum(residual^2)
  tss <- sum((df_encoded[[response]] - mean(df_encoded[[response]]))^2)
  r2_org <- 1 - sse / tss

  # Initialize permutation results
  df_encoded$y_s_r <- df_encoded[[response]]
  r2_per <- numeric(n)

  # Fit the model without one factor
  re_N <- MainModel(df_encoded, main_var_inp, interaction_var_inp, response, learning_rate = learning_rate, patience = 15)

  # Compute residuals for reduced model
  all_vars <- c(main_var_inp, interaction_var_inp)
  pp_N <- predicted_values(df_encoded, all_vars, re_N)
  residual_N <- df_encoded[[response]] - pp_N

  for (i in 1:n) {
    # Shuffle residuals and create a new response variable
    df_encoded$y_s_r <- pp_N + sample(residual_N, replace = FALSE)

    # Fit the full model again with permuted response values
    re_s_r <- MainModel(df_encoded, main_vars_Full, interaction_vars_Full, "y_s_r", learning_rate = learning_rate, patience = 15)
    pp_N2 <- predicted_values(df_encoded, all_vars, re_s_r)

    # Compute R² for permuted model
    sse <- sum((df_encoded$y_s_r - pp_N2)^2)
    tss <- sum((df_encoded$y_s_r - mean(df_encoded$y_s_r))^2)
    r2_per[i] <- 1 - sse / tss

    cat(sprintf("\r🔄 Processing permutation %d/%d...", i, n))
    flush.console()
  }

  # Compute p-value
  p <- (sum(r2_per > r2_org) + 1) / (n + 1)

  return(p)
}
