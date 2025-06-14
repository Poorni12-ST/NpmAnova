#' Create Main Model and Compute P-Values
#'
#' @description
#'
#' This function takes user input (dataset, variables, and hyperparameters),
#' validates the inputs, create a model, and computes p-values for feature importance.
#'
#' @param df_encoded A data frame containing the dataset with one hot encoded(0 and 1) variables.
#' @param main_vars A character vector of main effect variable names.
#' @param response A string specifying the name of the response variable.
#' @param learning_rate A numeric vector for the learning rates (if NULL, (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7) will be is used).
#' @param max_iter_main An integer specifying the max iterations for main effects (default = 1000).
#' @param max_iter_interaction An integer specifying the max iterations for interaction effects (default = 1000).
#' @param patience An integer specifying early stopping patience (default = 5).
#' @param n_permutations An integer for the number of permutations in p-value calculation (default = 100).
#' @param compute_p_values A logical value (TRUE/FALSE) indicating whether to run permutation testing (default = TRUE).
#' @return A list containing:
#'   \item{model}{The trained model.}
#'   \item{p_values}{A list of p-values for main and interaction variables (if `compute_p_values = TRUE`).}
#' @export
#'
#' @examples
#'library(caret)
#'library(NonParaTwoWayANOVA)

#'# Create Example Data set
#'data_generation <- function(n) {
#'  # Data Generation Function
#'
#'  # Two categorical variables with 2 levels each
#'  factorA <- factor(sample(c("A1", "A2"), n, replace = TRUE))
#'  factorB <- factor(sample(c("B1", "B2"), n, replace = TRUE))
#'
#'  # Generate interaction effects and response variable
#'  mu <- 50  # Global mean
#'  effectA <- c(A1 = 1, A2 = 0)
#'  effectB <- c(B1 = 3, B2 = -10)
#'  interaction_effect <- matrix(c(2, -2, -2, 2), nrow = 2, byrow = TRUE,
#'                               dimnames = list(c("A1", "A2"), c("B1", "B2")))
#'
#'  # Generate response variable (dependent variable)
#'  y <- mu + effectA[factorA] + effectB[factorB] +
#'    mapply(function(a, b) interaction_effect[a, b], factorA, factorB) +
#'    rexp(n,0.2)  # Adding random noise
#'
#'  # Create data frame
#'  df <- data.frame(factorA, factorB, y)
#'  # Apply one-hot encoding
#'  dummies <- dummyVars(~ factorA + factorB , data = df)
#'  encoded_data <- predict(dummies, newdata = df)
#'  head(encoded_data)
#'
#'  # Convert to data frame
#'  response ="y"
#'  encoded_df <- as.data.frame(encoded_data)
#'  df_encoded <- cbind(df[response], encoded_df)
#'
#'  df_encoded=df_encoded[,c("y", "factorA.A1", "factorB.B1")]
#'
#'  head(df_encoded)
#'
#'  main_vars <- c( "factorA.A1", "factorB.B1")
#'
#'  ncol(as.data.frame(df_encoded[, main_vars]))
#'  return(df_encoded)
#'}
#'
#'df_encoded <- data_generation(n = 100)
#'main_vars <- c("factorA.A1", "factorB.B1")
#'response <- "y"
#'set.seed(100)
#'# call NpmAnova with proper arguments
#'result <- NpmAnova(df_encoded,
#'                   main_vars,
#'                   response,
#'                   compute_p_values = TRUE)
#'print(result)


NpmAnova <- function(df_encoded, main_vars, response,
                      learning_rate = NULL, max_iter_main = 1000, max_iter_interaction = 1000,
                      patience = 5, n_permutations = 100, compute_p_values = TRUE) {

  # 1. Validate dataset
  if (!is.data.frame(df_encoded)) stop("Error: df_encoded must be a data frame.")
  if (!is.character(response)) stop(paste("Error: Response variable should be a string specifying the name of the response variable.."))
  if (!(response %in% names(df_encoded))) stop(paste("Error: Response variable", response, "not found in dataset."))
  if (!is.numeric(df_encoded[[response]])) stop("Error: Response variable must be numeric.")

  if (is.null(main_vars)) {
    stop(paste("Error: Main variables missing"))
  }
  if (ncol(as.data.frame(df_encoded[, main_vars]))<2) {
    stop(paste("Error: need at least 2 Main variables"))
  }

  # 2. Validate main variables
  if (!all(main_vars %in% names(df_encoded))) {
    missing_vars <- main_vars[!main_vars %in% names(df_encoded)]
    stop(paste("Error: The following main variables are missing from the dataset:", paste(missing_vars, collapse = ", ")))
  }

  # 3. Validate or Generate Interaction Variables
  interaction_vars = NULL
  if (!is.null(interaction_vars)) {
    incorrect_interactions <- interaction_vars[!grepl("^.+__.+$", interaction_vars)]
    if (length(incorrect_interactions) > 0) {
      stop(paste("Error: Interaction variables must be in 'var1__var2' format. Incorrect:", paste(incorrect_interactions, collapse = ", ")))
    }
  } else {
    interaction_vars <- c()
    for (i in 1:(length(main_vars) - 1)) {
      for (j in (i + 1):length(main_vars)) {
        new_var <- paste(main_vars[i], main_vars[j], sep = "__")
        df_encoded[[new_var]] <- df_encoded[[main_vars[i]]] * df_encoded[[main_vars[j]]]
        interaction_vars <- c(interaction_vars, new_var)
      }
    }
    message("✅ Interaction variables generated: ", paste(interaction_vars, collapse = ", "))
  }

  # 4. Validate that categorical variables are one-hot encoded
  for (var in c(main_vars, interaction_vars)) {
    unique_vals <- unique(df_encoded[[var]])
    if (!all(unique_vals %in% c(0, 1))) {
      stop(paste("Error: Variable", var, "is not one-hot encoded. It should contain only 0s and 1s."))
    }
  }

  # Check if dataset is large enough for 5-fold cross-validation
  if (nrow(df_encoded) < 50) {  # At least 2 observations per fold (5-fold)
    stop("Error: Not enough observations")
  }

  # 5. Determine Learning Rate Using Cross-Validation (if not provided)
  if (is.null(learning_rate)) {
    learning_rate <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
  }
  message("🔍 Running 5-fold cross-validation to determine the best learning rate...")

  # Learning rates to test
  best_rate <- learning_rate[1]
  best_loss <- Inf

  # Create 5-fold indices
  # Generate folds with a retry mechanism
  set.seed(42)  # For reproducibility
  max_attempts <- patience  # Limit attempts to prevent infinite looping
  attempt <- 1

  repeat {
    folds <- sample(rep(1:5, length.out = nrow(df_encoded)))  # Randomly assign folds

    if (validate_folds(df_encoded, folds, c(main_vars, interaction_vars), max_attempts, attempt)) {
      break  # Valid folds found, exit loop
    }
    if (attempt >= max_attempts) {
      stop("Error: Could not generate valid 5-fold splits. The dataset is too small or unbalanced.")
    }
    attempt <- attempt + 1
  }
  avg_losses <- numeric(length(learning_rate))  # Store average loss for each rate

  # progress bar for CV
  pb <- txtProgressBar(min = 0, max = length(learning_rate), style = 3)

  # Loop through each learning rate
  for (lr_idx in seq_along(learning_rate)) {
    rate <- learning_rate[lr_idx]
    fold_losses <- numeric(5)

    # Perform 5-fold cross-validation
    for (fold in 1:5) {
      train_data <- df_encoded[folds != fold, ]
      valid_data <- df_encoded[folds == fold, ]

      # Train model on training fold
      model_temp <- MainModel(train_data, main_vars, interaction_vars, response,
                              learning_rate = rate, max_iter_main, max_iter_interaction, patience)

      # Compute validation loss
      valid_pred <- predicted_values(valid_data, c(main_vars, interaction_vars), model_temp)
      fold_losses[fold] <- mean((valid_data[[response]] - valid_pred)^2)  # MSE loss

      setTxtProgressBar(pb, lr_idx)

    }

    # Compute average loss across folds
    avg_losses[lr_idx] <- mean(fold_losses)

  }
  close(pb)  # Close progress bar

  # Select the best learning rate
  best_rate <- learning_rate[which.min(avg_losses)]
  learning_rate <- best_rate

  message("✅ Best learning rate selected: ", learning_rate)

  # 6. Train the Main Model
  message("🚀 Training the Main Model...")
  model <- MainModel(df_encoded, main_vars, interaction_vars, response,
                     learning_rate, max_iter_main, max_iter_interaction, patience)
  message("✅ Model training complete!")

  # 7. Compute p-values using Permutation Testing (Optional)
  p_values <- list()
  if (compute_p_values) {
    message("🔬 Computing p-values via permutation testing...")

    # Permutation test for main variables
    for (k in seq_along(main_vars)) {
      main_var_inp <- main_vars[-k]  # Remove one main variable
      interaction_var_inp <- interaction_vars  # Keep all interactions
      p_values[[main_vars[k]]] <- perform_permutation_test(model, df_encoded, response, main_vars, interaction_vars, main_var_inp, interaction_var_inp, n_permutations, learning_rate)
    }

    # Permutation test for interaction variables
    for (k in seq_along(interaction_vars)) {
      interaction_var_inp <- interaction_vars[-k]  # Remove one interaction variable
      main_var_inp <- main_vars  # Keep all main effects
      p_values[[interaction_vars[k]]] <- perform_permutation_test(model, df_encoded, response, main_vars, interaction_vars, main_var_inp, interaction_var_inp, n_permutations, learning_rate)
    }
    message("✅ P-values computed!")
  } else {
    message("⚠️ Skipping permutation testing (p-value calculation disabled).")
  }

  # Return Model and P-values
  return(list(model = model, p_values = if (compute_p_values) p_values else NULL))
}

# Function to check if each variable in every fold contains both 0 and 1
#' @keywords internal
validate_folds <- function(data, folds, variables, max_attempts, attempt) {
  for (fold in unique(folds)) {
    fold_data <- data[folds == fold, ]
    for (var in variables) {
      unique_vals <- unique(fold_data[[var]])
      if (length(unique_vals) < 2) {
        return(FALSE)  # Fold does not contain both 0 and 1
      }
    }
  }
  return(TRUE)  # All folds are valid
}
