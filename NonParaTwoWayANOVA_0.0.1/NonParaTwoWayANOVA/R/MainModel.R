#' Train a model with main and interaction effects
#'
#' @param data_encoded A data frame containing the dataset.
#' @param main_vars A vector of column names representing the main effect variables.
#' @param interaction_vars A vector of column names representing the interaction effect variables.
#' @param response A string specifying the name of the response variable.
#' @param learning_rate A numeric value controlling the learning rate.
#' @param max_iter_main An integer specifying the maximum iterations for main effects (default = 1000).
#' @param max_iter_interaction An integer specifying the maximum iterations for interaction effects (default = 1000).
#' @param patience An integer specifying the early stopping patience (default = 5).
#' @return A data frame containing the estimated effects for main and interaction variables.
#' @keywords internal

MainModel = function(data_encoded,
                        main_vars,
                        interaction_vars,
                        response,
                        learning_rate = NULL,
                        max_iter_main = 1000,
                        max_iter_interaction = 1000,
                        patience = 5){

  # Split data into training and validation sets
  train_indices <- sample(1:nrow(data_encoded), 0.7 * nrow(data_encoded))
  train_data <- data_encoded[train_indices, ]
  valid_data <- data_encoded[-train_indices, ]
  train_pred <- rep(mean(train_data[[response]]), nrow(train_data))
  valid_pred <- rep(mean(valid_data[[response]]), nrow(valid_data))

  best_loss <- Inf
  best_iter <- 0
  best_iter_interaction <- 0
  main_effects_L1 <- setNames(rep(0, length(main_vars)), main_vars)
  main_effects_L0 <- setNames(rep(0, length(main_vars)), main_vars)
  interaction_effects_L1 <- setNames(rep(0, length(interaction_vars)), interaction_vars)
  interaction_effects_L0 <- setNames(rep(0, length(interaction_vars)), interaction_vars)

  loss_history_interaction <- c()
  loss_history_main <- c()

  main_effects_L1_history <- data.frame(iter = integer(), variable = character(), effect = numeric())
  main_effects_L0_history <- data.frame(iter = integer(), variable = character(), effect = numeric())

  interaction_effects_L1_history <- data.frame(iter = integer(), variable = character(), effect = numeric())
  interaction_effects_L0_history <- data.frame(iter = integer(), variable = character(), effect = numeric())

  # Boosting for main effects
  for (iter in 1:max_iter_main) {
    z <- compute_derivatives(train_data[[response]], train_pred)$z
    best_main_var <- select_best_variable(train_data, main_vars, z)

    tree <- fit_tree(train_data, best_main_var, z)
    train_pred <- update_predictions(train_data, train_pred, tree, learning_rate)
    valid_pred <- update_predictions(valid_data, valid_pred, tree, learning_rate)

    main_effects_L1[best_main_var] <- main_effects_L1[best_main_var] + learning_rate * tree$split_1_mean
    main_effects_L0[best_main_var] <- main_effects_L0[best_main_var] + learning_rate * tree$split_0_mean

    main_effects_L1_history <- rbind(main_effects_L1_history, data.frame(iter = iter, variable = best_main_var, effect = main_effects_L1[best_main_var]))
    main_effects_L0_history <- rbind(main_effects_L0_history, data.frame(iter = iter, variable = best_main_var, effect = main_effects_L0[best_main_var]))

    current_loss <- mean((valid_data[[response]] - valid_pred)^2)
    loss_history_main <- c(loss_history_main, current_loss)
    if (current_loss < best_loss){
      best_loss <- current_loss
      best_iter_main <- iter
      best_loss_main <- best_loss
    }else if (iter - best_iter >= patience) break
  }

  # Boosting for interaction effects
  if(length(interaction_vars)!=0){
    for (iter in 1:max_iter_interaction) {
      z <- compute_derivatives(valid_data[[response]], valid_pred)$z
      best_interaction_var <- select_best_variable(valid_data, interaction_vars, z)
      tree <- fit_tree(valid_data, best_interaction_var, z)
      valid_pred <- update_predictions(valid_data, valid_pred, tree, learning_rate)

      interaction_effects_L1[best_interaction_var] <- interaction_effects_L1[best_interaction_var] + learning_rate * tree$split_1_mean
      interaction_effects_L0[best_interaction_var] <- interaction_effects_L0[best_interaction_var] + learning_rate * tree$split_0_mean

      interaction_effects_L1_history <- rbind(interaction_effects_L1_history,
                                              data.frame(iter = iter, variable = best_interaction_var, effect = interaction_effects_L1[best_interaction_var]))
      interaction_effects_L0_history <- rbind(interaction_effects_L0_history,
                                              data.frame(iter = iter, variable = best_interaction_var, effect = interaction_effects_L0[best_interaction_var]))

      parent_main_vars <- unlist(strsplit(best_interaction_var, "__"))
      if (parent_main_vars[1] %in% names(main_effects_L1)) {
        main_effects_L1[parent_main_vars[1]] <- main_effects_L1[parent_main_vars[1]] - learning_rate * tree$split_1_mean / 2
        main_effects_L0[parent_main_vars[1]] <- main_effects_L0[parent_main_vars[1]] - learning_rate * tree$split_0_mean / 2
        main_effects_L1_history <- rbind(main_effects_L1_history,
                                         data.frame(iter = iter, variable = parent_main_vars[1], effect = main_effects_L1[parent_main_vars[1]]))

        main_effects_L0_history <- rbind(main_effects_L0_history,
                                         data.frame(iter = iter, variable = parent_main_vars[1], effect = main_effects_L0[parent_main_vars[1]]))
      }
      #else{ cat("main_effect_",parent_main_vars[1],"not_updated\n")}

      if(parent_main_vars[2] %in% names(main_effects_L1)){
        main_effects_L1[parent_main_vars[2]] <- main_effects_L1[parent_main_vars[2]] - learning_rate * tree$split_1_mean / 2
        main_effects_L0[parent_main_vars[2]] <- main_effects_L0[parent_main_vars[2]] - learning_rate * tree$split_0_mean / 2
        main_effects_L1_history <- rbind(main_effects_L1_history,
                                         data.frame(iter = iter, variable = parent_main_vars[2], effect = main_effects_L1[parent_main_vars[2]]))
        main_effects_L0_history <- rbind(main_effects_L0_history,
                                         data.frame(iter = iter, variable = parent_main_vars[2], effect = main_effects_L0[parent_main_vars[2]]))
      }
      #else{ cat("main_effect_",parent_main_vars[2],"not_updated\n")}

      current_loss <- mean((valid_data[[response]] - valid_pred)^2)
      loss_history_interaction <- c(loss_history_interaction, current_loss)

      if (current_loss < best_loss){
        best_loss <- current_loss
        best_iter_interaction <- iter
        best_loss_interaction <- best_loss
      }else if (iter - best_iter_interaction >= patience) break
    }
    # Final model
    final_model <- data.frame(
      row.names = c("0", "1"),
      sapply(main_vars, function(var) c(main_effects_L0[var], main_effects_L1[var])),
      sapply(interaction_vars, function(var) c(interaction_effects_L0[var], interaction_effects_L1[var])),
      global_mean = c(mean(train_data[[response]]), NA)
    )

  }else{
    # Final model
    final_model <- data.frame(
      row.names = c("0", "1"),
      sapply(main_vars, function(var) c(main_effects_L0[var], main_effects_L1[var])),
      global_mean = c(mean(train_data[[response]]), NA)
    )
  }
  return(final_model)
}


# model01<- MainModel(df_encoded, main_vars, interaction_vars, response, learning_rate = 0.3, patience = 15)

