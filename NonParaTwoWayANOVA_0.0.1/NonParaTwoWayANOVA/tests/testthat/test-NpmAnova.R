library(testthat)
library(caret)
library(NonParaTwoWayANOVA)

test_that("NpmAnova runs without errors", {
  # Data Generation Function
  data_generation <- function(n) {

    # Two categorical variables with 2 levels each
    factorA <- factor(sample(c("A1", "A2"), n, replace = TRUE))
    factorB <- factor(sample(c("B1", "B2"), n, replace = TRUE))

    # Generate interaction effects and response variable
    mu <- 50  # Global mean
    effectA <- c(A1 = 1, A2 = 0)
    effectB <- c(B1 = 3, B2 = -10)
    interaction_effect <- matrix(c(2, -2, -2, 2), nrow = 2, byrow = TRUE,
                                 dimnames = list(c("A1", "A2"), c("B1", "B2")))

    # Generate response variable (dependent variable)
    y <- mu + effectA[factorA] + effectB[factorB] +
      mapply(function(a, b) interaction_effect[a, b], factorA, factorB) +
      rexp(n,0.2)  # Adding random noise

    # Create data frame
    df <- data.frame(factorA, factorB, y)
    # Apply one-hot encoding
    dummies <- dummyVars(~ factorA + factorB , data = df)
    encoded_data <- predict(dummies, newdata = df)
    head(encoded_data)

    # Convert to data frame
    response ="y"
    encoded_df <- as.data.frame(encoded_data)
    df_encoded <- cbind(df[response], encoded_df)

    df_encoded=df_encoded[,c("y", "factorA.A1", "factorB.B1")]

    head(df_encoded)

    main_vars <- c( "factorA.A1", "factorB.B1")

    ncol(as.data.frame(df_encoded[, main_vars]))
    return(df_encoded)
  }

  df_encoded <- data_generation(n = 100)
  main_vars <- c("factorA.A1", "factorB.B1")
  response <- "y"
  set.seed(100)
  # call NpmAnova with proper arguments
  result <- NpmAnova(df_encoded,
                     main_vars,
                     response,
                     compute_p_values = TRUE)

  # Assertions
  expect_type(result, "list")
  expect_true("model" %in% names(result))
})
