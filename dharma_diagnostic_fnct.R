##Function to run DHARMA diagnostics


check_lmer_assumptions_dharma <- function(model, model_name = "Model") {
  
  # Load required libraries
  library(DHARMa)
  library(lme4)
  
  cat("=== DHARMa MODEL ASSUMPTION CHECKS ===\n")
  cat("Model:", model_name, "\n\n")
  
  # Create DHARMa residuals
  simulated_residuals <- simulateResiduals(fittedModel = model, plot = FALSE)
  
  # Plot all diagnostic plots
  plot(simulated_residuals, main = paste("DHARMa Residuals:", model_name))
  
  # Individual tests
  cat("\n=== STATISTICAL TESTS ===\n")
  
  # Test for uniformity (overall model fit)
  uniformity_test <- testUniformity(simulated_residuals, plot = FALSE)
  cat("Uniformity Test (KS test):\n")
  cat("D =", round(uniformity_test$statistic, 4), ", p-value =", 
      format.pval(uniformity_test$p.value), "\n")
  if(uniformity_test$p.value < 0.05) {
    cat("** Model may have poor fit **\n\n")
  } else {
    cat("Model fit appears adequate\n\n")
  }
  
  # Test for outliers
  outlier_test <- testOutliers(simulated_residuals, plot = FALSE)
  cat("Outlier Test:\n")
  cat("p-value =", format.pval(outlier_test$p.value), "\n")
  if(outlier_test$p.value < 0.05) {
    cat("** Potential outliers detected **\n\n")
  } else {
    cat("No significant outliers detected\n\n")
  }
  
  # Test for dispersion
  dispersion_test <- testDispersion(simulated_residuals, plot = FALSE)
  cat("Dispersion Test:\n")
  cat("p-value =", format.pval(dispersion_test$p.value), "\n")
  if(dispersion_test$p.value < 0.05) {
    cat("** Over/under-dispersion detected **\n\n")
  } else {
    cat("Dispersion appears appropriate\n\n")
  }
  
  # Test for zero-inflation (if applicable)
  tryCatch({
    zero_inflation_test <- testZeroInflation(simulated_residuals, plot = FALSE)
    cat("Zero-inflation Test:\n")
    cat("p-value =", format.pval(zero_inflation_test$p.value), "\n")
    if(zero_inflation_test$p.value < 0.05) {
      cat("** Zero-inflation detected **\n\n")
    } else {
      cat("No zero-inflation detected\n\n")
    }
  }, error = function(e) {
    cat("Zero-inflation test not applicable for this model\n\n")
  })
  
  # Return DHARMa object for further analysis
  invisible(simulated_residuals)
}