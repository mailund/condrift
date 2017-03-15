#' Construct a joint drift model from the mean and covariance matrix
#' of allele frequencies.
#'
#' @param mu    The mean allele frequency for each population.
#' @param Sigma The covariance matrix of drift.
#' @return A wrapper around the distribution
#' @export
joint_drift_model <- function(mu, Sigma) {
  structure(list(mu = mu, Sigma = Sigma),
            class = "drift_model")
}

#' Computes the drift model conditional on the allele frequency of one population
#'
#' @param model The shared drift model
#' @param index The population to condition on
#' @param freq  The frequency to fix the allele frequency in \code{index} to
#' @return A shared drift model for the remaining populations
#' @export
condition <- function(model, index, freq) {
  if (is.character(index))
    index <- which(index == rownames(model$Sigma))

  mu_1 <- model$mu[-index]
  mu_2 <- model$mu[index]
  Sigma_11 <- model$Sigma[-index, -index]
  Sigma_12 <- model$Sigma[-index, index]
  Sigma_21 <- model$Sigma[index, -index]
  Sigma_22 <- model$Sigma[index, index]

  X <- Sigma_12 %*% solve(Sigma_22)
  new_mu <- mu_1 + X %*% (freq - mu_2)
  new_Sigma <- Sigma_11 - X %*% Sigma_21

  names(new_mu) <- names(mu_1)
  rownames(new_Sigma) <- colnames(new_Sigma) <- rownames(Sigma_11)

  joint_drift_model(new_mu, new_Sigma)
}
