#' Extract a covariance matrix from a fitted admixture graph
#'
#' From a graph constructed from the \code{admixturegraph} package, and an environment
#' containing edge lengths -- as a fit using the \code{admixturegraph} package would provide
#' -- this function extracts an expected covariance matrix.
#'
#' The construction of the covariance matrix follows Felsenstein's Inferring Phylogenies (2004)
#' chapter 23 and needs a leaf to root the graph in. This must be provided by a parameter.
#'
#' @param graph A graph specified using the \code{admixturegraph} package.
#' @param env   An environment mapping edge names to edge lengths.
#' @param root  A population we use to root the graph.
#' @return The expected covariance matrix given the graph.
#' @export
extract_covariance_matrix <- function(graph, env, root) {

  populations <- graph$leaves[graph$leaves != root]
  Sigma <- matrix(NA, nrow = length(populations), ncol = length(populations))

  for (i in seq_along(populations)) {
    for (j in seq_along(populations)) {
      if (i != j)
        Sigma[i, j] <- eval(admixturegraph::sf3(graph, root, populations[i], populations[j]), env)
    }
  }

  for (i in seq_along(populations)) {
    Sigma[i, i]  <- eval(admixturegraph::sf2(graph, root, populations[i]), env)
  }

  rownames(Sigma) <- colnames(Sigma) <- populations
  return(Sigma)
}
