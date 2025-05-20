#' Compute Resultant Ground Reaction Force
#'
#' Calculates the Euclidean resultant magnitude of ground reaction force vectors.
#'
#' @param df A data.frame or tibble with numeric columns \code{Fx}, \code{Fy}, \code{Fz}.
#' @return A numeric vector of resultant magnitudes \eqn{\sqrt{Fx^2 + Fy^2 + Fz^2}}.
#' @export
compute_resultant <- function(df) {
  if (!all(c("Fx","Fy","Fz") %in% colnames(df))) {
    stop("Input must contain Fx, Fy, Fz columns.")
  }
  sqrt(df$Fx^2 + df$Fy^2 + df$Fz^2)
}
