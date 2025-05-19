#' Launch the Shiny App
#'
#' @export
run_app <- function() {
  shiny::runApp(system.file("app", package = "shinyForce3D"))
}
