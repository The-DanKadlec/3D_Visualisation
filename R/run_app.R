#' Launch the Shiny App
#'
#' This function launches the shinyForce3D application.
#'
#' @export
run_app <- function() {
  source(system.file("R", "shinyForce3D.R", package = "shinyForce3D"), local = TRUE)
  shiny::shinyApp(ui = ui, server = server)
}