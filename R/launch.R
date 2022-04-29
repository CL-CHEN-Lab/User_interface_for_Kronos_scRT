#' Launch Kronos scRT
#'
#' Main function to launch Kronos scRT in your browser.
#'
#' @return Launches the shiny application
#' @importFrom shiny runApp
#'
#' @export
#'
#'
#'
Processing <- function() {
  # Run the application
  dir_app = system.file('Kronos_processing', package = 'Kronos.scRT')
  shiny::runApp(appDir = dir_app, launch.browser = T)
}

#' Launch Kronos scRT post-processing
#'
#' Main function to launch Kronos scRT post-processing in your browser.
#'
#' @return Launches the shiny application
#' @importFrom shiny runApp
#'
#' @export
#'
#'
#'
Post_processing <- function() {
  # Run the application
  dir_app = system.file('Kronos_post_processing', package = 'Kronos.scRT')
  shiny::runApp(appDir = dir_app, launch.browser = T)
}


#' Launch Kronos scRT pre-processing
#'
#' Main function to launch Kronos scRT pre-processing in your browser.
#'
#' @return Launches the shiny application
#' @importFrom shiny runApp
#'
#' @export
#'
#'
Pre_processing <- function() {
  # Run the application
  dir_app = system.file('Kronos_pre_processing', package = 'Kronos.scRT')
  shiny::runApp(appDir = dir_app, launch.browser = T)
}
