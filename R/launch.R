#' Launch Kronos scRT processing
#'
#' Main function to launch Kronos scRT processing in your browser.
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

#' Launch Kronos scRT tools
#'
#' Main function to launch Kronos  RT tools in your browser.
#'
#' @return Launches the shiny application
#' @importFrom shiny runApp
#'
#' @export
#'
#'
#'
scRT_tools <- function() {
  # Run the application
  dir_app = system.file('Kronos_scRT_tools', package = 'Kronos.scRT')
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

#' Launch Kronos RT tools
#'
#' Main function to launch Kronos RT tools in your browser.
#'
#' @return Launches the shiny application
#' @importFrom shiny runApp
#'
#' @export
#'
#'
RT_tools <- function() {
  # Run the application
  dir_app = system.file('Kronos_RT_tools', package = 'Kronos.scRT')
  shiny::runApp(appDir = dir_app, launch.browser = T)
}


#' Info to set Kronos scRT for command line
#'
#'@param os
#'
#' @export
#'
#'

shell_interface <-function(os=c('auto','unix','window')){
  #find dir Kronos
  dir_app = system.file(file.path('Kronos_scRT','Kronos.R'), package = 'Kronos.scRT')
  #check platform and suggest operation
  print("add the following line to your profile:")
  if((.Platform$OS.type=='unix' & os[1]=='auto') | os[1]=='unix'){
        print(paste0("alias Kronos='Rscript ", dir_app,"'"))
  }else if((.Platform$OS.type=='windows' & os[1]=='auto') | os[1]=='windows'){
    print(paste('New-Alias -Name Kronos -Value Rscript.exe',dir_app))
  }
}

