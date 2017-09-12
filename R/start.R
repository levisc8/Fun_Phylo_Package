#' @importFrom stats rnorm
.onLoad <- function(libname, pkgname){
  packageStartupMessage(paste0('test', round(stats::rnorm(1), 3)))
}