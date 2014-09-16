#' Some extra info on the VB-ppm fit
#' 
#' @import vblogistic
#' @export
vbsummary <- function(x, ...) {
  y <- x$internal$glmfit
  summary(y)
  invisible(y)
}
