#' Formula interface for vblogit used by the logi.engine

vblogit.fmla <- function(formula, data, subset, weights, verbose=FALSE, epsilon=0.01, ...) {
  #' mimic glm fit?
  #environment(formula) <- environment()
  #f0 <- glm(formula, data = data, family = binomial(), subset = subset, weights=weights, method = "model.frame")
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  yX <- eval(mf, parent.frame())
  n <- ncol(yX)
  nams <- colnames(yX)
  offn <- "offset(-log(.logi.B))"
  yn <- ".logi.Y"
  offset <- yX[[offn]]
  y <- yX$.logi.Y
  varnames <- setdiff(nams, c(offn, yn))
  X <- cbind(1, as.matrix(yX[varnames]))
  colnames(X)[1] <- "(Intercept)"
  Vnames <- colnames(X)
  
  #' then we fit:
  fit <- vblogit(y, X, offset, verb=verbose, eps=epsilon, ...)
  #'
  names(fit$coefficients) <- names(fit$coef) <- Vnames
  #' add some variables to conform to summary.ppm
  fit$se <- sqrt(diag(as.matrix(fit$S)))
  fit$formula <- formula
  fit$method <- "vblogit"
  fit$terms <- attr(yX, "terms")
  fit
}