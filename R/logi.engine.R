#' Override spatstat's logi.engine
#' 
#' @import spatstat vblogistic
#' @export

logi.engine_with_vb <- function (Q, trend = ~1, interaction, ..., covariates = NULL, 
          correction = "border", rbord = reach(interaction), covfunargs = list(), 
          allcovar = FALSE, vnamebase = c("Interaction", "Interact."), 
          vnameprefix = NULL, justQ = FALSE, savecomputed = FALSE, 
          precomputed = NULL, VB=FALSE) {
  if (is.null(trend)) 
    trend <- ~1
  if (is.null(interaction)) 
    interaction <- Poisson()
  want.trend <- !identical.formulae(trend, ~1)
  want.inter <- !is.poisson(interaction)
  correction <- pickoption("correction", correction, c(border = "border", 
                                                       periodic = "periodic", isotropic = "isotropic", Ripley = "isotropic", 
                                                       trans = "translate", translate = "translate", translation = "translate", 
                                                       none = "none"))
  if (correction == "border") {
    check.1.real(rbord, "In ppm")
    explain.ifnot(rbord >= 0, "In ppm")
  }
  else rbord <- 0
  if (!missing(vnamebase)) {
    if (length(vnamebase) == 1) 
      vnamebase <- rep.int(vnamebase, 2)
    if (!is.character(vnamebase) || length(vnamebase) != 
          2) 
      stop("Internal error: illegal format of vnamebase")
  }
  if (!is.null(vnameprefix)) {
    if (!is.character(vnameprefix) || length(vnameprefix) != 
          1) 
      stop("Internal error: illegal format of vnameprefix")
  }
  if (inherits(Q, "ppp")) {
    Xplus <- Q
    Q <- quadscheme.logi(Xplus, ...)
    D <- Q$dummy
    Dinfo <- Q$param
  }
  else if (checkfields(Q, c("data", "dummy"))) {
    Xplus <- Q$data
    D <- Q$dummy
    Dinfo <- Q$param
    if (is.null(Dinfo)) {
      Dinfo <- list(how = "given", rho = npoints(D)/(area.owin(D) * 
                                                       markspace.integral(D)))
    }
    Q <- quadscheme.logi(Xplus, D)
  }
  else stop("Format of object Q is not understood")
  if (justQ) 
    return(Q)
  extraargs <- list(covfunargs = covfunargs, allcovar = allcovar, 
                    vnamebase = vnamebase, vnameprefix = vnameprefix)
  extraargs <- append(extraargs, list(...))
  if (correction == "border" && Dinfo$how == "grid") {
    Dbord <- D[bdist.points(D) >= rbord]
    Dinfo$rho <- npoints(Dbord)/(eroded.areas(as.owin(Dbord), 
                                              rbord) * markspace.integral(Dbord))
  }
  rho <- Dinfo$rho
  B <- list(...)$Barker
  if (is.null(B)) 
    B <- 1
  B <- B * rho
  Dinfo <- append(Dinfo, list(B = B))
  Dinfo <- append(Dinfo, list(extraargs = extraargs))
  Wplus <- as.owin(Xplus)
  nXplus <- npoints(Xplus)
  U <- superimpose(Xplus, D, W = Wplus, check = FALSE)
  E <- cbind(1:nXplus, 1:nXplus)
  computed <- if (savecomputed) 
    list(X = Xplus, Q = Q, U = U)
  else list()
  if (want.trend) {
    tvars <- variablesinformula(trend)
    wantxy <- c("x", "y") %in% tvars
    wantxy <- wantxy | rep.int(allcovar, 2)
    cvdf <- data.frame(x = U$x, y = U$y)[, wantxy, drop = FALSE]
    if (!is.null(covariates)) {
      df <- mpl.get.covariates(covariates, U, "quadrature points", 
                               covfunargs)
      cvdf <- cbind(cvdf, df)
    }
    wantmarks <- "marks" %in% tvars
    if (wantmarks) 
      cvdf <- cbind(cvdf, marks = marks(U))
  }
  else cvdf <- NULL
  if (!is.null(ss <- interaction$selfstart)) 
    interaction <- ss(Xplus, interaction)
  V <- evalInteraction(Xplus, U, E, interaction, correction, 
                       precomputed = precomputed, savecomputed = savecomputed)
  if (!is.matrix(V)) 
    stop("evalInteraction did not return a matrix")
  if (savecomputed) 
    computed <- append(computed, attr(V, "computed"))
  IsOffset <- attr(V, "IsOffset")
  if (is.null(IsOffset)) 
    IsOffset <- rep.int(FALSE, ncol(V))
  if (ncol(V) > 0) {
    Vnames <- colnames(V)
    if (is.null(Vnames)) {
      nc <- ncol(V)
      Vnames <- if (nc == 1) 
        vnamebase[1]
      else paste(vnamebase[2], 1:nc, sep = "")
      colnames(V) <- Vnames
    }
    else if (!is.null(vnameprefix)) {
      Vnames <- paste(vnameprefix, Vnames, sep = "")
      colnames(V) <- Vnames
    }
  }
  else Vnames <- character(0)
  glmdata <- as.data.frame(V)
  if (!is.null(cvdf)) 
    glmdata <- cbind(glmdata, cvdf)
  ok <- if (correction == "border") 
    (bdist.points(U) >= rbord)
  else rep.int(TRUE, npoints(U))
  KEEP <- if (ncol(V) > 0) 
    matrowall(V != -Inf)
  else rep.int(TRUE, npoints(U))
  ok <- ok & KEEP
  wei <- c(rep.int(1, npoints(Xplus)), rep.int(B/rho, npoints(D)))
  resp <- c(rep.int(1, npoints(Xplus)), rep.int(0, npoints(D)))
  glmdata <- cbind(glmdata, .logi.Y = resp, .logi.B = B, .logi.w = wei, 
                   .logi.ok = ok)
  trendpart <- paste(as.character(trend), collapse = " ")
  fmla <- paste(".logi.Y ", trendpart)
  if (want.inter) {
    VN <- Vnames
    if (any(IsOffset)) 
      VN[IsOffset] <- paste("offset(", VN[IsOffset], ")", 
                            sep = "")
    fmla <- paste(c(fmla, VN), collapse = "+")
  }
  fmla <- paste(fmla, "offset(-log(.logi.B))", sep = "+")
  fmla <- as.formula(fmla)
  .logi.B <- B
  .logi.w <- wei
  .logi.ok <- ok
  .logi.Y <- resp
  fit <- if(VB) vblogit.fmla(fmla, data = glmdata, subs = .logi.ok, w = .logi.w, ...)
         else glm(fmla, data = glmdata, family = binomial(), subset = .logi.ok, 
                  weights = .logi.w)
  co <- coef(fit)
  fitin <- fii(interaction, co, Vnames, IsOffset)
  maxlogpl <- logLik(fit) + sum(ok * resp * log(B))
  spv <- package_version(versionstring.spatstat())
  the.version <- list(major = spv$major, minor = spv$minor, 
                      release = spv$patchlevel, date = "$Date: 2014/04/17 08:45:16 $")
  fit <- list(method = "logi", fitter = ifelse(VB, "vblogit", "glm"), projected = FALSE, 
              coef = co, trend = trend, interaction = interaction, 
              Q = Q, correction = correction, rbord = rbord, version = the.version, 
              fitin = fitin, maxlogpl = maxlogpl, internal = list(Vnames = Vnames, se=if(VB)sqrt(diag(as.matrix(fit$S)))else NULL,
                                                                  IsOffset = IsOffset, glmdata = glmdata, glmfit = fit, 
                                                                  logistic = Dinfo, computed = computed))
  class(fit) <- "ppm"
  
  #' add some variables to conform to summary.ppm
  return(fit)
}


