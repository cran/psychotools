## workhorse function
mptmodel <- function(y, weights = NULL, spec, treeid = NULL,
  optimargs = list(control = list(reltol = .Machine$double.eps^(1/1.2),
                                  maxit = 1000),
                   init = NULL),
  start = NULL, vcov = TRUE, estfun = FALSE, ...)
{
  ## main arguments
  if(missing(y)) stop("response missing")
  stopifnot(class(spec) == "mptspec")

  ## sanity checking of data
  if(is.null(dim(y))) y <- matrix(y, nrow=1L, dimnames=list(NULL, names(y)))
  y <- as.matrix(y)
  nsubj <- nrow(y)
  # mob() adds names to y
  # if(!is.null(dnam <- colnames(y)) & !is.null(snam <- names(spec$prob)))
  #   if(!all(snam == dnam)) warning("variable names do not match")
  if(NCOL(y) != length(spec$prob))
    stop("number of response categories and model equations do not match")

  ## weights
  if(is.null(weights)) weights <- 1L
  weights <- rep(weights, length.out = nsubj)
  nsubj <- sum(weights > 0L)

  ## treeid
  tid <- if(length(treeid) == NCOL(y))
           treeid
         else if(length(spec$treeid) == NCOL(y))            # read from spec
           spec$treeid
         else if(!is.null(colnames(y))) {                   # before 1st dot
           tid <- gsub("([^.]+)\\..*", "\\1", colnames(y))
           factor(tid, levels = unique(tid))                # keep order
         } else
           rep(1, NCOL(y))
  tid <- factor(tid)  # make sure treeid is always factor

  ## for fitting only sums are needed
  ysum <- colSums(y * weights)
  
  ## determine number of parameters and starting values
  if(is.null(start)) {
    start <- spec$par[is.na(spec$par)]  # FIX ME: is.na still necessary?
    start[] <- 0
  } else {
    ## do sanity checking of starting values/names/etc.
    if(is.null(names(start))) names(start) <- names(spec$par[is.na(spec$par)])
    start <- qlogis(start)  # logit transform
  }

  ## set up log-likelihood and gradient
  nll <- function(par) -sum(ysum * log(spec$par2prob(plogis(par))))
  grad <- function(par) {
    yp <- drop(ysum/spec$par2prob(plogis(par)))
    dp <- spec$par2deriv(plogis(par))$deriv
    -drop(dp %*% yp) * dlogis(par)  # FIX ME: dlogis(par) optional? Ask Z.
  }
  
  ## call optim()
  if (!is.null(optimargs$init)) {
    start[] <- qlogis(optimargs$init)  # Workaround: override start, for
  }                                    #   setting start vals from mpttree()
  optimargs$init <- NULL
  optArgs <- list(par=start, fn=nll, gr=grad, method="BFGS")
  optArgs <- c(optArgs, as.list(optimargs))
  opt <- do.call(optim, optArgs)
  coef <- opt$par
  loglik <- -opt$value
  pcat <- spec$par2prob(plogis(coef))
 
  snam <- if(!is.null(names(spec$prob))) names(spec$prob)       # from spec
          else if(!is.null(colnames(y))) colnames(y)            # from data
          else ave(as.character(tid), tid,                      # guess
                   FUN = function(a) paste(a, seq_along(a), sep="."))
  ncat    <- table(tid)
  ntrees  <- length(ncat)
# n       <- setNames(tapply(ysum, tid, sum)[as.character(tid)], snam)
  nbytree <- tapply(ysum, tid, sum)
  n       <- setNames(nbytree[as.character(tid)], snam)
  fitted  <- n*pcat
  G2      <- 2*sum(ysum*log(ysum/fitted), na.rm=TRUE)
  df      <- sum(ncat - 1) - length(coef)
  gof     <- c(G2=G2, df=df, pval = pchisq(G2, df, lower.tail=FALSE))
  if(df == 0 && abs(G2) < sqrt(.Machine$double.eps)) {
    gof["G2"] <- 0
    gof["pval"] <- 1
  }

  rval <- list(
    y = y,
    coefficients = coef,
  # loglik = sum(lfactorial(nbytree)) - sum(lfactorial(ysum)) + loglik,
    loglik = loglik,
    nobs = nsubj,
    npar = length(start),
    weights = if(isTRUE(all.equal(weights, rep(1, nsubj)))) NULL else weights,
    fitted = fitted,
    goodness.of.fit = gof,
    ntrees = ntrees,
    n = n,
    pcat = setNames(pcat, snam),
    treeid = tid,
    spec = spec,
    optim = opt,
    ysum = setNames(ysum, snam)
  )
  class(rval) <- "mptmodel"
  if(estfun) rval$estfun <- estfun.mptmodel(rval)
  return(rval)
}


logLik.mptmodel <- function(object, ...)
  structure(object$loglik, df = object$npar, class = "logLik")


coef.mptmodel <- function(object, logit = FALSE, ...){
  coef <- object$coefficients
  if (logit) {
    nm <- paste0("logit(", names(coef), ")")
    setNames(coef, nm)
  } else {
    plogis(coef)
  }
}


estfun.mptmodel <- function(x, logit = TRUE, ...)
{
  dp <- t(x$spec$par2deriv(coef(x, logit=FALSE))$deriv)
  # ef <- t(sapply(seq_len(nrow(x$y)),   # does not work for single-par models
  #   function(i) colSums(x$y[i, ]/x$pcat * dp)))
  ef <- x$y %*% (dp/x$pcat)
  if (logit) ef <- t(dlogis(coef(x, logit=TRUE)) * t(ef))
  dimnames(ef) <- list(rownames(x$y), names(coef(x, logit=logit)))
  if(!is.null(x$weights)) ef <- ef * weights  
  ef
}


vcov.mptmodel <- function(object, logit = TRUE, what = c("vcov", "fisher"),
                     ...){
  what <- match.arg(what)
  coef <- coef(object, logit=logit)

  # Negative Hessian (estimated information) on probability scale.
  # %*% is slightly faster than sum( * )
  H <- function(par, y = object$ysum, spec = object$spec){
    pp  <- spec$par2prob(par)
    yp  <- drop(y/pp)
    dp  <- spec$par2deriv(par)$deriv
    npar <- length(par)
    H <- matrix(NA, npar, npar)
    for (i in seq_len(npar))
      for (j in i:npar)
          H[i, j] <- yp %*% (dp[i, ]*dp[j, ]/pp)
        # H[i, j] <- sum(y*dp[i, ]*dp[j, ]/pp^2)
        # H[i, j] <- sum(yp*dp[i, ]*dp[j, ]/pp)
        # H[i, j] <- sum(yp/pp * dp[i, ]*dp[j, ])
    H[lower.tri(H)] <- t(H)[lower.tri(H)]
    dimnames(H) <- list(names(par), names(par))
    H
  }

  if (logit) {  # delta method
    # Note that plogis(x)*(1 - plogis(x)) == dlogis(x)
    # Automatic dimnames for tcrossprod(...)
    # Conjecture: the smaller the matrix the more reliable its inverse
    # dhessian <- diag(1/(plogis(coef) * (1 - plogis(coef))), length(coef))
    # dhessian %*% solve(H(plogis(coef))) %*% dhessian
    # solve(tcrossprod(dlogis(coef)) * H(plogis(coef)))  # possibly singular
    if (what == "vcov") 1/tcrossprod(dlogis(coef)) * solve(H(plogis(coef)))
    else                  tcrossprod(dlogis(coef)) *       H(plogis(coef))
  } else {
    if (what == "vcov") solve(H(coef))
    else                      H(coef)
  }
}


## Based on stats::confint.default
confint.mptmodel <- function(object, parm, level = 0.95, logit = TRUE,
                             ...){
  cf <- coef(object, logit=logit)
  pnames <- names(cf)
  if (missing(parm)) 
      parm <- pnames
  else if (is.numeric(parm)) 
      parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100*a, trim=TRUE, scientific=FALSE, digits=3), "%")
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object, logit=logit)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}


print.mptmodel <- function(x, digits = max(3, getOption("digits") - 3),
                      logit=FALSE, ...){
  cat("\nMultinomial processing tree (MPT) models\n\n")
  cat("Parameter estimates:\n")
  print.default(format(coef(x, logit=logit), digits=digits), print.gap=2,
                quote = FALSE)
  G2   <- x$goodness.of.fit[1]
  df   <- x$goodness.of.fit[2]
  pval <- x$goodness.of.fit[3]
  cat("\nGoodness of fit (2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits), ", p = ",
      format(pval, digits=digits), "\n", sep="")
  cat("\n")
  invisible(x)
}


## basic plot of parameter estimates
plot.mptmodel <- function(x, logit = FALSE, xlab = "Parameter",
                          ylab = "Estimate", ...){
  coef <- coef(x, logit=logit)
  plot(coef, axes=FALSE, xlab = xlab, ylab = ylab, ...)
  axis(1, seq_along(coef), names(coef))
  axis(2)
  box()
}


nobs.mptmodel <- function(object, ...) object$nobs


summary.mptmodel <- function(object, ...){
  x <- object
  coef <- coef(x, logit=TRUE)
  pcoef <- coef(x, logit=FALSE)

  ## Catch vcov error, so there are at least some NA's in the summary
  s.err <- tryCatch(sqrt(diag(vcov(x, logit=TRUE))),
                    error = function(e) rep(NA, length(coef)))

  tvalue <- coef / s.err
  pvalue <- 2 * pnorm(-abs(tvalue))
  dn <- c("Estimate", "Logit Estim.", "Std. Error")
  coef.table <- cbind(pcoef, coef, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(pcoef), c(dn, "z value", "Pr(>|z|)"))

  aic <- AIC(x)
  ans <- list(ntrees=x$ntrees, coefficients=coef.table, aic=aic,
              gof=x$goodness.of.fit)  # , X2=sum(resid(x, "pearson")^2))
  class(ans) <- "summary.mptmodel"
  return(ans)
}


print.summary.mptmodel <- function(x, digits = max(3, getOption("digits") - 3),
                                   cs.ind = 2:3, ...){
  cat("\nCoefficients:\n")
  printCoefmat(x$coef, digits=digits, cs.ind=cs.ind, ...)
  cat("\nLikelihood ratio G2:", format(x$gof[1], digits=digits), "on",
    x$gof[2], "df,", "p-value:", format(x$gof[3], digits=digits), "\n")
  cat(# "Pearson X2: ", format(x$X2, digits=digits), ",    ",
      "AIC: ", format(x$aic, digits=max(4, digits + 1)), sep="")
  cat("\n")
  cat("Number of trees:", x$ntrees, "\n")
  invisible(x)
}


deviance.mptmodel <- function(object, ...) object$goodness.of.fit["G2"]


predict.mptmodel <- function(object, newdata = NULL, ...){
  if(is.null(newdata)) fitted(object)
  else {
    stopifnot(length(newdata) == length(object$n))
    object$pcat * newdata  # newdata: total response count per tree
  }
}


## model specification
mptspec <- function(..., .replicates = NULL, .restr = NULL)
{
  ## non-standard evaluation of arguments
  spec <- match.call()

  replicates <- spec$.replicates
  if(!is.null(replicates)) {
    if(!is.numeric(replicates)) stop(".replicates must be a single number")
  }

  restr <- spec$.restr
  if(!is.null(restr)) {
    if(as.character(restr[[1L]]) != "list") stop(".restr must be list")
    restr1 <- restr
    restrcl <- sapply(restr1[-1L], class)
  # restr <- sapply(restr1[-1L], deparse)
    restr <- sapply(restr1[-1L], function(s) paste(deparse(s), collapse=""))
    restr <- paste(names(restr), " = ",
                   ifelse(restrcl == "numeric", "", "expression("),
                   restr,
                   ifelse(restrcl == "numeric", "", ")[[1L]]"),
                   collapse = ", ")
  }

  # Remove .replicates and .restr from call (if included) and
  # (further below, after pre-specified models) turn into list of characters
  spec$.replicates <- NULL
  spec$.restr <- NULL
  spec <- as.list(spec[-1L])                  # exclude function name
  treeid <- NULL

  if (is.character(whichmod <- spec[[1]])) {  # pre-specified models
    prespec <- switch(EXPR = whichmod,
      "1HT" = expression(
        "1.1" = r + (1 - r)*b,
        "1.2" = (1 - r)*(1 - b),
        "2.1" = b,
        "2.2" = 1 - b
      ),
      "2HT" = expression(
        "1.1" = r + (1 - r)*b,
        "1.2" = (1 - r)*(1 - b),
        "2.1" = (1 - d)*b,
        "2.2" = (1 - d)*(1 - b) + d
      ),
      "PairAsso" = expression(
        "1.1" = p*q*r,
        "1.2" = p*q*(1 - r),
        "1.3" = p*(1 - q)*r,
        "1.4" = (1 - p) + p*(1 - q)*(1 - r)
      ),
      "proCNI" = expression(
        "1.1" = C + (1 - C)*(1 - N)*(1 - I),
        "1.2" = (1 - C)*N + (1 - C)*(1 - N)*I,
        "2.1" = C + (1 - C)*(1 - N)*J,
        "2.2" = (1 - C)*N + (1 - C)*(1 - N)*(1 - J),
        "3.1" = (1 - C)*(1 - N)*(1 - I),
        "3.2" = C + (1 - C)*N + (1 - C)*(1 - N)*I,
        "4.1" = (1 - C)*(1 - N)*J,
        "4.2" = C + (1 - C)*N + (1 - C)*(1 - N)*(1 - J)
      ),
      "prospec" = expression(
        "1.1" = C1*P*(1 - M1)*(1 - g) + C1*(1 - P) +
                (1 - C1)*P*(1 - M1)*(1 - g)*c + (1 - C1)*(1 - P)*c,
        "1.2" = (1 - C1)*P*(1 - M1)*(1 - g)*(1 - c) +
                (1 - C1)*(1 - P)*(1 - c),
        "1.3" = C1*P*M1 + C1*P*(1 - M1)*g + (1 - C1)*P*M1 +
                (1 - C1)*P*(1 - M1)*g,
        "2.1" = (1 - C2)*P*(1 - M1)*(1 - g)*c + (1 - C2)*(1 - P)*c,
        "2.2" = C2*P*(1 - M1)*(1 - g) + C2*(1 - P) +
                (1 - C2)*P*(1 - M1)*(1 - g)*(1 - c) +
                (1 - C2)*(1 - P)*(1 - c),
        "2.3" = C2*P*M1 + C2*P*(1 - M1)*g + (1 - C2)*P*M1 +
                (1 - C2)*P*(1 - M1)*g,
        "3.1" = C1*P*M2 + C1*P*(1 - M1)*(1 - g) + C1*(1 - P) +
                (1 - C1)*P*M2*c + (1 - C1)*P*(1 - M1)*(1 - g)*c +
                (1 - C1)*(1 - P)*c,
        "3.2" = (1 - C1)*P*M2*(1 - c) + (1 - C1)*P*(1 - M2)*(1 - g)*(1 - c) +
                (1 - C1)*(1 - P)*(1 - c),
        "3.3" = C1*P*(1 - M2)*g + (1 - C1)*P*(1 - M2)*g,
        "4.1" = (1 - C2)*P*M2*c + (1 - C2)*P*(1 - M2)*(1 - g)*c +
                (1 - C2)*(1 - P)*c,
        "4.2" = C2*P*M2 + C2*P*(1 - M2)*(1 - g) + C2*(1 - P) +
                (1 - C2)*P*M2*(1 - c) + (1 - C2)*P*(1 - M2)*(1 - g)*(1 - c) +
                (1 - C2)*(1 - P)*(1 - c),
        "4.3" = C2*P*(1 - M2)*g + (1 - C2)*P*(1 - M2)*g
      ),
      "rmodel" = expression(
        "1.1" = b,
        "1.2" = 1 - b,
        "2.1" = g,
        "2.2" = 1 - g,
        "3.1" = r*a + (1 - r)*b*a,
        "3.2" = r*(1 - a) + (1 - r)*(1 - b)*(1 - a),
        "3.3" = (1 - r)*(1 - b)*a,
        "3.4" = (1 - r)*b*(1 - a)
      ),
      "SourceMon" = expression(
        "1.1" = D1*d1 + D1*(1 - d1)*g + (1 - D1)*b*g,
        "1.2" = D1*(1 - d1)*(1 - g) + (1 - D1)*b*(1 - g),
        "1.3" = (1 - D1)*(1 - b),
        "2.1" = D2*(1 - d2)*g + (1 - D2)*b*g,
        "2.2" = D2*d2 + D2*(1 - d2)*(1 - g) + (1 - D2)*b*(1 - g),
        "2.3" = (1 - D2)*(1 - b),
        "3.1" = b*g,
        "3.2" = b*(1 - g),
        "3.3" = 1 - b
      ),
      "SR" = expression(
        "1.1" = c*r,
        "1.2" = (1 - c)*u^2,
        "1.3" = 2*(1 - c)*u*(1 - u),
        "1.4" = c*(1 - r) + (1 - c)*(1 - u)^2,
        "2.1" = u,
        "2.2" = 1 - u
      ),
      "SR2" = expression(
        "1.1" = c*r,
        "1.2" = (1 - c)*u^2,
        "1.3" = 2*(1 - c)*u*(1 - u),
        "1.4" = c*(1 - r) + (1 - c)*(1 - u)^2
      ),
      "WST" = expression(
        ".0000" = (1 - a)*(1 - P)*(1 - p)*(1 - Q)*(1 - q),
        ".0001" = a*c*(1 - d)*(1 - sb)*i +
                  (1 - a)*(1 - P)*(1 - p)*(1 - Q)*q,
        ".0010" = a*c*(1 - d)*sb*i +
                  (1 - a)*(1 - P)*(1 - p)*Q*(1 - q),
        ".0011" = a*(1 - c)*(1 - x)*(1 - d)*i +
                  (1 - a)*(1 - P)*(1 - p)*Q*q,
        ".0100" = a*c*d*(1 - sf)*i +
                  (1 - a)*(1 - P)*p*(1 - Q)*(1 - q),
        ".0101" = a*(1 - c)*x*(1 - sfb)*i +
                  (1 - a)*(1 - P)*p*(1 - Q)*q,
        ".0110" = a*c*d*(1 - sf)*(1 - i) +
                  a*c*(1 - d)*sb*(1 - i) +
                  (1 - a)*(1 - P)*p*Q*(1 - q),
        ".0111" = (1 - a)*(1 - P)*p*Q*q,
        ".1000" = a*c*d*sf*i +
                  (1 - a)*P*(1 - p)*(1 - Q)*(1 - q),
        ".1001" = a*c*d*sf*(1 - i) +
                  a*c*(1 - d)*(1 - sb)*(1 - i) +
                  (1 - a)*P*(1 - p)*(1 - Q)*q,
        ".1010" = a*(1 - c)*x*sfb*i +
                  (1 - a)*P*(1 - p)*Q*(1 - q),
        ".1011" = (1 - a)*P*(1 - p)*Q*q,
        ".1100" = a*(1 - c)*(1 - x)*d*i +
                  (1 - a)*P*p*(1 - Q)*(1 - q),
        ".1101" = (1 - a)*P*p*(1 - Q)*q,
        ".1110" = (1 - a)*P*p*Q*(1 - q),
        ".1111" = a*(1 - c)*x*sfb*(1 - i) +
                  a*(1 - c)*x*(1 - sfb)*(1 - i) +
                  a*(1 - c)*(1 - x)*d*(1 - i) +
                  a*(1 - c)*(1 - x)*(1 - d)*(1 - i) +
                  (1 - a)*P*p*Q*q
      ),
      NULL  # model not available
    )
    if(is.null(prespec))
      stop("'...' has to be either an expression or one of:\n",
           "  '1HT', '2HT', 'PairAsso', 'proCNI', 'prospec', 'rmodel',\n",
           "  'SourceMon', 'SR', 'SR2', 'WST'.\n")
  }

  if(exists("prespec", inherits = FALSE)) spec <- prespec

  ## Get treeid from names
  if(!is.null(names(spec))) {
    nm <- do.call(rbind, strsplit(names(spec), "\\."))  # treeid.catid
    if(any(nm[, 1] == "")) nm[, 1] <- "1"
    treeid <- nm[, 1]  # check only if there are .replicates
  }

  ## Replicates?
  if (!is.null(replicates) && replicates > 1) {
    if(is.null(names(spec)) || ncol(nm) != 2 || nrow(nm) != length(spec))
      stop('.replicates require expressions named "x.y"',
           ' where "x" is the treeid \n',
           '  See help("mptspec") for examples.\n')
    ntrees <- length(unique(treeid))
    treeid <- if(anyNA(suppressWarnings(as.numeric(treeid)))) {
      paste0(rep(treeid, replicates),                  # character treeid
             rep(seq_len(replicates), each = length(treeid)))
    } else {
      rep(as.numeric(treeid), replicates) +            # numeric treeid
      rep(seq(0, ntrees*(replicates - 1), ntrees),
          each = nrow(nm))
    }
    pd <- getParseData(parse(text=spec, keep.source=TRUE))
    pat <- paste0("(",
            paste(unique(pd$text[pd$token == "SYMBOL"]), collapse="|"), ")")
    newspec <- NULL
    for (i in seq_len(replicates))
      newspec <- c(newspec, gsub(pat, paste0("\\1", i), spec))
    spec <- setNames(parse(text=newspec),
                     paste(treeid, nm[, 2], sep="."))
  }

  ## list of strings
  spec <- lapply(spec, function(s) paste(deparse(s), collapse=""))

  ## substitute restrictions
  if(!is.null(restr)) {
    spec <- lapply(spec, function(s) {
      s <- sprintf("substitute(%s, list(%s))", s, restr)
      deparse(eval(parse(text = s)))
    })
  }

  ## parsed expressions (list of expressions)
  if(!is.null(restr)) restr <- lapply(restr1[-1L], as.expression)
  prob <- lapply(spec, function(s) parse(text=s, keep.source=TRUE))
  if(is.null(treeid) && !is.null(names(prob)))
    treeid <- gsub("([^.]+)\\..*", "\\1", names(prob))  # guess from names

  ## extract the parameters
  pars <- unique(unlist(lapply(prob, function(e) {
    pd <- getParseData(e)
    pd$text[pd$token == "SYMBOL"]                     # get parameter names
  })))
  pars <- structure(rep.int(NA_real_, length(pars)), .Names = pars)

  # ## use .pars to fix parameters or starting values or so
  # if(!is.null(.pars)) {
  #   if(is.list(.pars)) .pars <- do.call("c", .pars)
  #   if(is.null(names(.pars))) stop(".pars must be named list or vector")
  #   pars[names(.pars)] <- .pars
  # }

  ## compute class probabilities
  par2prob <- function(par) {
    ## get all parameters via lexical scoping
    pars <- pars
    
    ## replace NA parameters
    if(sum(is.na(pars)) != length(par))
      stop("numbers of parameters do not match")
    pars[is.na(pars)] <- par
    pars <- as.list(pars)
    
    ## compute probabilities
    rval <- sapply(prob, eval, pars)
    names(rval) <- names(prob)
    return(rval)
  }

  ## derivatives, deriv3() instead of deriv() for second derivatives
  deriv <- lapply(prob, deriv3, names(pars))
  names(deriv) <- names(prob)

  par2deriv <- function(par) {
    ## get all parameters via lexical scoping
    pars <- pars
    
    ## replace NA parameters: FIX ME still needed?
    na_pars <- is.na(pars)
    if(sum(na_pars) != length(par))
      stop("numbers of parameters do not match")
    pars[na_pars] <- par
    pars <- as.list(pars)
    
    ## compute first derivatives
    deriv1 <- rbind(
               sapply(deriv, function(ex) attr(eval(ex, pars), "gradient")))
    rownames(deriv1) <- names(pars)
    deriv1 <- deriv1[na_pars, , drop = FALSE]  # Jacobian

    ## compute second derivatives
    # deriv2 <- lapply(deriv, function(ex) attr(eval(ex, pars), "hessian"))
    # deriv2 <- array(unlist(deriv2),
    #                 c(length(pars), length(pars), length(prob)), 
    #                 list(names(pars), names(pars), names(prob)))
    # deriv2 <- deriv2[na_pars, na_pars, , drop = FALSE]

    # list(deriv = deriv1, deriv2 = deriv2)  # return 1st and 2nd derivatives
    list(deriv = deriv1)                     # return 1st derivative
  }

  retval <- list(
    par2prob = par2prob,
    par2deriv = par2deriv,
    prob = prob,
    deriv = deriv,
    par = pars,
    replicates = replicates,
    restr = restr,
    treeid = if(!is.null(treeid))
      factor(treeid, levels = unique(treeid))  # keep order
  )
  class(retval) <- "mptspec"
  retval
}


## Apply replicates/restrictions to existing mptspec object
update.mptspec <- function(object, .replicates = NULL, .restr = NULL,
                           ...){
  spec <- match.call()
  replicates <- spec$.replicates
  restr <- spec$.restr

  spec <- unlist(object$prob)
  if(!is.null(replicates)) {
    if(!is.numeric(replicates)) stop(".replicates must be a single number")
    spec$.replicates <- replicates
  }
  if(!is.null(restr)){
    if(as.character(restr[[1L]]) != "list") stop(".restr must be list")
    spec$.restr <- restr
  }
  do.call(mptspec, spec)
}


## Print model equations
print.mptspec <- function(x, ...){
  print(unlist(x$prob), ...)
}

