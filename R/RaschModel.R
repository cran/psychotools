## workhorse fitting function
RaschModel.fit <- function(y, weights = NULL, start = NULL, reltol = 1e-10, 
  deriv = c("sum", "diff", "numeric"), hessian = TRUE, maxit = 100L,
  full = TRUE, gradtol = reltol, iterlim = maxit, ...)
{
  ## argument matching
  if(missing(reltol) && !missing(gradtol) && !is.null(gradtol)) reltol <- gradtol
  if(missing(maxit) && !missing(iterlim) && !is.null(iterlim)) maxit <- iterlim
  deriv <- match.arg(deriv)

  ## original data
  y <- as.matrix(y)
  k <- k_orig <- ncol(y)
  n <- nrow(y)
  if(is.null(colnames(y))) colnames(y) <- paste("Item", gsub(" ", "0", format(1:k)), sep = "")

  ## weights processing
  if(is.null(weights)) weights <- rep.int(1L, n)
  ## data and weights need to match
  stopifnot(length(weights) == n)

  ## omit zero weights
  weights_orig <- weights
  y_orig <- y
  y <- y[weights > 0, , drop = FALSE]
  weights <- weights[weights > 0]
  n <- nrow(y)

  ## all parameters identified?
  if(n < 2) stop("not enough observations")
  cm <- colMeans(y, na.rm = TRUE)
  status <- as.character(cut(cm, c(-Inf, 1/(2 * n), 1 - 1/(2 * n), Inf), labels = c("0", "0/1", "1")))
  status[is.na(status)] <- "NA"
  status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  ident <- status == "0/1"
  names(status) <- colnames(y)

  ## just estimate identified parameters
  y_orig <- y_orig[,ident, drop = FALSE]
  y <- y[,ident, drop = FALSE]
  k <- ncol(y)
  y_na <- is.na(y)
  any_y_na <- any(y_na)

  if(!any_y_na) {
    ## compute likelihood/gradient/hessian on aggregated data
  
    ## data statistics
    cs <- colSums(y * weights)
    rs <- rowSums(y)
    rf <- as.vector(tapply(weights, factor(rs, levels = 0:k), sum))
    rf[is.na(rf)] <- 0

    ## starting values
    ## contrast: set parameter 1 to zero
    if(is.null(start)) {
      start <- log(sum(weights) - cs) - log(cs) #previously:# -qlogis(cs/sum(weights))
      start <- start[-1] - start[1]
    }
    rf <- rf[-1]
    cs <- cs[-1]

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {
      ## obtain esf and apply contrast
      esf <- elementary_symmetric_functions(c(0, par), order = 0, diff = deriv == "diff")
      g <- esf[[1]][-1]

      ## conditional log-likelihood
      cll <- sum(-cs * par) - sum(rf * log(g))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      return(-cll)
    }

    ## analytical gradient
    agrad <- function(par) {

        ## calculate esf
        esf <- elementary_symmetric_functions(c(0, par), order = 1, diff = deriv == "diff")

        ## calculate gradient
        - colSums(weights * (- y + esf[[2]][rs + 1, , drop = FALSE] / esf[[1]][rs + 1])[,-1, drop = FALSE])
    }

    ## analytical hessian
    ahessian <- function(par, esf) {
      ## obtain esf and apply contrast
      g <- esf[[1]][-1]
      g1 <- esf[[2]][-1, -1, drop = FALSE]
      g2 <- esf[[3]][-1, -1, -1, drop = FALSE]

      ## hessian
      hess <- matrix(0, ncol = k-1, nrow = k-1)
      g1s <- g1/g
      for (q in 1:(k-1)) hess[q,] <- colSums(rf * (g2[,q,]/g - (g1[,q]/g) * g1s))
    
      return(hess)
    }

  } else {
    ## compute likelihood/gradient/hessian on individual data

    ## process NA patterns and calculate static things once
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    na_i <- wi_i <- wi2_i <- cs_i <- rs_i <- rf_i <- k_i <- vector("list", nlevels(na_patterns))
    na_pattern_levels <- levels(na_patterns)

    for (i in seq_along(na_pattern_levels)) {

       ## parse NA pattern
       na_level_i <- na_pattern_levels[i]
       wi_i[[i]] <- as.integer(strsplit(na_level_i, "\r")[[1]])
       wi2_i[[i]]<- if (length(wi_i[[i]]) < 1) 1:k else (1:k)[-wi_i[[i]]]
       k_i[[i]] <- length(wi2_i[[i]])

       ## select subset
       na_i[[i]] <- which(na_patterns == na_level_i)
       if(length(wi_i[[i]]) < 1) y_i <- y[na_i[[i]], , drop = FALSE]
       else y_i <- y[na_i[[i]], -wi_i[[i]], drop = FALSE]
       weights_i <- weights[na_i[[i]]]
       cs_i[[i]] <- colSums(y_i * weights_i)
       rs_i[[i]] <- rowSums(y_i)
       rf_i[[i]] <- as.vector(tapply(weights_i, factor(rs_i[[i]], levels = 0:ncol(y_i)), sum))
       rf_i[[i]][is.na(rf_i[[i]])] <- 0
    }

    ## starting values
    if(is.null(start)) {
      cs <- colSums(y * weights, na.rm = TRUE)
      ws <- colSums(!y_na * weights)
      start <- log(ws - cs) - log(cs) #previously:# -qlogis(cs/ws)
      start <- start[-1] - start[1]
    }

    ## convenience function
    zero_fill <- function(obj, at) {
      if(length(at) < 1) return(obj)
      if(is.null(dim(obj))) {
        rval <- rep.int(0, length(obj) + length(at))
	rval[-at] <- obj
      } else {
        rval <- matrix(0, ncol = ncol(obj) + length(at), nrow = nrow(obj) + length(at))
	rval[-at,-at] <- obj      
      }
      return(rval)
    }

    ## conditional log-likelihood function for NA pattern i
    cll_i <- function (cs_i, par_i, rf_i) {
        sum(-cs_i * par_i) - sum(rf_i * log(elementary_symmetric_functions(par_i, order = 0, diff = deriv == "diff")[[1]]))
    }

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {

      ## initialize return values and extract esf parameters
      cll <- 0
      par_i <- lapply(wi_i, function(x) if (length(x) < 1) c(0, par) else c(0, par)[-x])
      
      ## conditional log-likelihood
      cll <- sum(mapply(cll_i, cs_i, par_i, rf_i, SIMPLIFY = TRUE, USE.NAMES = FALSE))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      ## collect and return
      return(-cll)
    }
 
    ## analytical gradient
    agrad <- function(par) {

      ## initialize return value and esf parameters
      rval <- matrix(0, nrow = n, ncol = k)
      par_i <- lapply(wi_i, function(x) if (length(x) < 1) c(0, par) else c(0, par)[-x])
      esf_i <- mapply(elementary_symmetric_functions, par = par_i, MoreArgs = list(order = 1, diff = deriv == "diff"), SIMPLIFY = FALSE)

      ## loop over observed NA patterns	 
      for(i in seq_along(levels(na_patterns))) {
          rval[na_i[[i]], wi2_i[[i]]] <- weights[na_i[[i]]] * (- y[na_i[[i]], wi2_i[[i]], drop = FALSE] +
          esf_i[[i]][[2]][rs_i[[i]] + 1, , drop = FALSE] / esf_i[[i]][[1]][rs_i[[i]] + 1])
      }
    
      return(- colSums(rval[, -1, drop = FALSE]))
    }

    ## analytical hessian
    ahessian <- function(par, esf) {

      ## set up return value    
      rval <- matrix(0, ncol = k-1, nrow = k-1)

      ## loop over observed NA patterns      
      for(i in seq_along(levels(na_patterns))) {
	
        ## obtain esf
        g_i <- esf[[i]][[1]]
        g1_i <- esf[[i]][[2]]
        g2_i <- esf[[i]][[3]]

        ## hessian
	hess <- matrix(0, nrow = k_i[[i]], ncol = k_i[[i]])
        for (q in 1:k_i[[i]]) hess[q,] <- colSums(rf_i[[i]] * (g2_i[,q,]/g_i - (g1_i[,q]/g_i) * g1_i/g_i))
        rval <- rval + zero_fill(hess, wi_i[[i]])[-1, -1]
      }

      return(rval)
    }

  }
  
  ## optimization
  if(maxit > 0L) {
  opt <- optim(par = start, fn = cloglik, gr = agrad, method = "BFGS",
               hessian = (deriv == "numeric") & hessian, control = list(reltol = reltol, maxit = maxit, ...))
  } else {
    opt <- list(
      estimate = start,
      minimum = cloglik(start),
      hessian = if(deriv != "numeric") ahessian(start, esf) else NULL, ## no numeric Hessian available here
      iterations = 0,
      code = 4)
  }
  
  ## collect and annotate results
  cf <- opt$par
  names(cf) <- colnames(y)[-1]
  if(full) {
    esf <- if(any_y_na) {
      lapply(levels(na_patterns), function(z) {
        wi <- as.integer(strsplit(z, "\r")[[1]])
        cfi <- if(length(wi) < 1) c(0, cf) else c(0, cf)[-wi]
        elementary_symmetric_functions(cfi,
          order = 2 - (deriv == "numeric" | !hessian), 
	  diff = deriv == "diff")
      })
    } else {
      elementary_symmetric_functions(c(0, cf),
        order = 2 - (deriv == "numeric" | !hessian),
        diff = deriv == "diff")
    }
    if(any_y_na) names(esf) <- levels(na_patterns)
  
    if(hessian) {
      vc <- if(deriv == "numeric") opt$hessian else ahessian(cf, esf)
      vc <- solve(vc)
    } else {
      vc <- matrix(NA, nrow = length(cf), ncol = length(cf))
    }
    rownames(vc) <- colnames(vc) <- names(cf)
  } else {
    esf <- NULL
    vc <- NULL
  }

  ## collect, class, and return
  rval <- list(
    coefficients = cf,
    vcov = vc,
    loglik = -opt$value,
    df = k-1,
    data = y_orig,
    weights = if(identical(as.vector(weights_orig), rep(1L, nrow(y_orig)))) NULL else weights_orig,
    n = sum(weights_orig > 0),
    items = status,
    na = any_y_na,
    elementary_symmetric_functions = esf,
    code = opt$convergence,
    iterations = tail(na.omit(opt$counts), 1L),
    reltol = reltol,
    deriv = deriv        
  )
  class(rval) <- "RaschModel"
  return(rval)
}

## methods
coef.RaschModel <- function(object, ...) object$coefficients

vcov.RaschModel <- function(object, ...) object$vcov

logLik.RaschModel <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

weights.RaschModel <- function(object, ...) if(!is.null(object$weights)) object$weights else rep(1, nrow(object$data))

print.RaschModel <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Rasch model difficulty parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

worth.RaschModel <- function(object, difficulty = TRUE, ...) {
  cf <- c(0, object$coefficients)
  if(!difficulty) cf <- -cf
  cf <- cf - mean(cf)
  rval <- structure(rep(NA, length(object$items)), .Names = names(object$items))
  rval[object$items == "0/1"] <- cf
  rval[object$items == "0"] <- if(!difficulty) -Inf else Inf
  rval[object$items == "1"] <- if(!difficulty) Inf else -Inf
  return(rval)
}

summary.RaschModel <- function(object, vcov. = NULL, ...)
{
  ## coefficients
  cf <- coef(object)

  ## covariance matrix
  if(is.null(vcov.)) 
      vc <- vcov(object)
  else {
      if(is.function(vcov.)) vc <- vcov.(object)
        else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2 * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf      
  class(object) <- "summary.RaschModel"
  return(object)
}

print.summary.RaschModel <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat("\nRasch model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if(any(x$items != "0/1")) cat("Excluded items:",
    paste(names(x$items)[x$items != "0/1"], collapse = ", "), "\n\n")

  cat("Difficulty parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in BFGS optimization:", x$iterations, "\n\n")
  invisible(x)
}

plot.RaschModel <- function(x, difficulty = TRUE,
  center = TRUE, index = TRUE, names = NULL, abbreviate = FALSE, ref = TRUE,
  col = cbind("lightgray", "black"), refcol = "lightgray", linecol = "black", lty = 2,
  cex = 1, pch = cbind(19, 1), type = NULL, ylim = NULL, xlab = "Items", ylab = NULL, ...)
{
  ## parameters to be plotted
  cf <- worth(x, difficulty = difficulty)
  cf_ident <- is.finite(cf) & !is.na(cf)
  cf_inf <- cf >= Inf
  cf_ninf <- cf <= -Inf
  if(!center) cf <- cf - (cf[cf_ident])[1]
  cf_ref <- mean(cf[cf_ident])
  ncf <- length(cf)

  ## labeling
  if(is.null(names)) names <- !index
  if(is.character(names)) {
    names(cf) <- names
    names <- TRUE
  }
  if(!names & index) {
    lab <- rep(NA, ncf)
    lab[c(1, ncf)] <- c(1, ncf)
    pr <- pretty(1:ncf)
    pr <- pr[pr > 1 & pr < ncf]
    lab[pr] <- pr    
    names(cf) <- lab
  }

  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(names(cf)))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  names(cf) <- abbreviate(names(cf), abbreviate)

  ## graphical parameter processing  
  if(is.null(type)) type <- if(index) "b" else "p"

  if(NCOL(pch) == 2) {
    pch2 <- pch[,2]
    pch <- pch[,1]
  } else {
    pch2 <- NULL
  }
  if(NCOL(col) == 2) {
    col2 <- col[,2]
    col <- col[,1]
  } else {
    col2 <- NULL
  }
  pch <- rep(pch, length.out = ncf)
  col <- rep(col, length.out = ncf)
  cex <- rep(cex, length.out = ncf)
  pch[!cf_ident] <- NA
  pch2 <- rep(pch2, length.out = ncf)
  col2 <- rep(col2, length.out = ncf)
  if(!is.null(pch2)) pch2[!cf_ident] <- NA
  
  if(is.null(ylim)) ylim <- range(cf[cf_ident])
  ylim <- rep(ylim, length.out = 2)
  if(any(!is.finite(cf))) {
    ydiff <- diff(ylim) * 0.7
    if(index & any(cf_ninf)) ylim[1] <- ylim[1] - ydiff
    if(index & any(cf_inf))  ylim[2] <- ylim[2] + ydiff
  }

  ## substitute non-identified parameters with plottable values
  cf[is.na(cf)] <- cf_ref
  if(index) {
    cf[cf_ninf] <- ylim[1]
    cf[cf_inf] <- ylim[2]
  }

  if(is.null(ylab)) ylab <- paste(if(center) "Centered item" else "Item",
    if(difficulty) "difficulty" else "easiness", "parameters")

  ## raw plot
  ix <- if(index) seq(along = cf) else rep(0, ncf)
  plot(ix, cf, xlab = xlab, ylab = ylab, type = "n", axes = FALSE, ylim = ylim, ...)
  if(ref) abline(h = cf_ref, col = refcol)
  axis(2)
  box()  

  ## actual data
  if(!index & names) {
    text(names(cf), x = ix, y = cf, col = if(!is.null(col2)) col2 else col, ...)
    if(any(!cf_ident)) {
      legend("topright", c("Not identified:", names(cf)[!cf_ident]), bty = "n")
    }
  } else {
    if(type %in% c("l", "b", "o")) lines(ix, cf, type = type, lty = lty, pch = NA, col = linecol)
    lines(ix, cf, type = type, lty = 0, pch = pch, col = col, cex = cex)
    if(type %in% c("b", "p", "o")) points(ix, cf, pch = pch2, col = col2, cex = cex)
    if(index) axis(1, at = ix, labels = names(cf))
    if(type %in% c("l", "b", "o")) {
      if(any(cf_ninf)) for(i in which(cf_ninf)) lines(c(ix[i], ix[i]), c(ylim[1], ylim[1] - 10 * ydiff), type = "l", lty = lty)
      if(any(cf_inf)) for(i in which(cf_inf))   lines(c(ix[i], ix[i]), c(ylim[2], ylim[2] + 10 * ydiff), type = "l", lty = lty)
    }    
  }
}

itempar.RaschModel <- function (object, ref = NULL, vcov = TRUE, ...) {

  ## extract cf and labels, include restricted parameter
  cf <- c(0.00, coef(object))
  m <- length(cf)
  lbs <- names(cf)
  lbs[1] <- ifelse(lbs[1] == "Item02", "Item01", colnames(object$data)[1])

  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else stop("Argument 'ref' can only be a character vector with item labels or a numeric vector with item indices.")

  ## create contrast matrix
  D <- diag(m)
  D[, ref] <- D[, ref] - 1/length(ref)

  ## impose restriction
  cf <- as.vector(D %*% cf)
  names(cf) <- lbs

  ## create adjusted vcov if requested
  if (vcov) {
    vc <- D %*% rbind(0, cbind(0, vcov(object))) %*% t(D)
    colnames(vc) <- rownames(vc) <- lbs
  }

  ## return results
  rval <- structure(cf, class = "itempar", model = "RM", ref = ref, vcov = if (vcov) vc else NULL)
  return(rval)
}

threshold.RaschModel <- function (object, type = c("mode", "median", "mean", "unmodified"), ref = NULL, ...) {
  ## check type, extract item parameters
  type <- match.arg(type)
  ip <- itempar.RaschModel(object, ref = ref, vcov = FALSE)

  ## return requested threshold parameters // in RM: mode = median = mean = unmodified
  class(ip) <- "threshold"
  attr(ip, "type") <- type
  return(ip)
}

personpar.RaschModel <- function (object, ref = NULL, vcov = TRUE, start = NULL, tol = 1e-6, ...) {

  ## extract item parameters and relevant informations #FIXME: constrained parameter included?
  ipar <- itempar.RaschModel(object, ref = ref, vcov = FALSE)
  m <- length(ipar)
  rng <- 1:(m-1)

  ## fetch data, weight processing only to select data, rest is in predict, calculate scores
  if (vcov) {
    y <- object$data[weights(object) > 0, , drop = FALSE]
    cs <- colSums(y)
    rs <- rowSums(y)
    rf <- tabulate(rs, nbins = m - 1)

    ## remove unidentified parameters
    rs <- rs[rf != 0]
    rng <- rng[rf != 0]
    rf <- rf[rf != 0]
  }

  ## start values / range of person parameters
  if(is.null(start)) start <- if(vcov) qlogis(rng/m) else c(-1, 1) * qlogis(1/m * 1e-3) #FIXME: 1e3 enough?

  if (vcov) {
    ## objective function
    cloglik <- function (ppar) {
      - sum(rf * rng * ppar) + sum(cs * ipar) + sum(rf * rowSums(log(1 + outer(exp(ppar), exp(-ipar)))))
    }
    
    ## optimization & target variables
    opt <- nlm(cloglik, start, gradtol = tol, hessian = TRUE, check.analyticals = FALSE, ...)
    ppar <- opt$estimate
    vc <- solve(opt$hessian)
  } else {
    ## iterate over raw scores (from 1 until rmax-1)
    ppar <- sapply(rng, function(rawscore) {
      uniroot(function(ppar) rawscore - sum(plogis(ppar - ipar)), interval = start, tol = tol, ...)$root
    })
  }

  ## return value
  return(structure(ppar, .Names = rng, class = "personpar", model = "RM", vcov = if (vcov) vc else NULL))
}

predict.RaschModel <- function (object, newdata = NULL, type = c("probability", "cumprobability", "mode", "median", "mean"),
                                ref = NULL, ...) {
  
  ## check type, process newdata, if NULL, use person parameters of given model object
  type <- match.arg(type)
  if (is.null(newdata)) {
    rs <- rowSums(object$data, na.rm = TRUE)
    rs <- rs[0 < rs & rs < ncol(object$data)]
    newdata <- personpar.RaschModel(object, ref = ref, vcov = FALSE)[rs]
    names(newdata) <- NULL
  }

  ## calculate probabilities with internal function
  probs <- prm(theta = newdata, beta = itempar.RaschModel(object, ref = ref, vcov = FALSE))

  ## add category zero probabilities for consistency with other predict functions
  if (type %in% c("probability", "cumprobability")) {
    clnms <- colnames(probs)
    rwnms <- rownames(probs)
    nc <- ncol(probs)
    probs0 <- matrix(0, ncol = 2 * nc, nrow = nrow(probs))
    probs0[, seq(from = 2, by = 2, length.out = nc)] <- probs
    probs0[, seq(from = 1, by = 2, length.out = nc)] <- 1 - probs
    if (type == "cumprobability") {
      probs0[, seq(from = 1, by = 2, length.out = nc)] <- probs0[, seq(from = 1, by = 2, length.out = nc)] + probs
    }
    probs <- probs0
    rownames(probs) <- rwnms
    colnames(probs) <- as.vector(t(outer(clnms, c("C0", "C1"), paste, sep = if (type == "probability") "-" else ">=")))
  }

  ## return as requested in type, for RM mode, median, mean is the same
  switch(type,
         "probability" = probs,
         "cumprobability" = probs,
         "mode" = round(probs),
         "median" = round(probs),
         "mean" = round(probs))
}



### misc. internal functions

## prm: calculate response probabilities for given thetas and betas under the RM.
prm <- function(theta = NULL, beta = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, prm, beta = beta))
  if (is.list(beta)) return(lapply(beta, prm, theta = theta))

  ## calculate probabilities
  return(plogis(outer(theta, beta, "-")))
}

## rrm: calculate response matrices for given thetas and betas under the RM.
rrm <- function(theta = NULL, beta = NULL, return_setting = TRUE)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, rrm, beta = beta, return_setting = return_setting))
  if (is.list(beta)) return(lapply(beta, rrm, theta = theta, return_setting = return_setting))

  ## calculate response probabilities and responses (randomized cutpoint like in eRm:::sim.rasch)
  n <- length(theta)
  m <- length(beta)
  probs <- prm(theta = theta, beta = beta)
  resp <- (matrix(runif(n * m), nrow = n, ncol = m) < probs) + 0

  ## return
  if (return_setting) list(theta = theta, beta = beta, data = resp) else resp
}
