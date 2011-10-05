## workhorse fitting function
RaschModel.fit <- function(y, weights = NULL, start = NULL, gradtol = 1e-6, 
  deriv = c("sum", "diff", "numeric"), hessian = TRUE, ...)
{
  ## argument matching
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
      start <- -qlogis(cs/sum(weights))
      start <- start[-1] - start[1]
    }
    rf <- rf[-1]
    cs <- cs[-1]

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {
      ## obtain esf and apply contrast
      esf <- elementary_symmetric_functions(c(0, par),
        order = 1 - (deriv == "numeric"), 
	diff = deriv == "diff")
      g <- esf[[1]][-1]

      ## conditional log-likelihood
      cll <- sum(-cs * par) - sum(rf * log(g))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      ## gradient
      if(deriv != "numeric") {
        g1 <- esf[[2]][-1, -1, drop = FALSE]
        grad <- - cs + colSums(rf * g1/g)
        attr(cll, "gradient") <- - grad
      }

      return(-cll)
    }

    ## analytical gradient -> now moved to estfun() method
    ## agrad <- function(par, esf)
    ##  weights * (- y + esf[[2]][rs + 1, , drop = FALSE] / esf[[1]][rs + 1])[,-1, drop = FALSE]

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

    ## all NA patterns
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))

    ## starting values
    if(is.null(start)) {
      start <- -qlogis(colSums(y * weights, na.rm = TRUE)/colSums(!y_na * weights))
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

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {

      ## initialize return values    
      cll <- 0
      grad <- rep.int(0, k-1)
    
      ## loop over observed NA patterns
      for(i in levels(na_patterns)) {
        ## parse NA pattern
        wi_i <- as.integer(strsplit(i, "\r")[[1]])
        
        ## select subset
	na_i <- which(na_patterns == i)
        if(length(wi_i) < 1) {
	  y_i <- y[na_i, , drop = FALSE]
	  par_i <- c(0, par)
	} else {
	  y_i <- y[na_i, -wi_i, drop = FALSE]
	  par_i <- c(0, par)[-wi_i]
        }
	weights_i <- weights[na_i]
        cs_i <- colSums(y_i * weights_i)
        rs_i <- rowSums(y_i)
        rf_i <- as.vector(tapply(weights_i, factor(rs_i, levels = 0:ncol(y_i)), sum))
        rf_i[is.na(rf_i)] <- 0

        ## obtain esf and apply contrast
        esf_i <- elementary_symmetric_functions(par_i,
          order = 1 - (deriv == "numeric"), 
  	  diff = deriv == "diff")
        g_i <- esf_i[[1]]

        ## conditional log-likelihood
        cll <- cll + (sum(-cs_i * par_i) - sum(rf_i * log(g_i)))

        ## gradient
        if(deriv != "numeric") {
          g1_i <- esf_i[[2]]
          grad <- grad + zero_fill(-cs_i + colSums(rf_i * g1_i/g_i), wi_i)[-1]
        }
      }

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      ## add gradient (if desired)
      if(deriv != "numeric") {
        attr(cll, "gradient") <- - grad
      }
      
      ## collect and return
      return(-cll)
    }
 
    ## ## analytical gradient -> now moved to estfun() method
    ## agrad <- function(par, esf) {
    ## 
    ##   ## set up return value    
    ##   rval <- matrix(0, nrow = n, ncol = k)
    ## 
    ##   ## loop over observed NA patterns	 
    ##   for(i in seq_along(levels(na_patterns))) {
    ##     ## parse NA pattern
    ##     lev_i <- levels(na_patterns)[i]
    ##     na_i <- which(na_patterns == lev_i)
    ##     wi_i <- as.integer(strsplit(lev_i, "\r")[[1]])
    ##     wi_i <- if(length(wi_i) < 1) 1:k else (1:k)[-wi_i]
    ## 
    ##     ## compute gradient per pattern
    ##     esf_i <- esf[[i]]
    ##     rs_i <- rowSums(y[na_i, wi_i, drop = FALSE])
    ##     rval[na_i, wi_i] <- weights[na_i] * (- y[na_i, wi_i, drop = FALSE] +
    ##       esf_i[[2]][rs_i + 1, , drop = FALSE] / esf_i[[1]][rs_i + 1])
    ##   }
    ## 
    ##   return(rval[, -1, drop = FALSE])
    ## }

    ## analytical hessian
    ahessian <- function(par, esf) {

      ## set up return value    
      rval <- matrix(0, ncol = k-1, nrow = k-1)

      ## loop over observed NA patterns      
      for(i in seq_along(levels(na_patterns))) {
        ## parse NA pattern
	lev_i <- levels(na_patterns)[i]
	na_i <- which(na_patterns == lev_i)
        wi_i <- as.integer(strsplit(lev_i, "\r")[[1]])
        wi2_i <- if(length(wi_i) < 1) 1:k else (1:k)[-wi_i]
	
	## compute hessian per pattern
	rs_i <- rowSums(y[na_i, wi2_i, drop = FALSE])
	k_i <- length(wi2_i)
        rf_i <- as.vector(tapply(weights[na_i], factor(rs_i, levels = 0:k_i), sum))
        rf_i[is.na(rf_i)] <- 0

        ## obtain esf
        g_i <- esf[[i]][[1]]
        g1_i <- esf[[i]][[2]]
        g2_i <- esf[[i]][[3]]

        ## hessian
	hess <- matrix(0, nrow = k_i, ncol = k_i)
        g1s_i <- g1_i/g_i
        for (q in 1:k_i) hess[q,] <- colSums(rf_i * (g2_i[,q,]/g_i - (g1_i[,q]/g_i) * g1s_i))

        rval <- rval + zero_fill(hess, wi_i)[-1, -1]
      }

      return(rval)
    }

  }
  
  ## optimization
  opt <- nlm(cloglik, start, gradtol = gradtol, 
    hessian = (deriv == "numeric") & hessian, check.analyticals = FALSE)
  
  ## collect and annotate results
  cf <- opt$estimate
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
  names(cf) <- rownames(vc) <- colnames(vc) <- colnames(y)[-1]

  ## collect, class, and return
  rval <- list(
    coefficients = cf,
    vcov = vc,
    loglik = -opt$minimum,
    df = k-1,
    data = y_orig,
    weights = if(identical(as.vector(weights_orig), rep(1L, nrow(y_orig)))) NULL else weights_orig,
    n = sum(weights_orig > 0),
    items = status,
    na = any_y_na,
    elementary_symmetric_functions = esf,
    nlm_code = opt$code,
    iterations = opt$iterations,
    gradtol = gradtol
  )
  class(rval) <- "RaschModel"
  return(rval)
}

## elementary symmetric functions
elementary_symmetric_functions <- function(par, order = 0, log = TRUE, diff = FALSE) {
  ## Michelle Liou (1994). More on the Computation of Higher-Order
  ## Derivatives of the Elementary Symmetric Functions in the
  ## Rasch Model. Applied Psychological Measurement, 18(1), 53-62.

  ## Use difference algorithm for orders > 0?
  diff <- rep(diff, length.out = 2)

  ## derivatives up to order
  order <- round(order)[1]
  stopifnot(order %in% 0:2)
  rval <- list()[1:(order+1)]
  names(rval) <- 0:order

  ## transformations
  par <- as.vector(par)
  beta <- if(log) par else -log(par)
  eps <- exp(-beta)
  n <- length(eps)
  stopifnot(n > 1)
  
  ## Order: 0  
  ## initialization: gamma_1(eps_1) = eps_1
  gamma <- c(eps[1], numeric(n-1))

  ## recursion: Equation 3
  for(i in 2:n) gamma[1:i] <- c(eps[i] + gamma[1],
    eps[i] * gamma[(2:i) - 1] + gamma[2:i])
  
  ## gamma_0 = 1
  gamma <- c(1, gamma)

  ## return value
  rval[[1]] <- gamma
  if(order < 1) return(rval)

  ## Order: 1
  if(diff[1]) {
    ## initialization: gamma_1^(j) = 1
    gamma1 <- matrix(0, nrow = n+1, ncol = n)
    gamma1[2,] <- 1
    ## recursion: Equation 4
    for(q in 3:(n+1)) gamma1[q,] <- gamma[q-1] - eps * gamma1[q-1,]
  } else {
    ## re-call self, omitting i-th parameter
    gamma1 <- sapply(1:n, function(i)
      c(0, elementary_symmetric_functions(eps[-i], order = 0, log = FALSE)[[1]]))
  }
  ## if input on log scale: include inner derivative
  if(log) gamma1 <- exp(t(t(log(gamma1)) - beta))
  ## return value
  rval[[2]] <- gamma1
  if(order < 2) return(rval)

  ## Order: 2
  if(diff[2]) {
    ## initialization: gamma_2^(i,j) = 1
    gamma2 <- array(0, c(n+1, n, n))
    gamma2[3,,] <- 1
    ## auxiliary variables
    eps_plus_eps <- outer(eps, eps, "+")
    eps_times_eps <- exp(-outer(beta, beta, "+"))
    ## recursion: Jansen's Equation (Table 1, Forward, Second-Order)
    if(n > 2) for(q in 4:(n+1)) gamma2[q,,] <- gamma[q-2] -
      (eps_plus_eps * gamma2[q-1,,] + eps_times_eps * gamma2[q-2,,])
  } else {
    ## re-call self, omitting i-th and j-th parameter
    gamma2 <- array(0, c(n + 1, n, n))
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        gamma2[, i, j] <- gamma2[, j, i] <- c(0, 0,
	  elementary_symmetric_functions(eps[-c(i, j)], order = 0, log = FALSE)[[1]])
      }
    }
    if(log) eps_times_eps <- exp(-outer(beta, beta, "+"))  
  }
  ## if input on log scale: include inner derivative
  if(log) for(q in 1:(n+1)) gamma2[q,,] <- eps_times_eps * gamma2[q,,]
  ## each diagonal is simply first derivative
  for(q in 2:(n+1)) diag(gamma2[q,,]) <- gamma1[q,]
  ## return value
  rval[[3]] <- gamma2
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
  cat("Number of iterations in nlm optimization:", x$iterations, "\n\n")
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

