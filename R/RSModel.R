### Function to fit a rating scale modell via nlm.  See e.g. Andrich (1978), Andersen (1995)
RSModel.fit <- function (y, weights = NULL, start = NULL, reltol = 1e-10,
                         deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L,
                         full = TRUE, ...)
{
  ## argument matching
  deriv <- match.arg(deriv)

  ## process input and data.
  diff <- deriv == "diff"
  y <- as.matrix(y)
  n_org <- nrow(y)
  m <- ncol(y)
  ident_items <- rep.int(TRUE, m)

  ## process weights
  if (is.null(weights)) weights <- rep.int(1L, n_org) else stopifnot(length(weights) == n_org)
  y <- y[weights > 0, , drop = FALSE]
  weights_org <- weights
  weights <- weights[weights > 0]
  n <- nrow(y)  

  ## move categories to zero if necessary, get number of categories per item
  mincat <- apply(y, 2, min, na.rm = TRUE)
  if(any(mincat > 0)) {
    warning("Minimum score is not zero for all items (", paste("I", which(mincat > 0), sep = "", collapse = ", "), "). These items are scaled to zero.")
    y[, mincat > 0] <- scale.default(y[, mincat > 0], center = mincat[mincat > 0], scale = FALSE)
  }
  oj <- apply(y, 2, max, na.rm = TRUE)

  ## remove unidentified items (only one category used)
  unid <- oj == 0
  nunid <- sum(unid)
  if (nunid) {
    y <- y[, !unid, drop = FALSE]
    oj <- oj[!unid]
    m <- m - nunid
    ident_items[unid] <- FALSE
    warning("There were unidentified items (only one category used, ", paste("I", which(unid), sep = "", collapse = ", "), ").")
  }

  ## stop if not all oj equal or oj = 1 (dichotomous items).
  o <- max(oj)
  if (any(oj != o)) warning("Not all items have the same number of categories. Consider using PCModel.fit() for data with a variable number of categories per item.")
  oj <- rep.int(o, m)
  if (o == 1) stop("Rating scale model requires polytomous items. Use RaschModel.fit() for dichotomous items.")

  ## check for 'total' null categories (null for all items), if present, stop, because tau is assessable
  cat_freq <- apply(y + 1, 2, tabulate, nbins = o + 1)
  if (any(rowSums(cat_freq) == 0)) stop("There are categories which are null for all items. ",
                   "Please use PCModel.fit() instead (see argument nullcats in ?PCModel.fit).")

  ## check for missing
  y_na <- is.na(y)
  any_na <- any(y_na)

  ## set starting values
  if (is.null(start)) start <- rep.int(0, m + o - 2)

  ## small helper function: conversion of a given rating scale parameter vector (and additional information)
  ## to a pcm parameter vector (item-category-parameters) (expects that beta_0 = 0 and tau_0 = 0.)
  rsm2pcm <- function(rsm_par, item_par_pos, cat_par_pos, oj, cat_scores) rep.int(c(0, rsm_par[item_par_pos]), oj) * cat_scores + c(0, rsm_par[cat_par_pos])

  ## calculate positions and helper variables: ipar/cpar without and mv/ov with restricted parameters
  ipar <- 1:(m-1)
  cpar <- m:(m + o - 2)
  mv <- 1:m
  ov <- 1:o

  ## calculate item, category and person totals
  itot <- colSums(y * weights)[-1]                                   # without item 1
  cs <- t(apply(y, 1, tabulate, nbins = o))
  ctot <- colSums(cs * weights)[-1] # without category 0 and 1
  rs <- rowSums(y, na.rm = TRUE)
  ptot <- as.vector(tapply(weights, factor(rs, levels = 0:sum(oj)), sum))
  ptot[is.na(ptot)] <- 0
  
  ## build log-likelihood
  if(!any_na) {

    ## objective function: conditional log-likelihood
    cloglik <- function (par) {

      ## transform rsm par to pcm par list for esf function
      esf_par <- split(rsm2pcm(par, ipar, cpar, oj, ov), rep.int(mv, oj))

      ## finally: the cll
      cll <- - sum(itot * par[ipar]) - sum(ctot * par[cpar]) - sum(ptot * log(elementary_symmetric_functions(par = esf_par, order = 0, diff = diff)[[1]]))

      ## catch degenerated cases (typically caused by non-finite gamma)
      if (is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax
      
      return(-cll)
    }
    
    ## static helper variables for analytical gradient
    betaindex <- rep.int(mv,  rep.int(o, m))
    tauindex <- rep.int(ov, m)

    ## analytical gradient
    agrad <- function (par) {
        
        ## transform RSM parameters to PCM parameters and calculate actual ESF
        parx <- rsm2pcm(par, ipar, cpar, oj, ov)
        esf <- elementary_symmetric_functions(par = split(parx, rep.int(mv, oj)), order = 1, diff = diff)
        gamma0 <- esf[[1]][rs + 1]
        gamma1_pcm <- esf[[2]]

        ## calculate transformed derivatives, select relevant derivatives with ptot, drop unindentified parameters.
        gamma1 <- matrix(0, nrow = nrow(gamma1_pcm), ncol = m + o)
        for (j in mv) gamma1[, j] <- gamma1_pcm[, betaindex == j, drop = FALSE] %*% ov
        for (k in ov) gamma1[, k + m] <- rowSums(gamma1_pcm[, tauindex == k, drop = FALSE])
        gamma1 <- apply(gamma1, 2, "[", rs + 1)

        ## finally: the gradient
        agrad <- matrix(0, nrow = n, ncol = m + o)
        agrad[, mv] <- weights * (- y + (gamma1 / gamma0)[, mv, drop = FALSE])
        agrad[, m + ov] <- weights * (- cs + (gamma1 / gamma0)[, m + ov, drop = FALSE])

        ## return aggegrated gradient
        return(- colSums(agrad[, -c(1, m + 1), drop = FALSE]))
    }

  } else {

    ## fetch different NA-patterns like Achim does
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    lev_na_patterns <- levels(na_patterns)
    par_pos <- rep.int(mv, oj) ## calculate positions of pcm parameters (needed for ESF)

    ## setup na_pattern lists for various elements of loglik
    ptot_i <- itot_i <- ctot_i <- ipar_i <- oj_i <- pcm_par_i <- m_i <- mv_i <- vector("list", length(lev_na_patterns))
    na_obs_i <- rs_i <- weights_i <- vector("list", length(lev_na_patterns))

    ## loop over observed NA patterns, calculate fix things once and store in list
    for(i in seq_along(lev_na_patterns)) {
      
      ## from pattern i: get NA item(s), observations of and weights with this pattern
      na_items_i <- as.integer(strsplit(lev_na_patterns[i], "\r")[[1]])
      n_na_items_i <- length(na_items_i)
      na_obs_i[[i]] <- which(na_patterns == lev_na_patterns[i])
      weights_i[[i]] <- weights[na_obs_i[[i]]]
      
      ## select subset
      if(n_na_items_i < 1) {            # no missings
        y_i <- y[na_obs_i[[i]], , drop = FALSE]
        m_i[[i]] <- m
        mv_i[[i]] <- mv
        pcm_par_i[[i]] <- rep.int(TRUE, m*o)
        ipar_i[[i]] <- ipar
        oj_i[[i]] <- oj
      } else {                          # specific NA pattern
        y_i <- y[na_obs_i[[i]], -na_items_i, drop = FALSE]
        rs_i[[i]] <- rowSums(y_i)
        m_i[[i]] <- m - n_na_items_i
        mv_i[[i]] <- mv[-na_items_i]
        pcm_par_i[[i]] <- !(par_pos%in% na_items_i)
        ipar_i[[i]] <- if ((n_na_items_i == 1) && na_items_i == 1) ipar else ipar[-(setdiff(na_items_i, 1) - 1)] # beta_j is in col j-1
        oj_i[[i]] <- oj[-na_items_i]
    }

      ## calculate category and person totals for NA-group i
      itot_i[[i]] <- colSums(y_i * weights_i[[i]])
      if (!(1 %in% na_items_i)) itot_i[[i]] <- itot_i[[i]][-1]                           # remove item 1 if not already removed
      ctot_i[[i]] <- colSums(t(apply(y_i, 1, tabulate, nbins = o)) * weights_i[[i]])[-1] # without category 0 and 1
      rs_i[[i]] <- rowSums(y_i)
      ptot_i[[i]] <- as.vector(tapply(weights_i[[i]], factor(rs_i[[i]], levels = 0:sum(oj_i[[i]])), sum))
      ptot_i[[i]][is.na(ptot_i[[i]])] <- 0
    }

    ## fun for mapply (calculates cll contributions per na_group with variable par's)
    cll_i <- function (itot_i, item_par_i, ctot_i, ptot_i, esf_par_i, cat_par) {
      - sum(itot_i * item_par_i) - sum(ctot_i * cat_par) - sum(ptot_i * log(elementary_symmetric_functions(par = esf_par_i, order = 0, diff = diff)[[1]]))
    }

    ## objective function: conditional log-likelihood (build incrementally by looping over the NA-patterns
    cloglik <- function (par) {

      ## fetch parameter vectors
      esf_par <- rsm2pcm(par, ipar, cpar, oj, ov)
      esf_par_i <- lapply(pcm_par_i, function(x) esf_par[x])
      esf_par_i <- mapply(split, esf_par_i, mapply(function (x, y) rep.int(1:x, y), m_i, oj_i))
      item_par_i <- lapply(ipar_i, function(x) par[x])
      cat_par <- par[cpar]

      ## conditional log-likelihood
      cll <- sum(mapply(cll_i, itot_i, item_par_i, ctot_i, ptot_i, esf_par_i,
                        MoreArgs = list(cat_par = cat_par), SIMPLIFY = TRUE, USE.NAMES = FALSE))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      return(-cll)
    }

    ## static helper variables for analytical gradient
    betaindex_i <- mapply(function (x, y, ...) rep.int(x, rep.int(o, y)), mv_i, m_i, MoreArgs = list(o = o), SIMPLIFY = FALSE)
    tauindex_i <- lapply(m_i, function (x) rep.int(ov, x))
    
    ## analytical gradient
    agrad <- function (par) {

        ## transform RSM parameters to PCM parameters
        parx <- rsm2pcm(par, ipar, cpar, oj, ov)
        esf_par_i <- lapply(pcm_par_i, function (x) parx[x])
        esf_par_i <- mapply(split, esf_par_i, mapply(function (x, y) rep.int(1:x, y), m_i, oj_i))
        esf_i <- mapply(elementary_symmetric_functions, par = esf_par_i, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)

        ## loop over obseved NA patterns and gradually built up analytical gradient
        grad <- matrix(0, nrow = n, ncol = m + o)
        for (i in seq_along(lev_na_patterns)) {

            ## fetch ESF with PCM parametrization
            gamma0_i <- esf_i[[i]][[1]][rs_i[[i]] + 1]
            gamma1_pcm_i <- esf_i[[i]][[2]]

            ## transform derivatives via delta rule to RSM parametrization
            gamma1_i <- matrix(0, nrow = nrow(gamma1_pcm_i), ncol = m_i[[i]] + o)
            for (j in 1:m_i[[i]]) gamma1_i[, j] <- gamma1_pcm_i[, betaindex_i[[i]] == mv_i[[i]][j], drop = FALSE] %*% ov
            for (k in ov) gamma1_i[, k + m_i[[i]]] <- rowSums(gamma1_pcm_i[, tauindex_i[[i]] == k, drop = FALSE])
            gamma1_i<- apply(gamma1_i, 2, "[", rs_i[[i]] + 1)
            if (!is.matrix(gamma1_i)) gamma1_i<- matrix(gamma1_i, nrow = 1)

            ## finally: the gradient for NA group i
            grad[na_obs_i[[i]], mv_i[[i]]] <- weights_i[[i]] * (- y[na_obs_i[[i]], mv_i[[i]], drop = FALSE] + (gamma1_i/gamma0_i)[, 1:m_i[[i]], drop = FALSE])
            grad[na_obs_i[[i]], m + ov] <- weights_i[[i]] * (- cs[na_obs_i[[i]], , drop = FALSE] + (gamma1_i/gamma0_i)[, m_i[[i]] + ov, drop = FALSE])

        }

        return(- colSums(grad[, -c(1, m + 1), drop = FALSE]))
    }
  }

  ## optimization
  opt <- optim(par = start, fn = cloglik, gr = agrad, method = "BFGS",
               hessian = hessian, control = list(reltol = reltol, maxit = maxit, ...))
 
  ## final estimates ...
  est <- opt$par
  names(est) <- if (is.null(colnames(y))) {
    c(paste("I", setdiff(which(ident_items), 1), sep = ""), paste("C", 2:o, sep = ""))
  } else c(colnames(y)[setdiff(which(ident_items), 1)], paste("C", 2:o, sep = ""))
    

  ## ... and (if requested) esf of these ...
  if (full) {
    parx <- rsm2pcm(est, ipar, cpar, oj, ov)
    esf <- if (any_na) {
      esf_par_i <- lapply(pcm_par_i, function (x) parx[x])
      esf_par_i <- mapply(split, esf_par_i, mapply(function (x, y) rep.int(1:x, y), m_i, oj_i))
      mapply(elementary_symmetric_functions, par = esf_par_i, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)
    } else {
      elementary_symmetric_functions(par = split(parx, rep.int(mv, oj)), order = 1, diff = diff)
    }
    
    ## ... as well as variance-bcovariance matrix
    if (hessian) {
      vc <- opt$hessian
      vc <- solve(vc)
    } else {
      vc <- matrix(NA, nrow = length(est), ncol = length(est))
    }
    rownames(vc) <- colnames(vc) <- names(est)
  } else {
    esf <- NULL
    vc <- NULL
  }


  ## collect results, set class, and return
  res <- list(coefficients = est,
              vcov = vc,
              data = y,
              items = ident_items,
              categories = oj,
              n = sum(weights_org > 0),
              n_org = n_org,
              weights = if (identical(as.vector(weights_org), rep.int(1L, n_org))) NULL else weights_org,
              na = any_na,
              esf = esf,
              loglik = -opt$value,
              df = length(est),
              code = opt$convergence,
              iterations = tail(na.omit(opt$counts), 1L),
              reltol = reltol)
  class(res) <- c("RSModel", "PCModel")
  return(res )
}

print.RSModel <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("RSM item and threshold parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

coef.RSModel <- function (object, ...) object$coefficients

vcov.RSModel <- function (object, ...) object$vcov

logLik.RSModel <- function (object, ...) structure(object$loglik, df = object$df, class = "logLik")

weights.RSModel <- function (object, ...) if (is.null(object$weights)) rep.int(1L, object$n_org) else object$weights

summary.RSModel <- function (object, vcov. = NULL, ...) {
  ## coefficients
  cf <- coef(object)

  ## covariance matrix
  if (is.null(vcov.)) 
    vc <- vcov(object)
  else {
    if (is.function(vcov.)) vc <- vcov.(object)
    else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2 * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf      
  class(object) <- "summary.RSModel"
  return(object)
}


print.summary.RSModel <- function (x, digits = max(3, getOption("digits") - 3),
                                   signif.stars = getOption("show.signif.stars"), ...)
{
  if (is.null(x$call)) {
    cat("\nRating scale model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if (any(!x$items)) cat("Excluded items:",
                         paste(names(x$items)[!x$items], collapse = ", "), "\n\n")

  cat("Item location and threshold parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
      "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in BFGS optimization:", x$iterations, "\n\n")
  invisible(x)
}

plot.RSModel <- function (x, pattern = TRUE, names = NULL, ref = NULL, main = NULL,
                          ylab = "Latent trait", ylim = NULL, col = c("gray", "gray"), pch = c(21, 22),
                          lty = c(2, 3), refline = TRUE, reflinecol = "lightgray", ...)
{
  ## pattern = FALSE -> goto plot.PCModel()
  if (!pattern) plot.PCModel(x, ...)
  else {

    ## check input
    stopifnot(length(col) == 2)
    stopifnot(length(pch) == 2)
    stopifnot(length(lty) == 2)
    
    ## setup basic data
    m <- sum(x$items)
    mvec <- which(x$items)
    o <- x$categories[1]

    ## setup parameters
    cf <- itempar.RSModel(x, ref = ref, vcov = FALSE, simplify = FALSE)
    ipar <- cf[[1]]
    cpar <- cf[[2]]

    ## setup y-axis and labels
    if (is.null(ylim)) ylim <- extendrange(c(ipar, cpar), f = 0.25)
    if (is.null(names)) names <- c(paste("Item", mvec), paste("Tau", 1:o)) else stopifnot(length(names) == m + o)

    ## setup plotting parameters and par
    par(mar = c(2.5, 4.5, 2.5, 1.5))
    ix <- 1:(m + o)

    ## setup plotting window and (if requested) reference line
    plot(ix, unlist(cf, use.names = FALSE), main = main, ylab = ylab, type = "n", axes = FALSE, ylim = ylim)
    if (refline) lines(x = c(0, ix[m] + 0.49), y = rep(mean(ipar), 2), col = reflinecol, ...)
    axis(2)
    axis(1, at = ix, labels = names)
    box()

    ## plot item parameters
    lines(x = ix[1:m], y = ipar, type = "l", lty = lty[1], ...)
    points(x = ix[1:m], y = ipar, pch = pch[1], bg = col[1], ...)

    ## plot seperation line and threshold parameters
    abline(v = ix[m] + c(0.49, 0.51))
    lines(x = ix[(m + 1):(m + o)], y = cpar, type = "l", lty = lty[2], ...)
    points(x = ix[(m + 1):(m + o)], y = cpar, pch = pch[2], bg = col[2], ...)
  }
}

itempar.RSModel <- function (object, ref = NULL, vcov = TRUE, simplify = TRUE, ...) {

  ## extract cf and labels, include restricted parameters
  cf <- coef(object)
  k <- length(cf)
  m <- sum(object$items)
  icf <- c(0.00, cf[1:(m-1)])
  ccf <- c("C1" = 0.00, cf[m:k])
  ilbs <- names(icf)
  clbs <- names(ccf)
  names(cf)[1] <- names(icf)[1] <- ilbs[1] <- if (ilbs[2] == "I2") "I1" else colnames(object$data)[1]
  
  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.character(ref)) {
    stopifnot(all(ref %in% ilbs))
    ref <- which(ilbs %in% ref)
  } else if (is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else stop("Argument 'ref' can only be a character vector with item labels or a numeric vector with item indices.")

  ## create contrast matrix
  D <- diag(m)
  D[, ref] <- D[, ref] - 1/length(ref)
  
  ## impose restriction
  icf <- D %*% icf
  names(icf) <- ilbs

  ## adjust vcov if requested
  if (vcov) {
    k <- k + 2
    vc <- matrix(0.00, nrow = k, ncol = k)
    vc[c(2:m, (m+2):k), c(2:m, (m+2):k)] <- vcov(object)

    ## create new D to include covariances of ip and taus in transformation
    D2 <- diag(k)                       
    D2[1:m, 1:m] <- D
    vc <- D2 %*% vc %*% t(D2)
    colnames(vc) <- rownames(vc) <- c(ilbs, clbs)
  }

  ## create list if requested
  cf <- if (simplify) c(icf, ccf) else list(icf, ccf)

  ## return results
  rval <- structure(cf, class = "itempar", model = "RSM", ref = ref, vcov = if (vcov) vc else NULL)
  return(rval)
}

threshold.RSModel <- function (object, type = c("mode", "median", "mean", "unmodified"), ref = NULL, simplify = TRUE, ...) {

  ## check type, extract item parameters, add names, backup attributes
  type <- match.arg(type)
  ip <- itempar.RSModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  lbs <- lapply(names(ip[[1]]), paste, names(ip[[2]]), sep = "-")
  mdl <- attr(ip, "model")
  ref <- attr(ip, "ref")
  ip <- lapply(as.list(ip[[1]]), function (beta) diff(0:length(ip[[2]]) * beta + c(0, ip[[2]]))) ## convert to pcm parameters, rest as before
  ip <- mapply("names<-", ip, lbs, SIMPLIFY = FALSE)
  names(ip) <- NULL

  ## calculate threshold parameters (if necessary)
  if (type == "mode") {
   ## check if all threshold parameters are in order, if not, calculate sorted ones
    us <- sapply(ip, is.unsorted)
    if (any(us)) {
      usj <- which(us)
      for (j in usj) {
        ipj <- ip[[j]]
        nj <- length(ipj)
        
        ## check if there is a point with a biggest parameter, if yes, take mean
        for (i in 1:nj) {
          if (all(ipj[i] > ipj[(i+1):nj])) {
            ipj[i] <- mean(ipj[i:nj])
            ipj <- ipj[-(i+1:nj)]
            break
          }
        }
        
        ## recursive sorting if there is still unorder (e.g. 4, 2, 3, 1)
        while(is.unsorted(ipj)) {
          uo_pos <- which(diff(ipj) < 0)                             # locate unordered parameters, returns position of the first
          ipj[uo_pos] <- (ipj[uo_pos] + ipj[uo_pos + 1]) / 2 # replace first with location of intersection of ccc curves (= (eps1 + eps2)/ 2)
          ipj <- ipj[-(uo_pos + 1)]                              # remove second
        }

        ip[[j]] <- ipj
      }
    }
  } else {
    oj <- sapply(ip, length)

    if (type == "median") {
      ## function to find locations on theta axis
      zmedian <- function (theta = NULL, delta = NULL, geq = NULL, ncat = NULL) {
        rowSums(ppcm(theta = theta, delta = delta)[, (geq + 1):ncat, drop = FALSE]) - 0.5
      }
      
      ## loop though items and find locations by means of zmedian() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(1:oj[j], function (geq) uniroot(f = zmedian, interval = c(-10, 10), delta = ip[[j]], geq = geq, ncat = oj[j] + 1)$root)
      }
    } else if (type == "mean") {
      ## function to find locations on theta axis
      xpct <- lapply(oj, function (oj) 1:oj - 0.5)
      zexpct <- function (theta = NULL, delta = NULL, expct = NULL) ppcm(theta = theta, delta = delta) %*% 0:length(delta) - expct
      
      ## loop though items and find locations by means of zexpct() and uniroot()
      for (j in seq_along(ip)) {
        ip[[j]] <- sapply(xpct[[j]], function (xp) uniroot(f = zexpct, interval = c(-10, 10), delta = ip[[j]], expct = xp)$root)
      }
    }

    ## add labels again
    ip <- mapply("names<-", ip, lbs, SIMPLIFY = FALSE)
  }
  
  ## backup attributes, simplify if requested, then set attributes again
  if (simplify) ip <- unlist(ip)
  return(structure(ip, "class" = "threshold", "model" = mdl, "ref" = ref, "type" = type))
}

personpar.RSModel <- function (object, ref = NULL, vcov = TRUE, start = NULL, tol = 1e-6, ...) {

  ## extract item parameters and relevant informations #FIXME: constrained parameter included?
  ipar <- itempar.RSModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  m <- length(ipar[[1]])
  o <- length(ipar[[2]])
  os <- 1:o
  osv <- rep.int(os, m)
  ipar <- lapply(as.list(ipar[[1]]), function (beta) os * beta + ipar[[2]]) ## convert to pcm parameters
  rng <- 1:(m * o - 1)
  
  ## if vcov is requested, fetch relevant data for loglik
  if (vcov) {
    y <- object$data[weights(object) > 0, , drop = FALSE]
    rs <- rowSums(y)
    rf <- tabulate(rs, nbins = m * o - 1)
    cs <- apply(y, 2, tabulate, nbins = o)

    ## remove unidentified parameters
    rs <- rs[rf != 0]
    rng <- rng[rf != 0]
    rf <- rf[rf != 0]

    ## transform ipar to matrix (o x m)
    ipar <- do.call("cbind", ipar)
  }

  ## start values / range of person parameters
  if(is.null(start)) start <- if (vcov) numeric(length(rng)) else c(-1, 1) * qlogis(1/m * 1e-3) #FIXME: 1e3 enough?

  if (vcov) {
    ## objective function
    cloglik <- function (ppar) {
      pparx <- outer(os, ppar) # l * theta_i
      - sum(rf * rng * ppar) + sum(cs * ipar) + sum(rf * colSums(log(1 + apply(pparx, 2, function (x) colSums(exp(x - ipar))))))
    }
    ## optimization & target variables
    opt <- nlm(cloglik, start, gradtol = tol, hessian = TRUE, check.analyticals = FALSE, ...)
    ppar <- opt$estimate
    vc <- solve(opt$hessian)
  } else {
    ## iterate over raw scores (from 1 until rmax-1)
    ppar <- sapply(1:(m * o - 1), function(rawscore) {
      uniroot(function(ppar) rawscore - sum(osv * exp(osv * ppar - unlist(ipar)) / unlist(lapply(ipar, function (beta) rep.int(1 + sum(exp(os * ppar - beta)), o)))), interval = start, tol = tol)$root
    })
  }
  
  ## return named vector of thetas
  return(structure(ppar, .Names = rng, class = "personpar", model = "RSM"))
}

predict.RSModel <- function (object, newdata = NULL, type = c("probability", "cumprobability", "mode", "median", "mean"),
                             ref = NULL, ...) {
  
  ## check type, process newdata, if NULL, use person parameters of given model object
  type <- match.arg(type)
  if (is.null(newdata)) {
    rs <- rowSums(object$data, na.rm = TRUE)
    rs <- rs[0 < rs & rs < ncol(object$data)]
    newdata <- personpar.RSModel(object, ref = ref, vcov = FALSE)[rs]
    names(newdata) <- NULL
  }
  nms <- names(newdata)

  ## itempars and raw probabilities
  ip <- itempar.RSModel(object, ref = ref, vcov = FALSE, simplify = FALSE)
  ilbs <- names(ip[[1]])
  os <- 0:length(ip[[2]])
  probs <- prsm(theta = newdata, beta = ip[[1]], tau = ip[[2]])

  ## return as requested in type, for RM mode, median, mean is the same
  switch(type,
         "probability" = {
           ## construct labels
           lbs <- as.vector(t(outer(ilbs, os, paste, sep = "-C")))
           ## create (named) matrix with probabilities
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, lbs)
           probs
           },
         "cumprobability" = {
           ## construct labels
           lbs <- as.vector(t(outer(ilbs, os, paste, sep = ">=C")))
           ## create (named) matrix with probabilities
           probs <- lapply(probs, function (probsj) t(apply(probsj, 1, function (x) rev(cumsum(rev(x))))))
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, lbs)
           probs
           },
         "mode" = sapply(probs, apply, 1, which.max) - 1,
         "median" = sapply(probs, function (probsj) apply(probsj, 1, which.max) - 1),
         "mean" = round(sapply(probs, function (probsj) apply(probsj, 1, function (j) sum(j * 0:(length(j) - 1))))),round(probs))
}



### misc. internal functions

## prsm: calculate response probabilities for given thetas, betas and taus under the RSM.
prsm <- function(theta = NULL, beta = NULL, tau = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta) && !is.null(tau))

  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, prsm, beta = beta, tau = tau))
  if (length(beta) > 1) return(lapply(as.list(beta), prsm, theta = theta, tau = tau))

  ## calculate probabilities
  num <- cbind(0, outer(theta - beta, 1:length(tau), "*"))         # k * (theta_i - beta_j ) for all i, j and k = 1:o, add 0 for cat 0
  num <- exp(t(apply(num, 1, function (x) x - cumsum(c(0, tau))))) # - sum(tau_k)
  denom <- rowSums(num)                                            # denominator: sum over all numerators
  return(num/denom)
}

## rpcm: calculate response matrices for given thetas, betas and taus under the RSM.
rrsm <- function(theta = NULL, beta = NULL, tau = NULL, nullcats = FALSE, return_setting = TRUE)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta) && !is.null(tau))
  if (is.list(beta)) stopifnot(is.list(tau) && (length(beta) == length(tau)))

  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, rrsm, beta = beta, tau = tau))
  if (is.list(beta)) return(mapply(rrsm, beta, tau, MoreArgs =
                                   list(theta = theta, nullcats = nullcats, return_setting = return_setting), SIMPLIFY = FALSE))
  
  ## calculate response probabilities
  probs <- prsm(theta = theta, beta = beta, tau = tau)
  if(!is.list(probs)) probs <- list(probs)

  ## calculate response matrices for given set of item
  rsp_item_i <- function (probmat_i) {  #inline function to calculate responses per item
    oj <- ncol(probmat_i)               #calculate rsp vector once
    rsp <- apply(probmat_i, 1, function (p) sample.int(n = oj, size = 1, prob = p))
    if(!nullcats) {                     #recalculate if null categories not allowed
      while (length(unique.default(rsp)) != oj) {
        rsp <- apply(probmat_i, 1, function (p) sample.int(n = oj, size = 1, prob = p))
      }
    }
    return(rsp-1)
  }
  res <- lapply(probs, rsp_item_i)
  res <- do.call(cbind, res)
  
  if (return_setting)
    return(list(theta = theta, beta = beta, tau = tau, data = res))
  else
    return(res)
}
